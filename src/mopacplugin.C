/*
Mopac VMD plugin
Authors: Maximilian Scheurer, Marcelo Melo, May 2017
 */

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <cctype>
#include "qmplugin.h"
#include "unit_conversion.h"
#include "periodic_table.h"

typedef std::vector<std::vector<std::vector<float>>> MoCoeff;
typedef std::vector<std::vector<float>> CoeffRowBlock;

#define DEBUGGING 1
#ifdef DEBUGGING
#define PRINTERR fprintf(stderr, "\n In file %s, line %d: \n %s \n \n", \
                           __FILE__, __LINE__, strerror(errno))
#else
#define PRINTERR (void)(0)
#endif

#define ANGSTROM 0
#define BOHR     1

/*
 * Error reporting macro for the multiple fgets calls in
 * the code
 */
#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE

#define NOTFOUND 0
#define FOUND    1
#define STOPPED  2


typedef struct {
  int version;
} mopacdata;

// Checks if loaded file is really an Mopac output file
static int have_mopac(qmdata_t *data, mopacdata* mopac);
static int parse_static_data(qmdata_t *data, int* natoms);
static void* open_mopac_read(const char* filename, const char* filetype, int *natoms);


static int have_mopac(qmdata_t *data, mopacdata* mopac) {
  int programLine;
  int versionLine;
  char buffer[BUFSIZ];
  int mainVersion, secondDigit, thirdDigit;
  buffer[0] = '\0';
  programLine = goto_keyline(data->file, "**                                MOPAC2016                                  **", NULL);
  if (programLine != 1) {
    return FALSE;
  } else {
    mopac->version = 2016;
  }
  return TRUE;
}

static int parse_static_data(qmdata_t *data, int* natoms) {
  mopacdata *mopac = (mopacdata *)data->format_specific_data;
  // if (!get_job_info(data)) return FALSE;

  // if (!get_input_structure(data, mopac)) return FALSE;

  // if (!get_basis(data)) return FALSE;

  // if (!analyze_traj(data, mopac)) {
    // printf("mopacplugin) WARNING: Truncated or abnormally terminated file!\n\n");
  // }

  // *natoms = data->numatoms;

  // read_first_frame(data);

  // print_input_data(data);

  return FALSE;
}



static void* open_mopac_read(const char* filename, const char* filetype, int *natoms) {
  FILE* fd;
  qmdata_t *data = NULL;

  #ifdef DEBUGGING
    printf("DEBUG: Open Mopac Read called: %s\n", filename);
  #endif

  mopacdata* mopac;

  // open the mopac output files
  fd = fopen(filename, "rb");
  if (!fd) {
    PRINTERR;
    return NULL;
  }

  // initialize/allocate qm data (qmplugin.h)
  data = init_qmdata();
  if (data == NULL) {
    PRINTERR;
    return NULL;
  }

  // mopac specific information alloc
  mopac = (mopacdata *) calloc(1, sizeof(mopacdata));
  mopac->version = 0;
  data->format_specific_data = mopac;

  // filename in qm data struct
  data->file = fd;

  if(have_mopac(data, mopac)) {
    std::cout << mopac->version << std::endl;
    if (mopac->version == 2016) {
      printf("mopacplugin) Mopac version: %d\n", mopac->version);
    } else {
      printf("mopacplugin) Mopac version not supported: %d\n", mopac->version);
      return NULL;
    }
    if (parse_static_data(data, natoms) == FALSE) {
      return NULL;
    }
  } else {
    printf("mopacplugin) This is not an Mopac output file!\n");
    return NULL;
  }

  return data;
}

static int read_mopac_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  qmdata_t *data = (qmdata_t *)mydata;
  qm_atom_t *cur_atom;
  molfile_atom_t *atom;
  int i = 0;
  *optflags = MOLFILE_ATOMICNUMBER;

  cur_atom = data->atoms;

  for(i=0; i<data->numatoms; i++) {
    atom = atoms+i;
    strncpy(atom->name, cur_atom->type, sizeof(atom->name));
    strncpy(atom->type, cur_atom->type, sizeof(atom->type));
    strncpy(atom->resname,"", sizeof(atom->resname));
    atom->resid = 1;
    atom->chain[0] = '\0';
    atom->segid[0] = '\0';
    atom->atomicnumber = cur_atom->atomicnum;
    #ifdef DEBUGGING
    printf("mopacplugin) atomicnum[%d] = %d\n", i, atom->atomicnumber);
    #endif

    /* if (data->have_mulliken)
    atom->charge = data->qm_timestep->mulliken_charges[i];
    */
    cur_atom++;
  }

  return MOLFILE_SUCCESS;
}



/*************************************************************
 *
 * plugin registration
 *
 **************************************************************/
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "mopac";
  plugin.prettyname = "Mopac";
  plugin.author = "Maximilian Scheurer";
  plugin.majorv = 0;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "mopac";
  plugin.open_file_read = open_mopac_read;
  // plugin.read_structure = read_mopac_structure;
  // plugin.close_file_read = close_mopac_read;
  //
  // plugin.read_qm_metadata = read_mopac_metadata;
  // plugin.read_qm_rundata  = read_mopac_rundata;

#if vmdplugin_ABIVERSION > 11
  // plugin.read_timestep_metadata    = read_timestep_metadata;
  // plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;
  // plugin.read_timestep = read_timestep;
#endif

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}



/**********************************************************
 *
 * close file and free memory
 *
 **********************************************************/
 // move to top later...
static void close_mopac_read(void *mydata) {
  printf("Freeing memory.\n");

  qmdata_t *data = (qmdata_t *)mydata;
  int i, j;
  fclose(data->file);

  free(data->atoms);
  free(data->basis);
  free(data->shell_types);
  free(data->atomicnum_per_basisatom);
  free(data->num_shells_per_atom);
  free(data->num_prim_per_shell);
  free(data->bonds);
  free(data->angles);
  free(data->dihedrals);
  free(data->impropers);
  free(data->internal_coordinates);
  free(data->bond_force_const);
  free(data->angle_force_const);
  free(data->dihedral_force_const);
  free(data->improper_force_const);
  free(data->inthessian);
  free(data->carthessian);
  free(data->wavenumbers);
  free(data->intensities);
  free(data->normal_modes);
  free(data->imag_modes);
  free(data->angular_momentum);
  data->angular_momentum = NULL;
  free(data->filepos_array);

  if (data->basis_set) {
    for(i=0; i<data->num_basis_atoms; i++) {
      // printf("Freeing basis set of atom %d\n", i);
      for (j=0; j<data->basis_set[i].numshells; j++) {
        // printf("Freeing shell primitives %d\n", j);
        // printf("--- Address of prim is %p\n", (void *)data->basis_set[i].shell[j].prim);
        free(data->basis_set[i].shell[j].prim);
	      data->basis_set[i].shell[j].prim = NULL;
      }
      // printf("- Address of shell is %p\n", (void *)data->basis_set[i].shell);
      free(data->basis_set[i].shell);
      data->basis_set[i].shell = NULL;
      // printf("- Address of shell is %p\n", (void *)data->basis_set[i].shell);
    }
    free(data->basis_set);
    data->basis_set = NULL;
  }

  for (i=0; i<data->num_frames; i++) {
    free(data->qm_timestep[i].scfenergies);
    free(data->qm_timestep[i].gradient);
    free(data->qm_timestep[i].mulliken_charges);
    free(data->qm_timestep[i].lowdin_charges);
    free(data->qm_timestep[i].esp_charges);
    for (j=0; j<data->qm_timestep[i].numwave; j++) {
      free(data->qm_timestep[i].wave[j].wave_coeffs);
      free(data->qm_timestep[i].wave[j].orb_energies);
      free(data->qm_timestep[i].wave[j].orb_occupancies);
    }
    free(data->qm_timestep[i].wave);
  }
  free(data->qm_timestep);
  free(data->format_specific_data);
  free(data);
}



static void print_input_data(qmdata_t *data) {
  int i, j, k;
  int primcount=0;
  int shellcount=0;

  printf("\nDATA READ FROM FILE:\n\n");
  printf(" %10s WORDS OF MEMORY AVAILABLE\n", data->memory);
  printf("\n");
  printf("     BASIS OPTIONS\n");
  printf("     -------------\n");
  printf("%s\n", data->basis_string);
  printf("\n\n\n");
  printf("     RUN TITLE\n");
  printf("     ---------\n");
  printf(" %s\n", data->runtitle);
  printf("\n");
  printf(" THE POINT GROUP OF THE MOLECULE IS %s\n", "XXX");
  printf(" THE ORDER OF THE PRINCIPAL AXIS IS %5i\n", 0);
  printf("\n");
  printf(" YOUR FULLY SUBSTITUTED Z-MATRIX IS\n");
  printf("\n");
  printf(" THE MOMENTS OF INERTIA ARE (AMU-ANGSTROM**2)\n");
  printf(" IXX=%10.3f   IYY=%10.3f   IZZ=%10.3f\n", 0.0, 0.0, 0.0);
  printf("\n");
  printf(" ATOM      ATOMIC                      COORDINATES (BOHR)\n");
  printf("           CHARGE         X                   Y                   Z\n");
  for (i=0; i<data->numatoms; i++) {
    printf(" %-8s %6d", data->atoms[i].type, data->atoms[i].atomicnum);

    printf("%17.10f",   ANGS_TO_BOHR*data->atoms[i].x);
    printf("%20.10f",   ANGS_TO_BOHR*data->atoms[i].y);
    printf("%20.10f\n", ANGS_TO_BOHR*data->atoms[i].z);
  }
  printf("\n");
  printf("     ATOMIC BASIS SET\n");
  printf("     ----------------\n");
  printf(" THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED\n");
  printf(" THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY\n");
  printf("\n");
  printf("  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)\n");
  printf("\n");

#if 0
  for (i=0; i<data->numatoms; i++) {
    printf("%-8s\n\n", data->atoms[i].type);
    printf("\n");
    printf("nshells=%d\n", data->num_shells_per_atom[i]);

    for (j=0; j<data->num_shells_per_atom[i]; j++) {
      printf("nprim=%d\n", data->num_prim_per_shell[shellcount]);

      for (k=0; k<data->num_prim_per_shell[shellcount]; k++) {
        printf("%6d   %d %7d %22f%22f\n", j, data->shell_types[shellcount],
               primcount+1, data->basis[2*primcount], data->basis[2*primcount+1]);
        primcount++;
      }

      printf("\n");
      shellcount++;
    }
  }
#endif
  printf("mopacplugin) =================================================================\n");
  for (i=0; i<data->num_basis_atoms; i++) {
    printf("%-8s (%10s)\n\n", data->atoms[i].type, data->basis_set[i].name);
    printf("\n");

    for (j=0; j<data->basis_set[i].numshells; j++) {

      for (k=0; k<data->basis_set[i].shell[j].numprims; k++) {
        printf("%6d   %d %7d %22f%22f\n", j,
               data->basis_set[i].shell[j].type,
               primcount+1,
               data->basis_set[i].shell[j].prim[k].exponent,
               data->basis_set[i].shell[j].prim[k].contraction_coeff);
        primcount++;
      }

      printf("\n");
      shellcount++;
    }
  }
  printf("\n");
  printf(" TOTAL NUMBER OF BASIS SET SHELLS             =%5d\n", data->num_shells);
  printf(" NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =%5d\n", data->wavef_size);
  printf(" NUMBER OF ELECTRONS                          =%5d\n", data->num_electrons);
  printf(" CHARGE OF MOLECULE                           =%5d\n", data->totalcharge);
  printf(" SPIN MULTIPLICITY                            =%5d\n", data->multiplicity);
  printf(" NUMBER OF OCCUPIED ORBITALS (ALPHA)          =%5d\n", data->num_occupied_A);
  printf(" NUMBER OF OCCUPIED ORBITALS (BETA )          =%5d\n", data->num_occupied_B);
  printf(" TOTAL NUMBER OF ATOMS                        =%5i\n", data->numatoms);
  printf("\n");
}
