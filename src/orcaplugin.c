/*
Orca VMD plugin
Authors: Maximilian Scheurer, Marcelo Melo, April 2017
Inspired from gamessplugin.c
 */

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include "qmplugin.h"
#include "unit_conversion.h"
#include "periodic_table.h"

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

#define ORCA4 4
#define ORCA3 3

/* ######################################################## */
/*                 orca specific struct                     */
/*
  might be extended in the future...
 */
/* ######################################################## */

typedef struct {
  int version; /* version = 0 : not supported, do not read the file!
                * version = 1 : version 3.0.0, 3.0.1, 3.0.3
                * version = 2 : version 4.0.0
                */

  int digits[3]; /* The three digits of the orca version number
                  * e.g. 4.0.0 (i.e. digits[0].digits[1].digits[2])
                  * */
} orcadata;


/* ######################################################## */
/*                    static functions                      */
/* ######################################################## */

// Checks if loaded file is really an Orca output file
static int have_orca(qmdata_t *data, orcadata* orca);

static int read_orca_structure(void *mydata, int *optflags, molfile_atom_t *atoms);

// Freeing memory
static void close_orca_read(void *mydata);

// Function for reading timestep independent information: Main Parser
static int parse_static_data(qmdata_t *, int *);

// atom definitions & geometry
static int get_input_structure(qmdata_t *data, orcadata *orca);

// reading coord block
static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit, int *numatoms);

/*************************************************************
 *
 * MAIN ORCA CODE PART
 *
 **************************************************************/

static void* open_orca_read(const char* filename, const char* filetype, int *natoms) {
  FILE* fd;
  qmdata_t *data = NULL;

  #ifdef DEBUGGING
    printf("DEBUG: Open Orca Read called: %s\n", filename);
  #endif

  orcadata* orca;

  // open the orca output files
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

  // orca specific information alloc
  orca = (orcadata *) calloc(1, sizeof(orcadata));
  orca->version = 0;
  data->format_specific_data = orca;

  // filename in qm data struct
  data->file = fd;

  if(have_orca(data, orca)) {
    if (orca->version != 0) {
      printf("orcaplugin) Orca version: %d.%d.%d \n", orca->digits[0],
                                            orca->digits[1],orca->digits[2]);
    } else {
      printf("orcaplugin) Orca version not supported: %d.%d.%d \n", orca->digits[0],
                                            orca->digits[1],orca->digits[2]);
      return NULL;
    }
    if (parse_static_data(data, natoms) == FALSE) {
      return NULL;
    }
  } else {
    printf("orcaplugin) This is not an Orca output file!\n");
    return NULL;
  }

  return data;
}

static int read_orca_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  qmdata_t *data = (qmdata_t *)mydata;
  qm_atom_t *cur_atom;
  molfile_atom_t *atom;
  int i = 0;
  return MOLFILE_SUCCESS;
}



static int have_orca(qmdata_t *data, orcadata* orca) {
  int programLine;
  int versionLine;
  char buffer[BUFSIZ];
  int mainVersion, secondDigit, thirdDigit;
  buffer[0] = '\0';
  programLine = goto_keyline(data->file, "O   R   C   A", NULL);
  if (programLine != 1) {
    return FALSE;
  }

  versionLine = goto_keyline(data->file, "Program Version", NULL);
  // thisline(data->file);
  GET_LINE(buffer, data->file);
  if (strstr(buffer,"Version") != NULL) {
    sscanf(buffer, "%*s %*s %d.%d.%d", &mainVersion, &secondDigit, &thirdDigit);
    #ifdef DEBUGGING
      printf("DEBUG: build: %d.%d.%d\n", mainVersion, secondDigit, thirdDigit);
    #endif
    int build[3] = { mainVersion, secondDigit, thirdDigit };
    for (size_t i = 0; i < 3; i++) {
        orca->digits[i] = build[i];
    }
    switch (mainVersion) {
      case ORCA4:
        orca->version = 2;
        break;
      case ORCA3:
        orca->version = 1;
        break;
      default:
        orca->version = 0;
        break;
    }
  } else {
    PRINTERR;
    return FALSE;
  }



  return TRUE;
}

static int parse_static_data(qmdata_t *data, int* natoms) {
  orcadata *orca = (orcadata *)data->format_specific_data;

  get_input_structure(data, orca);

  *natoms = data->numatoms;
  return TRUE;
}


/**********************************************************
 *
 * Read the input atom definitions and geometry
 *
 **********************************************************/
static int get_input_structure(qmdata_t *data, orcadata *orca) {
  char buffer[BUFSIZ];
  char units[BUFSIZ];
  int numatoms = -1;
  int bohr;
  long filepos;
  filepos = ftell(data->file);

  if (goto_keyline(data->file, "CARTESIAN COORDINATES (ANGSTROEM)", NULL)) {
    GET_LINE(buffer, data->file);
    thisline(data->file);
    bohr = 0;
    // sscanf()
  } else {
    printf("orcaplugin) No cartesian coordinates in ANGSTROEM found.\n");
    return FALSE;
  }

  // skip the ---- line
  eatline(data->file, 1);
  /* Read the coordinate block */
  if (get_coordinates(data->file, &data->atoms, bohr, &numatoms))
    data->num_frames_read = 0;
  else {
    printf("gamessplugin) Bad atom coordinate block!\n");
    return FALSE;
  }

  data->numatoms = numatoms;
  return TRUE;
}


static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit,
                           int *numatoms) {
  int i = 0;
  int growarray = 0;

  if (*numatoms<0) {
    *atoms = (qm_atom_t*)calloc(1, sizeof(qm_atom_t));
    growarray = 1;
  }

  /* Read in the coordinates until an empty line is reached.
   * We expect 5 entries per line */
  while (1) {
    char buffer[BUFSIZ];
    char atname[BUFSIZ];
    float atomicnum;
    float x,y,z, dum;
    int n;
    qm_atom_t *atm;

    GET_LINE(buffer, file);
    // thisline(file);

    /* For FMO there is an additional atom index in the
     * second column. Try both variants: */
    n = sscanf(buffer,"%s %f %f %f",atname,&x,&y,&z);
    // printf("%s\n", atname);
    if (n!=4) {
      // n = sscanf(buffer,"%s %f %f %f %f",atname,&atomicnum,&x,&y,&z);
      break;
    }
    // if (n!=5 && n!=6) break;

    if (growarray && i>0) {
      *atoms = (qm_atom_t*)realloc(*atoms, (i+1)*sizeof(qm_atom_t));
    }
    atm = (*atoms)+i;

    // just get the atomic number from periodic_table.h
    atomicnum = get_pte_idx(atname);

    strncpy(atm->type, atname, sizeof(atm->type));
    atm->atomicnum = floor(atomicnum+0.5); /* nuclear charge */
    printf("coor: %s %d %f %f %f\n", atm->type, atm->atomicnum, x, y, z);

    /* if coordinates are in Bohr convert them to Angstrom */
    if (unit==BOHR) {
      x *= BOHR_TO_ANGS;
      y *= BOHR_TO_ANGS;
      z *= BOHR_TO_ANGS;
    }

    atm->x = x;
    atm->y = y;
    atm->z = z;
    i++;
  }

  /* If file is broken off in the middle of the coordinate block
   * we cannot use this frame. */
  if (*numatoms>=0 && *numatoms!=i) {
    (*numatoms) = i;
    return FALSE;
  }

  (*numatoms) = i;
  return TRUE;
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
  plugin.name = "orca";
  plugin.prettyname = "Orca";
  plugin.author = "Maximilian Scheurer, Marcelo Melo";
  plugin.majorv = 0;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "orca";
  plugin.open_file_read = open_orca_read;
  plugin.read_structure = read_orca_structure;
  plugin.close_file_read = close_orca_read;
  //
  // plugin.read_qm_metadata = read_orca_metadata;
  // plugin.read_qm_rundata  = read_orca_rundata;

// #if vmdplugin_ABIVERSION > 11
//   plugin.read_timestep_metadata    = read_timestep_metadata;
//   plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;
//   plugin.read_timestep = read_timestep;
// #endif

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
static void close_orca_read(void *mydata) {

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
  free(data->filepos_array);

  if (data->basis_set) {
    for(i=0; i<data->num_basis_atoms; i++) {
      for (j=0; j<data->basis_set[i].numshells; j++) {
        free(data->basis_set[i].shell[j].prim);
      }
      free(data->basis_set[i].shell);
    }
    free(data->basis_set);
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
