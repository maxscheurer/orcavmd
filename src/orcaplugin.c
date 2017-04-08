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

// for VMD
static int read_orca_metadata(void *mydata, molfile_qm_metadata_t *metadata);
static int read_orca_rundata(void *mydata, molfile_qm_t *qm_data);

static int read_timestep(void *mydata, int natoms,
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
			 molfile_qm_timestep_t *qm_ts);

static int read_timestep_metadata(void *mydata, molfile_timestep_metadata_t *meta);

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
    printf("orcaplugin) atomicnum[%d] = %d\n", i, atom->atomicnumber);
    #endif

    /* if (data->have_mulliken)
    atom->charge = data->qm_timestep->mulliken_charges[i];
    */
    cur_atom++;
  }

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
    // thisline(data->file);
    // UNITS ARE ANGSTROEM
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
    printf("orcaplugin) Bad atom coordinate block!\n");
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
    // printf("coor: %s %d %f %f %f\n", atm->type, atm->atomicnum, x, y, z);

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
  plugin.read_qm_metadata = read_orca_metadata;
  plugin.read_qm_rundata  = read_orca_rundata;

#if vmdplugin_ABIVERSION > 11
  // plugin.read_timestep_metadata    = read_timestep_metadata;
//   plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;
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


#if vmdplugin_ABIVERSION > 11

/***********************************************************
 * Provide non-QM metadata for next timestep.
 * Required by the plugin interface.
 ***********************************************************/
static int read_timestep_metadata(void *mydata,
                                  molfile_timestep_metadata_t *meta) {
  meta->count = -1;
  meta->has_velocities = 0;

  return MOLFILE_SUCCESS;
}
/***********************************************************
 *
 * This function provides the data of the next timestep.
 * Here we actually don't read the data from file, that had
 * to be done already upon calling read_timestep_metadata().
 * Instead we copy the stuff from the local data structure
 * into the one's provided by VMD.
 *
 ***********************************************************/
static int read_timestep(void *mydata, int natoms,
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
			 molfile_qm_timestep_t *qm_ts)
{
  printf("READING TIMESTEP\n");
  qmdata_t *data = (qmdata_t *)mydata;
  qm_timestep_t *cur_ts;
  int offset;
  int i = 0;
  int num_charge_sets = 0;

  if (data->trajectory_done == TRUE) return MOLFILE_ERROR;

  printf("copying coords.\n");
  /* copy the coordinates */
  for (i=0; i<natoms; i++) {
    ts->coords[3*i  ] = data->atoms[i].x;
    ts->coords[3*i+1] = data->atoms[i].y;
    ts->coords[3*i+2] = data->atoms[i].z;
  }

    printf("ts pointer.\n");
  /* get a convenient pointer to the current qm timestep */
  cur_ts = data->qm_timestep+data->num_frames_sent;

  /* store the SCF energies */
  for (i=0; i<cur_ts->num_scfiter; i++) {
    qm_ts->scfenergies[i] = cur_ts->scfenergies[i];
  }

  /* store gradients */
  if (cur_ts->gradient) {
    for (i=0; i<3*natoms; i++) {
      qm_ts->gradient[i] = cur_ts->gradient[i];
    }
  }

  /* store charge sets*/
  if (cur_ts->have_mulliken) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->mulliken_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_MULLIKEN;
    num_charge_sets++;
  }

  if (cur_ts->have_lowdin) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->lowdin_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_LOWDIN;
    num_charge_sets++;
  }
  if (cur_ts->have_esp) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->esp_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_ESP;
    num_charge_sets++;
  }


  /* store the wave function and orbital energies */
  if (cur_ts->wave) {
    for (i=0; i<cur_ts->numwave; i++) {
      qm_wavefunction_t *wave = &cur_ts->wave[i];
      qm_ts->wave[i].type         = wave->type;
      qm_ts->wave[i].spin         = wave->spin;
      qm_ts->wave[i].excitation   = wave->exci;
      qm_ts->wave[i].multiplicity = wave->mult;
      qm_ts->wave[i].energy       = wave->energy;
      strncpy(qm_ts->wave[i].info, wave->info, MOLFILE_BUFSIZ);

      if (wave->wave_coeffs) {
        memcpy(qm_ts->wave[i].wave_coeffs, wave->wave_coeffs,
               wave->num_orbitals*data->wavef_size*sizeof(float));
      }
      if (wave->orb_energies) {
        memcpy(qm_ts->wave[i].orbital_energies, wave->orb_energies,
               wave->num_orbitals*sizeof(float));
      }
      if (wave->has_occup) {
        memcpy(qm_ts->wave[i].occupancies, wave->orb_occupancies,
               wave->num_orbitals*sizeof(float));
      }
    }
  }

  if (data->runtype == MOLFILE_RUNTYPE_ENERGY ||
      data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    /* We have only a single point */
    data->trajectory_done = TRUE;
  }

  data->num_frames_sent++;

  return MOLFILE_SUCCESS;
}
#endif

/*****************************************************
 *
 * provide VMD with the sizes of the QM related
 * data structure arrays that need to be made
 * available
 *
 *****************************************************/
static int read_orca_metadata(void *mydata,
    molfile_qm_metadata_t *metadata) {

  qmdata_t *data = (qmdata_t *)mydata;

  if (data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    metadata->ncart = (3*data->numatoms);
    metadata->nimag = data->nimag;

    if (data->have_internals) {
      metadata->nintcoords = data->nintcoords;
    } else {
      metadata->nintcoords = 0;
    }
  }
  else {
    metadata->ncart = 0;
    metadata->nimag = 0;
    metadata->nintcoords = 0;
  }

  /* orbital data */
  metadata->num_basis_funcs = data->num_basis_funcs;
  metadata->num_basis_atoms = data->num_basis_atoms;
  metadata->num_shells      = data->num_shells;
  metadata->wavef_size      = data->wavef_size;

#if vmdplugin_ABIVERSION > 11
  /* system and run info */
  metadata->have_sysinfo = 1;

  /* hessian info */
  metadata->have_carthessian = data->have_cart_hessian;
  metadata->have_inthessian  = data->have_int_hessian;

  /* normal mode info */
  metadata->have_normalmodes = data->have_normal_modes;
#endif

  return MOLFILE_SUCCESS;
}


/******************************************************
 *
 * Provide VMD with the static (i.e. non-trajectory)
 * data. That means we are filling the molfile_plugin
 * data structures.
 *
 ******************************************************/
static int read_orca_rundata(void *mydata,
                               molfile_qm_t *qm_data) {

  qmdata_t *data = (qmdata_t *)mydata;
  int i, j;
  int ncart;
  molfile_qm_hessian_t *hessian_data = &qm_data->hess;
  molfile_qm_basis_t   *basis_data   = &qm_data->basis;
  molfile_qm_sysinfo_t *sys_data     = &qm_data->run;

  /* fill in molfile_qm_hessian_t */
  if (data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    ncart = 3*data->numatoms;

    /* Hessian matrix in cartesian coordinates */
    if (data->have_cart_hessian) {
      for (i=0; i<ncart; i++) {
        for (j=0; j<=i; j++) {
          hessian_data->carthessian[ncart*i+j] = data->carthessian[ncart*i+j];
          hessian_data->carthessian[ncart*j+i] = data->carthessian[ncart*i+j];
        }
      }
    }

    /* Hessian matrix in internal coordinates */
    if (data->have_int_hessian) {
      for (i=0; i<(data->nintcoords)*(data->nintcoords); i++) {
        hessian_data->inthessian[i] = data->inthessian[i];
      }
    }

    /* wavenumbers, intensities, normal modes */
    if (data->have_normal_modes) {
      for (i=0; i<ncart*ncart; i++) {
        hessian_data->normalmodes[i] = data->normal_modes[i];
      }
      for (i=0; i<ncart; i++) {
        hessian_data->wavenumbers[i] = data->wavenumbers[i];
        hessian_data->intensities[i] = data->intensities[i];
      }
    }

    /* imaginary modes */
    for (i=0; i<data->nimag; i++) {
      /*printf("imag_modes[%d]=%d\n", i, data->imag_modes[i]);*/
      hessian_data->imag_modes[i] = data->imag_modes[i];
    }
  }

  /* fill in molfile_qm_sysinfo_t */
  sys_data->runtype = data->runtype;
  sys_data->scftype = data->scftype;
  sys_data->nproc   = data->nproc;
  sys_data->num_electrons  = data->num_electrons;
  sys_data->totalcharge    = data->totalcharge;
  sys_data->num_occupied_A = data->num_occupied_A;
  sys_data->num_occupied_B = data->num_occupied_B;
  sys_data->status         = data->status;


  strncpy(sys_data->basis_string, data->basis_string,
          sizeof(sys_data->basis_string));

  sys_data->memory = 0; /* XXX fixme */

  strncpy(sys_data->runtitle, data->runtitle, sizeof(sys_data->runtitle));
  strncpy(sys_data->geometry, data->geometry, sizeof(sys_data->geometry));
  strncpy(sys_data->version_string, data->version_string,
          sizeof(sys_data->version_string));

#if vmdplugin_ABIVERSION > 11
  /* fill in molfile_qm_basis_t */
  if (data->num_basis_funcs) {
    for (i=0; i<data->num_basis_atoms; i++) {
      basis_data->num_shells_per_atom[i] = data->num_shells_per_atom[i];
      basis_data->atomic_number[i] = data->atomicnum_per_basisatom[i];
    }

    for (i=0; i<data->num_shells; i++) {
      basis_data->num_prim_per_shell[i] = data->num_prim_per_shell[i];
      basis_data->shell_types[i] = data->shell_types[i];
    }

    for (i=0; i<2*data->num_basis_funcs; i++) {
      basis_data->basis[i] = data->basis[i];
    }

    for (i=0; i<3*data->wavef_size; i++) {
      basis_data->angular_momentum[i] = data->angular_momentum[i];
    }
  }
#endif

  return MOLFILE_SUCCESS;
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
