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

// analyze trajectory, getting number of frames, file positions etc.
static int analyze_traj(qmdata_t *data, orcadata *orca);

// atom definitions & geometry
static int get_input_structure(qmdata_t *data, orcadata *orca);

// reading coord block
static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit, int *numatoms);

// reading first trajectory frame
static int read_first_frame(qmdata_t *data);

int get_basis(qmdata_t *data);

// main routine for extracting "step"-dependent info from the output file
static int get_traj_frame(qmdata_t *data, qm_atom_t *atoms, int natoms);

static int get_scfdata(qmdata_t *data, qm_timestep_t *ts);

static int check_add_wavefunctions(qmdata_t *data, qm_timestep_t *ts);

static int fill_basis_arrays(qmdata_t *data);

static int shelltype_int(char type);

// for VMD
static int read_orca_metadata(void *mydata, molfile_qm_metadata_t *metadata);
static int read_orca_rundata(void *mydata, molfile_qm_t *qm_data);
static int read_qm_timestep_metadata(void *mydata, molfile_qm_timestep_metadata_t *meta);

static int read_timestep(void *mydata, int natoms,
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
			 molfile_qm_timestep_t *qm_ts);

static int read_timestep_metadata(void *mydata, molfile_timestep_metadata_t *meta);

static void print_input_data(qmdata_t *data);
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

  if (!get_input_structure(data, orca)) return FALSE;

  if (!get_basis(data)) return FALSE;

  if (!analyze_traj(data, orca)) {
    printf("orcaplugin) WARNING: Truncated or abnormally terminated file!\n\n");
  }

  *natoms = data->numatoms;

  read_first_frame(data);

  // print_input_data(data);

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

int get_basis(qmdata_t *data) {

  orcadata *orca = (orcadata *)data->format_specific_data;

  char buffer[BUFSIZ];
  char word[4][BUFSIZ];
  int i = 0;
  int success = 0;
  int numread, numshells;
  shell_t *shell;
  long filepos;

  // TODO: needed for PM3
  // if (!strcmp(data->gbasis, "MNDO") ||
  //     !strcmp(data->gbasis, "AM1")  ||
  //     !strcmp(data->gbasis, "PM3")) {
  //   /* Semiempirical methods are based on STOs.
  //    * The only parameter we need for orbital rendering
  //    * are the exponents zeta for S, P, D,... shells for
  //    * each atom. Since GAMESS doesn't print these values
  //    * we skip reading the basis set and hardcode the
  //    * parameters in tables in VMD. */
  //   return TRUE;
  // }

  /* Search for "ATOMIC BASIS SET" line */
  if (pass_keyline(data->file, "BASIS SET IN INPUT FORMAT", NULL) != FOUND ) {
    printf("orcaplugin) No basis set found!\n");
    return FALSE;
  }

  /* initialize buffers */
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';


  /* skip the next 3 lines */
  // eatline(data->file, 2);

  /* Allocate space for the basis for all atoms */
  /* When the molecule is symmetric the actual number atoms with
   * a basis set could be smaller */
  data->basis_set = (basis_atom_t*)calloc(data->numatoms, sizeof(basis_atom_t));

  filepos = ftell(data->file);
  i = 0; /* basis atom counter */
  int finished = FALSE;
  while (!finished) {
    printf("Trying to read bf. \n");
    if (pass_keyline(data->file, "Basis set for element", NULL) == FOUND ) {
      GET_LINE(buffer, data->file);
      numread = sscanf(buffer,"%s %s",&word[0][0], &word[1][0]);
      printf("New element found: %s\n", &word[1][0]);
      int elementCompleted = 0;

      prim_t *prim = NULL;

      shell = (shell_t*)calloc(1, sizeof(shell_t));
      numshells = 0;
      int readingShell = 0;
      int primcounter;

      // this is very sloppy at the moment...
      // for PM3 etc. Orca prints the bf per atom...
      // float exponent = 0.0;
      // float contract = 0.0;

      while(!elementCompleted) {
        GET_LINE(buffer, data->file);
        numread = sscanf(buffer,"%s %s %s",&word[0][0], &word[1][0],&word[2][0]);
        printf("numread: %d -- %s %s %s \n",numread, &word[0][0], &word[1][0],&word[2][0]);
        switch (numread) {
          case 1:
            if (strcmp(trimleft(trimright(&word[0][0])), "end")) {
              printf("Section ended. \n");
              elementCompleted = 1;
              break;
            }
          case 2:
            shell[numshells].numprims = atoi(trimleft(trimright(&word[1][0])));
            shell[numshells].type = shelltype_int(word[0][0]);
            printf("orcaplugin) Type: %d NPrims: %d\n", shell[numshells].type, shell[numshells].numprims);
            primcounter = 0;
            prim = (prim_t*)calloc(1, sizeof(prim_t));
            shell[numshells].prim = prim;
            numshells++;
            if (numshells) {
              shell = (shell_t*)realloc(shell, (numshells+1)*sizeof(shell_t));
            }
            break;
          case 3:
            printf("coeffients.\n");
            prim[primcounter].exponent = atof(&word[1][0]);
            prim[primcounter].contraction_coeff = atof(&word[2][0]);
            primcounter++;
            if (primcounter) {
              prim = (prim_t*)realloc(prim, (primcounter+1)*sizeof(prim_t));
            }
            break;
          default:
            printf("unkown line in bf. \n");
            elementCompleted = 1;
            break;
        }
      }
      printf("Number of shells: %d \n", numshells);
    } else {
      finished = TRUE;
      printf("orcaplugin) Reading basis set finished! \n");
    }
  }




  return FALSE;

  do {
    prim_t *prim = NULL;
    char shelltype;
    int numprim = 0;
    int icoeff = 0;
    filepos = ftell(data->file);
    GET_LINE(buffer, data->file);

    /* Count the number of relevant words in the line. */
    numread = sscanf(buffer,"%s %s %s %s",&word[0][0], &word[1][0],
           &word[2][0], &word[3][0]);

    switch (numread) {
      case 1:
        /* Next atom */
        strcpy(data->basis_set[i].name, &word[0][0]);

        /* skip initial blank line */
        eatline(data->file, 1);

        /* read the basis set for the current atom */
        shell = (shell_t*)calloc(1, sizeof(shell_t));
        numshells = 0;

        do {
          filepos = ftell(data->file);
          // TODO: edit for Orca
          // numprim = read_shell_primitives(data, &prim, &shelltype, icoeff, gms->have_pcgamess);

          if (numprim>0) {
            /* make sure we have eiter S, L, P, D, F or G shells */
            if ( (shelltype!='S' && shelltype!='L' && shelltype!='P' &&
                  shelltype!='D' && shelltype!='F' && shelltype!='G') ) {
              printf("orcaplugin) WARNING ... %c shells are not supported \n", shelltype);
            }

            /* create new shell */
            if (numshells) {
              shell = (shell_t*)realloc(shell, (numshells+1)*sizeof(shell_t));
            }
            shell[numshells].numprims = numprim;
            /* assign a numeric shell type */
            shell[numshells].type = shelltype_int(shelltype);
            shell[numshells].prim = prim;
            data->num_basis_funcs += numprim;

            /* We split L-shells into one S and one P-shell.
             * I.e. for L-shells we have to go back read the shell again
             * this time using the second contraction coefficients. */
            if (shelltype=='L' && !icoeff) {
              fseek(data->file, filepos, SEEK_SET);
              icoeff++;
            } else if (shelltype=='L' && icoeff) {
              shell[numshells].type = SP_P_SHELL;
              icoeff = 0;  /* reset the counter */
            }

            numshells++;
          }
        } while (numprim);

        /* store shells in atom */
        data->basis_set[i].numshells = numshells;
        data->basis_set[i].shell = shell;

        /* Update total number of basis functions */
        data->num_shells += numshells;
        i++;

        /* go back one line so that we can read the name of the
         * next atom */
        fseek(data->file, filepos, SEEK_SET);

        break;

      case 4:
        /* this is the very end of the basis set */
        // if(gms->have_pcgamess){
        //     if (!strcmp(&word[0][0],"TOTAL")  &&
        //         !strcmp(&word[1][0],"NUMBER") &&
        //         !strcmp(&word[2][0],"OF")     &&
        //         !strcmp(&word[3][0],"SHELLS")) {
        //       success = 1;
        //       /* go back one line so that get_basis_stats()
        //          can use this line as a keystring. */
        //       fseek(data->file, filepos, SEEK_SET);
        //     }
        // }
        // else {
        //     if (!strcmp(&word[0][0],"TOTAL")  &&
        //         !strcmp(&word[1][0],"NUMBER") &&
        //         !strcmp(&word[2][0],"OF")     &&
        //         !strcmp(&word[3][0],"BASIS")) {
        //       success = 1;
        //       /* go back one line so that get_basis_stats()
        //          can use this line as a keystring. */
        //       fseek(data->file, filepos, SEEK_SET);
        //     }
        // }
        break;
    }

  } while (!success);


  printf("orcaplugin) Parsed %d uncontracted basis functions for %d atoms.\n",
         data->num_basis_funcs, i);

  data->num_basis_atoms = i;


  /* allocate and populate flat arrays needed for molfileplugin */
  return fill_basis_arrays(data);
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

/* Read the first trajectory frame. */
static int read_first_frame(qmdata_t *data) {
  /* The angular momentum is populated in get_wavefunction
   * which is called by get_traj_frame(). We have obtained
   * the array size wavef_size already from the basis set
   * statistics */
  data->angular_momentum = (int*)calloc(3*data->wavef_size, sizeof(int));

  /* Try reading the first frame.
   * If there is only one frame then also read the
   * final wavefunction. */
  if (!get_traj_frame(data, data->atoms, data->numatoms)) {
    return FALSE;
  }

  return TRUE;
}



/******************************************************
 *
 * this function extracts the trajectory information
 * from the output file
 *
 * *****************************************************/
static int get_traj_frame(qmdata_t *data, qm_atom_t *atoms,
                          int natoms) {
  orcadata *orca = (orcadata *)data->format_specific_data;
  qm_timestep_t *cur_ts;
  char buffer[BUFSIZ];
  char word[BUFSIZ];
  int units;
  buffer[0] = '\0';
  word[0]   = '\0';

  printf("orcaplugin) Timestep %d:\n", data->num_frames_read);
  printf("orcaplugin) ============\n");


  // debugging the trajectory reading file positions
  // printf("nfread: %d \n", data->num_frames_read);
  // if (!data->filepos_array) {
  //   printf("filepos array empty!!!\n");
  //   return FALSE;
  // }

  fseek(data->file, data->filepos_array[data->num_frames_read], SEEK_SET);

  /*
  * distinguish between job types
  * at the moment, only Single Points will work
  * lines 2840 - 3122 in gamessplugin.c
   */

  /* Read the coordinate block */
  // if (data->runtype==MOLFILE_RUNTYPE_OPTIMIZE ||
  //     data->runtype==MOLFILE_RUNTYPE_SADPOINT) {
  //   goto_keyline(data->file, "COORDINATES OF ALL ATOMS", NULL);
  //   /* get the units */
  //   GET_LINE(buffer, data->file);
  //   sscanf(buffer, " COORDINATES OF ALL ATOMS ARE %s", word);
  //   units = !strcmp(word, "(BOHR)");
  //   eatline(data->file, 2);
  //
  //   if (!get_coordinates(data->file, &data->atoms, units, &natoms)) {
  //     printf("orcaplugin) Couldn't find coordinates for timestep %d\n", data->num_frames_read);
  //   }
  // }
  // else if (data->runtype==MOLFILE_RUNTYPE_SURFACE) {
  //   if (pass_keyline(data->file, "HAS ENERGY VALUE",
  //                    "...... END OF ONE-ELECTRON INTEGRALS ......")
  //       == FOUND) {
  //     /* Read the coordinate block following
  //      * ---- SURFACE MAPPING GEOMETRY ---- */
  //     int i, n;
  //     for (i=0; i<natoms; i++) {
  //       char atname[BUFSIZ];
  //       float x,y,z;
  //       GET_LINE(buffer, data->file);
  //       n = sscanf(buffer,"%s %f %f %f", atname, &x,&y,&z);
  //       if (n!=4 || strcmp(atname, data->atoms[i].type)) break;
  //       data->atoms[i].x = x;
  //       data->atoms[i].y = y;
  //       data->atoms[i].z = z;
  //     }
  //     if (i!=natoms) {
  //       printf("orcaplugin) Couldn't read surface mapping geometry for timestep %d\n", data->num_frames_read);
  //     }
  //   }
  //   else {
  //     /* Read the coordinate block following
  //      * ATOM      ATOMIC                      COORDINATES (BOHR) */
  //     goto_keyline(data->file, "ATOM      ATOMIC", NULL);
  //     /* get the units */
  //     GET_LINE(buffer, data->file);
  //     sscanf(buffer, " ATOM      ATOMIC                      COORDINATES %s", word);
  //     units = !strcmp(word, "(BOHR)");
  //     eatline(data->file, 1);
  //
  //     if (!get_coordinates(data->file, &data->atoms, units, &natoms)) {
  //       printf("orcaplugin) Couldn't find coordinates for timestep %d\n", data->num_frames_read);
  //     }
  //   }
  // }
  // /* XXX could merge this with OPTIMIZE/SADPOINT */
  // else if (data->runtype==MOLFILE_RUNTYPE_MEX) {
  //   int numuniqueatoms = natoms;
  //   goto_keyline(data->file, "COORDINATES OF SYMMETRY UNIQUE ATOMS", NULL);
  //   /* get the units */
  //   GET_LINE(buffer, data->file);
  //   sscanf(buffer, " COORDINATES OF SYMMETRY UNIQUE ATOMS ARE %s", word);
  //   units = !strcmp(word, "(BOHR)");
  //   eatline(data->file, 2);
  //   if (!get_coordinates(data->file, &data->atoms, units, &numuniqueatoms)) {
  //     printf("orcaplugin) Expanding symmetry unique coordinates for timestep %d\n", data->num_frames_read);
  //
  //     /* Create images of symmetry unique atoms so that we have
  //      * the full coordinate set. */
  //     symmetry_expand(&data->atoms, numuniqueatoms, natoms,
  //                     data->pointgroup, data->naxis);
  //   }
  // }

  /* get a convenient pointer to the current qm timestep */
  // cur_ts = data->qm_timestep + data->num_frames_read;

  /* read the SCF energies */
  if (get_scfdata(data, cur_ts) == FALSE) {
    printf("orcaplugin) Couldn't find SCF iterations for timestep %d\n",
           data->num_frames_read);
  }

  /* Try reading canonical alpha/beta wavefunction */
  check_add_wavefunctions(data, cur_ts);

  /* Read population analysis (Mulliken and Lowdin charges)
   * only if wasn't read already while parsing the final
   * property section. Otherwise we would potentially
   * overwrite the data with empty fields. */
  // if (!cur_ts->have_mulliken &&
  //     get_population(data, cur_ts)) {
  //   printf("orcaplugin) Mulliken/Loewdin charges found\n");
  // }


  /* Read the energy gradients (=forces on atoms) */
  // if (get_gradient(data, cur_ts)) {
  //   printf("orcaplugin) Energy gradient found.\n");
  // }

  /* If this is the last frame of the trajectory and the file
   * wasn't truncated and the program didn't terminate
   * abnormally then read the final wavefunction. */
  if ((data->runtype == MOLFILE_RUNTYPE_OPTIMIZE ||
       data->runtype == MOLFILE_RUNTYPE_SADPOINT) &&
      (data->num_frames_read+1 == data->num_frames &&
       (data->status == MOLFILE_QMSTATUS_UNKNOWN ||
        data->status == MOLFILE_QMSTATUS_OPT_CONV ||
        data->status == MOLFILE_QMSTATUS_OPT_NOT_CONV))) {

    /* We need to jump over the end of the trajectory because
     * this is also the keystring for get_wavefunction() to
     * bail out. */
    if (data->status == MOLFILE_QMSTATUS_OPT_CONV ||
        data->status == MOLFILE_QMSTATUS_OPT_NOT_CONV) {
      fseek(data->file, data->end_of_traj, SEEK_SET);
    }

    /* Try to read final wavefunction and orbital energies
     * A preexisting canonical wavefunction for this timestep
     * with the same characteristics (spin, exci, info) will
     * be overwritten by the final wavefuntion if it has more
     * orbitals. */
    check_add_wavefunctions(data, cur_ts);
  }

  data->num_frames_read++;

  return TRUE;
}

static int get_scfdata(qmdata_t *data, qm_timestep_t *ts) {
  return TRUE;
}

/*********************************************************
 *
 * Reads a set of wavefunctions for the current timestep.
 * These are typically the alpha and beta spin wavefunctions
 * or the MCSCF natural and optimized orbitals or the GVB
 * canonical orbitals and geminal pairs.
 *
 **********************************************************/
 // 3461
static int check_add_wavefunctions(qmdata_t *data, qm_timestep_t *ts) {
  return 0;
}

/******************************************************
 *
 * Populate the flat arrays containing the basis
 * set data.
 *
 ******************************************************/
static int fill_basis_arrays(qmdata_t *data) {
  orcadata *orca = (orcadata *)data->format_specific_data;
  int i, j, k;
  int shellcount = 0;
  int primcount = 0;

  float *basis;
  int *num_shells_per_atom;
  int *num_prim_per_shell;
  int *shell_types;
  int *atomicnum_per_basisatom;

  /* Count the total number of primitives which
   * determines the size of the basis array. */
  for(i=0; i<data->num_basis_atoms; i++) {
    for (j=0; j<data->basis_set[i].numshells; j++) {
      primcount += data->basis_set[i].shell[j].numprims;
    }
  }

  /* reserve space for pointer to array containing basis
   * info, i.e. contraction coeficients and expansion
   * coefficients; need 2 entries per basis function, i.e.
   * exponent and contraction coefficient; also,
   * allocate space for the array holding the orbital symmetry
   * information per primitive Gaussian.
   * Finally, initialize the arrays holding the number of
   * shells per atom and the number of primitives per shell*/
  basis = (float *)calloc(2*primcount,sizeof(float));

  /* make sure memory was allocated properly */
  if (basis == NULL) {
    PRINTERR;
    return FALSE;
  }

  shell_types = (int *)calloc(data->num_shells, sizeof(int));

  /* make sure memory was allocated properly */
  if (shell_types == NULL) {
    PRINTERR;
    return FALSE;
  }

  num_shells_per_atom = (int *)calloc(data->num_basis_atoms, sizeof(int));

  /* make sure memory was allocated properly */
  if (num_shells_per_atom == NULL) {
    PRINTERR;
    return FALSE;
  }

  num_prim_per_shell = (int *)calloc(data->num_shells, sizeof(int));

  /* make sure memory was allocated properly */
  if (num_prim_per_shell == NULL) {
    PRINTERR;
    return FALSE;
  }

  atomicnum_per_basisatom = (int *)calloc(data->num_basis_atoms, sizeof(int));

  /* make sure memory was allocated properly */
  if (atomicnum_per_basisatom == NULL) {
    PRINTERR;
    return FALSE;
  }


  /* store pointers in struct qmdata_t */
  data->basis = basis;
  data->shell_types = shell_types;
  data->num_shells_per_atom = num_shells_per_atom;
  data->num_prim_per_shell = num_prim_per_shell;
  data->atomicnum_per_basisatom = atomicnum_per_basisatom;

  /* Go through all basis set atoms and try to assign the
   * atomic numbers. The basis set atoms are specified by
   * name strings (the same as in the coordinate section,
   * except for FMO calcs.) and we try to match the names
   * from the two lists. The basis set atom list is symmetry
   * unique while the coordinate atom list is complete.*/
  primcount = 0;
  for (i=0; i<data->num_basis_atoms; i++) {
    int found = 0;

    /* For this basis atom find a matching atom from the
     * coordinate atom list. */
    for(j=0; j<data->numatoms; j++) {
      char basisname[BUFSIZ];
      strcpy(basisname, data->basis_set[i].name);

      /* for FMO calculations we have to strip the "-n" tail
       * of the basis atom name. */
      // if (gms->have_fmo) {
      //   *strchr(basisname, '-') = '\0';
      // }

      if (!strcmp(data->atoms[j].type, basisname)) {
        found = 1;
        break;
      }
    }
    if (!found) {
      printf("gamessplugin) WARNING: Couldn't find atomic number for basis set atom %s\n",
             data->basis_set[i].name);
      data->basis_set[i].atomicnum = 0;
      atomicnum_per_basisatom[i] = 0;
    } else {
      /* assign atomic number */
      data->basis_set[i].atomicnum = data->atoms[j].atomicnum;
      atomicnum_per_basisatom[i]   = data->atoms[j].atomicnum;
    }
    num_shells_per_atom[i] = data->basis_set[i].numshells;

    for (j=0; j<data->basis_set[i].numshells; j++) {
      shell_types[shellcount] = data->basis_set[i].shell[j].type;
      num_prim_per_shell[shellcount] = data->basis_set[i].shell[j].numprims;

      for (k=0; k<data->basis_set[i].shell[j].numprims; k++) {
        basis[2*primcount  ] = data->basis_set[i].shell[j].prim[k].exponent;
        basis[2*primcount+1] = data->basis_set[i].shell[j].prim[k].contraction_coeff;
        primcount++;
      }
      shellcount++;
    }
  }

  return TRUE;
}


/**************************************************
 *
 * Convert shell type from char to int.
 *
 ************************************************ */
static int shelltype_int(char type) {
  int shelltype;

  switch (type) {
    case 'L':
      /* SP_P shells are assigned in get_basis() */
      shelltype = SP_S_SHELL;
      break;
    case 'S':
      shelltype = S_SHELL;
      break;
    case 'P':
      shelltype = P_SHELL;
      break;
    case 'D':
      shelltype = D_SHELL;
      break;
    case 'F':
      shelltype = F_SHELL;
      break;
    case 'G':
      shelltype = G_SHELL;
      break;
    default:
      shelltype = UNK_SHELL;
      break;
  }

  return shelltype;
}

/* Analyze the trajectory.
 * Read the parameters controlling geometry search and
 * find the end of the trajectory, couinting the frames
 * on the way. Store the filepointer for the beginning of
 * each frame in *filepos_array. */
static int analyze_traj(qmdata_t *data, orcadata *orca) {
  char buffer[BUFSIZ], nserch[BUFSIZ];
  char *line;
  long filepos;
  filepos = ftell(data->file);

  data->filepos_array = (long* )calloc(1, sizeof(long ));

  /* currently, only one frame is supported!
   * lines 3130-3348 in gamessplugin.c
   */
  if (TRUE) {
    /* We have just one frame */
    data->num_frames = 1;
    pass_keyline(data->file, "Single Point Calculation", NULL);
    data->filepos_array[0] = ftell(data->file);

    /* Check wether SCF has converged */
    // if (pass_keyline(data->file,
    //                  "SCF IS UNCONVERGED, TOO MANY ITERATIONS",
    //                  "ENERGY COMPONENTS")==FOUND) {
    //   printf("orcaplugin) SCF IS UNCONVERGED, TOO MANY ITERATIONS\n");
    //   data->status = MOLFILE_QMSTATUS_SCF_NOT_CONV;
    // } else {
    //   data->status = MOLFILE_QMSTATUS_OPT_CONV;
    //   fseek(data->file, data->filepos_array[0], SEEK_SET);
    // }

    pass_keyline(data->file, "FINAL SINGLE POINT ENERGY", NULL);
    data->end_of_traj = ftell(data->file);

    /* Allocate memory for the frame */
    data->qm_timestep = (qm_timestep_t *)calloc(1, sizeof(qm_timestep_t));
    memset(data->qm_timestep, 0, sizeof(qm_timestep_t));

    return TRUE;
  }

  printf("orcaplugin) Analyzing trajectory...\n");
  data->status = MOLFILE_QMSTATUS_UNKNOWN;

  while (0) {
    if (!fgets(buffer, sizeof(buffer), data->file)) break;
    line = trimleft(buffer);

    /* at this point we have to distinguish between
     * pre="27 JUN 2005 (R2)" and "27 JUN 2005 (R2)"
     * versions since the output format for geometry
     * optimizations has changed */

    // if (gms->version==FIREFLY8POST6695){
    //   strcpy(nserch, "NSERCH=");
    // }
    // else if (gms->version==FIREFLY8PRE6695) {
    //   strcpy(nserch, "1NSERCH=");
    // }
    // else if (gms->version==GAMESSPRE20050627) {
    //   strcpy(nserch, "1NSERCH=");
    // }
    // else if (gms->version==GAMESSPOST20050627) {
    //   strcpy(nserch, "BEGINNING GEOMETRY SEARCH POINT NSERCH=");
    // }

    if (strstr(line, nserch) ||
        strstr(line, "---- SURFACE MAPPING GEOMETRY") ||
        strstr(line, "MINIMUM ENERGY CROSSING POINT SEARCH") ||
        (data->runtype==MOLFILE_RUNTYPE_MEX && strstr(line, "NSERCH=")==line)) {
      printf("orcaplugin) %s", line);

      if (data->num_frames > 0) {
        data->filepos_array = (long*)realloc(data->filepos_array,
                                (data->num_frames+1)*sizeof(long));
      }
      data->filepos_array[data->num_frames] = ftell(data->file);
      if (data->runtype==MOLFILE_RUNTYPE_SURFACE) {
        int ret = goto_keyline(data->file,
                               "ATOM      ATOMIC", "HAS ENERGY VALUE",
                               "---- SURFACE MAPPING GEOMETRY ----", NULL);
        if (ret>0 && ret<3 &&
            (have_keyline(data->file, "...... END OF ONE-ELECTRON INTEGRALS ......",
                          "---- SURFACE MAPPING GEOMETRY ----") ||
             have_keyline(data->file, "... DONE WITH POTENTIAL SURFACE SCAN",
                          "---- SURFACE MAPPING GEOMETRY ----"))) {
          data->num_frames++;
        }
      }
      else if (pass_keyline(data->file, "COORDINATES OF",
                            "BEGINNING GEOMETRY SEARCH POINT NSERCH=")==FOUND)
      {
        /* Make sure that we have at least a complete coordinate
           block in order to consider this a new frame. */
        if (have_keyline(data->file, "INTERNUCLEAR DISTANCES",
                         "1 ELECTRON INTEGRALS") ||
            have_keyline(data->file, "1 ELECTRON INTEGRALS",
                         "BEGINNING GEOMETRY SEARCH POINT NSERCH=")) {
          data->num_frames++;
        }
      }
    }
    else if (strstr(line, "***** EQUILIBRIUM GEOMETRY LOCATED") ||
             strstr(line, "... DONE WITH POTENTIAL SURFACE SCAN")) {
      printf("orcaplugin) ==== End of trajectory (%d frames) ====\n",
             data->num_frames);
      data->status = MOLFILE_QMSTATUS_OPT_CONV;
      break;
    }
    else if (strstr(line, "***** FAILURE TO LOCATE STATIONARY POINT,")) {
      printf("orcaplugin) %s\n", line);
      if (strstr(strchr(line, ','), "SCF HAS NOT CONVERGED")) {
        data->status = MOLFILE_QMSTATUS_SCF_NOT_CONV;
        break;
      }
      else if (strstr(strchr(line, ','), "TOO MANY STEPS TAKEN")) {
        data->status = MOLFILE_QMSTATUS_OPT_NOT_CONV;
        break;
      }
    }
  }

  data->end_of_traj = ftell(data->file);
  fseek(data->file, filepos, SEEK_SET);

  if (data->status == MOLFILE_QMSTATUS_UNKNOWN) {
    /* We didn't find any of the regular key strings,
     * the run was most likely broken off and we have an
     * incomplete file. */
    data->status = MOLFILE_QMSTATUS_FILE_TRUNCATED;
  }


  /* Allocate memory for all frames */
  data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames,
                                              sizeof(qm_timestep_t));
  memset(data->qm_timestep, 0, data->num_frames*sizeof(qm_timestep_t));


  if (data->status == MOLFILE_QMSTATUS_SCF_NOT_CONV ||
      data->status == MOLFILE_QMSTATUS_FILE_TRUNCATED) {
    return FALSE;
  }

  return TRUE;
}



/***********************************************************
 * Provide QM metadata for next timestep.
 * This actually triggers reading the entire next timestep
 * since we have to parse the whole timestep anyway in order
 * to get the metadata. So we store the read data locally
 * and hand them to VMD when requested by read_timestep().
 *
 ***********************************************************/
static int read_qm_timestep_metadata(void *mydata,
                                    molfile_qm_timestep_metadata_t *meta) {
  int have = 0;

  qmdata_t *data = (qmdata_t *)mydata;

  meta->count = -1; /* Don't know the number of frames yet */

  if (data->num_frames_read > data->num_frames_sent) {
    have = 1;
  }
  else if (data->num_frames_read < data->num_frames) {
    printf("orcaplugin) Probing timestep %d\n", data->num_frames_read);

    have = get_traj_frame(data, data->atoms, data->numatoms);
  }

  if (have) {
    int i;
    qm_timestep_t *cur_ts;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep+data->num_frames_sent;

    for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<cur_ts->numwave); i++) {
      meta->num_orbitals_per_wavef[i] = cur_ts->wave[i].num_orbitals;
      meta->has_occup_per_wavef[i]    = cur_ts->wave[i].has_occup;
      meta->has_orben_per_wavef[i]    = cur_ts->wave[i].has_orben;
    }
    meta->wavef_size      = data->wavef_size;
    meta->num_wavef       = cur_ts->numwave;
    meta->num_scfiter     = cur_ts->num_scfiter;
    meta->num_charge_sets = cur_ts->have_mulliken +
      cur_ts->have_lowdin + cur_ts->have_esp;
    if (cur_ts->gradient) meta->has_gradient = TRUE;

  } else {
    meta->has_gradient = FALSE;
    meta->num_scfiter  = 0;
    meta->num_orbitals_per_wavef[0] = 0;
    meta->has_occup_per_wavef[0] = FALSE;
    meta->num_wavef = 0;
    meta->wavef_size = 0;
    meta->num_charge_sets = 0;
    data->trajectory_done = TRUE;
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
  plugin.read_timestep_metadata    = read_timestep_metadata;
  plugin.read_qm_timestep_metadata = read_qm_timestep_metadata;
  plugin.read_timestep = read_timestep;
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
  qmdata_t *data = (qmdata_t *)mydata;
  qm_timestep_t *cur_ts;
  int offset;
  int i = 0;
  int num_charge_sets = 0;

  if (data->trajectory_done == TRUE) {
    printf("orcaplugin) Trajectory done.\n");
    return MOLFILE_ERROR;
  }

  /* copy the coordinates */
  for (i=0; i<natoms; i++) {
    ts->coords[3*i  ] = data->atoms[i].x;
    ts->coords[3*i+1] = data->atoms[i].y;
    ts->coords[3*i+2] = data->atoms[i].z;
    // printf("x: %f y: %f z: %f\n", data->atoms[i].x, data->atoms[i].y, data->atoms[i].z);
  }

  /* get a convenient pointer to the current qm timestep */
  // cur_ts = data->qm_timestep+data->num_frames_sent;
  //
  // /* store the SCF energies */
  // for (i=0; i<cur_ts->num_scfiter; i++) {
  //   qm_ts->scfenergies[i] = cur_ts->scfenergies[i];
  // }
  //
  // /* store gradients */
  // if (cur_ts->gradient) {
  //   for (i=0; i<3*natoms; i++) {
  //     qm_ts->gradient[i] = cur_ts->gradient[i];
  //   }
  // }
  //
  // /* store charge sets*/
  // if (cur_ts->have_mulliken) {
  //   offset = num_charge_sets*data->numatoms;
  //   for (i=0; i<data->numatoms; i++) {
  //     qm_ts->charges[offset+i] = cur_ts->mulliken_charges[i];
  //   }
  //   qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_MULLIKEN;
  //   num_charge_sets++;
  // }
  //
  // if (cur_ts->have_lowdin) {
  //   offset = num_charge_sets*data->numatoms;
  //   for (i=0; i<data->numatoms; i++) {
  //     qm_ts->charges[offset+i] = cur_ts->lowdin_charges[i];
  //   }
  //   qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_LOWDIN;
  //   num_charge_sets++;
  // }
  // if (cur_ts->have_esp) {
  //   offset = num_charge_sets*data->numatoms;
  //   for (i=0; i<data->numatoms; i++) {
  //     qm_ts->charges[offset+i] = cur_ts->esp_charges[i];
  //   }
  //   qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_ESP;
  //   num_charge_sets++;
  // }
  //
  //
  // /* store the wave function and orbital energies */
  // if (cur_ts->wave) {
  //   for (i=0; i<cur_ts->numwave; i++) {
  //     qm_wavefunction_t *wave = &cur_ts->wave[i];
  //     qm_ts->wave[i].type         = wave->type;
  //     qm_ts->wave[i].spin         = wave->spin;
  //     qm_ts->wave[i].excitation   = wave->exci;
  //     qm_ts->wave[i].multiplicity = wave->mult;
  //     qm_ts->wave[i].energy       = wave->energy;
  //     strncpy(qm_ts->wave[i].info, wave->info, MOLFILE_BUFSIZ);
  //
  //     if (wave->wave_coeffs) {
  //       memcpy(qm_ts->wave[i].wave_coeffs, wave->wave_coeffs,
  //              wave->num_orbitals*data->wavef_size*sizeof(float));
  //     }
  //     if (wave->orb_energies) {
  //       memcpy(qm_ts->wave[i].orbital_energies, wave->orb_energies,
  //              wave->num_orbitals*sizeof(float));
  //     }
  //     if (wave->has_occup) {
  //       memcpy(qm_ts->wave[i].occupancies, wave->orb_occupancies,
  //              wave->num_orbitals*sizeof(float));
  //     }
  //   }
  // }
  //
  // if (data->runtype == MOLFILE_RUNTYPE_ENERGY ||
  //     data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    /* We have only a single point */
    data->trajectory_done = TRUE;
  // }

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
  printf("orcaplugin) =================================================================\n");
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
