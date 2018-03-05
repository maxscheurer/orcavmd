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
#include <sstream>
#include "qmplugin.h"
#include "unit_conversion.h"
#include "periodic_table.h"
#include "stofit.h"


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

typedef std::vector<std::vector<std::vector<float>>> MoCoeff;

// Checks if loaded file is really an Mopac output file
static int have_mopac(qmdata_t *data, mopacdata* mopac);
static int parse_static_data(qmdata_t *data, int* natoms);
static void* open_mopac_read(const char* filename, const char* filetype, int *natoms);
static int get_job_info(qmdata_t *data);
static int get_input_structure(qmdata_t *data, mopacdata *mopac);
static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit, int *numatoms);
static int get_wavefunction(qmdata_t *data, qm_timestep_t *ts, qm_wavefunction_t *wf);
static int check_add_wavefunctions(qmdata_t *data, qm_timestep_t *ts);
static int analyze_traj(qmdata_t *data, mopacdata *mopac);

static int get_traj_frame(qmdata_t *data, qm_atom_t *atoms, int natoms);

static void close_mopac_read(void *mydata);
static void print_input_data(qmdata_t *data);
static int get_basis_for_element(std::string element, shell_t* shells);

static int read_mopac_rundata(void *mydata, molfile_qm_t *qm_data);
static int read_mopac_metadata(void *mydata, molfile_qm_metadata_t *metadata);


template<typename Out>
void split(const std::string &s, char delim, Out result);


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

std::string reduce(const std::string& str,
                   const std::string& fill = " ",
                   const std::string& whitespace = " \t\n")
{
    // trim first
    auto result = trim(str, whitespace);

    // replace sub ranges
    auto beginSpace = result.find_first_of(whitespace);
    while (beginSpace != std::string::npos)
    {
        const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
        const auto range = endSpace - beginSpace;

        result.replace(beginSpace, range, fill);

        const auto newStart = beginSpace + fill.length();
        beginSpace = result.find_first_of(whitespace, newStart);
    }

    return result;
}

static int get_job_info(qmdata_t *data) {
  long filepos;
  char buffer[BUFSIZ];
  filepos = ftell(data->file);
  if (goto_keyline(data->file, "CALCULATION DONE", NULL)) {
    if (!(goto_keyline(data->file, "*  VECTORS", NULL) && goto_keyline(data->file, "*  EIGEN", NULL))) {
      printf("mopacplugin) Specify the VECTORS and EIGEN keyword in your input file! Otherwise VMD cannot deal with MOPAC output!\n");
      return FALSE;
    }
  } else {
    printf("mopacplugin) MOPAC calculation did not succeed.\n");
    return FALSE;
  }
  data->runtype = MOLFILE_RUNTYPE_ENERGY;
  // go to top of file again and then to the info section
  rewind(data->file);
  goto_keyline(data->file, "CALCULATION DONE", NULL);
  eatline(data->file, 1);

  std::vector<std::string> inputFile;
  int endOfInput = 0;
  while(!endOfInput) {
    GET_LINE(buffer, data->file);
    const std::string test(buffer);
    int lineNumber = 0;
    const std::vector<std::string> jobInfos = split(reduce(test), ' ');
    if (test.find("********") != std::string::npos) {
		  endOfInput = 1;
      break;
    }
    if (jobInfos.size() > 1) {
      const std::string kw = jobInfos[1];
      if (kw.compare("QMMM")) {
        data->runtype = MOLFILE_RUNTYPE_GRADIENT;
      }
    }
	}
  data->multiplicity = 1;
  return TRUE;
}

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
    int list_idx;
    int n;
    qm_atom_t *atm;

    GET_LINE(buffer, file);
    // thisline(file);

    n = sscanf(buffer,"%d %s %f %f %f", &list_idx, atname, &x, &y, &z);
    // printf("%s\n", atname);
    if (n!=5) {
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


static int get_input_structure(qmdata_t *data, mopacdata *mopac) {
  char buffer[BUFSIZ];
  char units[BUFSIZ];
  int numatoms = -1;
  int bohr;
  long filepos;
  filepos = ftell(data->file);

  if (goto_keyline(data->file, "CARTESIAN COORDINATES", NULL)) {
    GET_LINE(buffer, data->file);
    thisline(data->file);
    // UNITS ARE ANGSTROEM
    bohr = 0;
    // sscanf()
  } else {
    printf("mopacplugin) No cartesian coordinates in ANGSTROEM found.\n");
    return FALSE;
  }

  // skip the blank line, line with header and another blank line
  eatline(data->file, 3);
  /* Read the coordinate block */
  if (get_coordinates(data->file, &data->atoms, bohr, &numatoms))
    data->num_frames_read = 0;
  else {
    printf("mopacplugin) Bad atom coordinate block!\n");
    return FALSE;
  }

  data->numatoms = numatoms;
  printf("Number of atoms: %d\n", numatoms);
  return TRUE;
}


static int check_add_wavefunctions(qmdata_t *data, qm_timestep_t *ts) {
  qm_wavefunction_t *wavef;
  int i, n=1;

  if (data->scftype==MOLFILE_SCFTYPE_UHF ||
      data->scftype==MOLFILE_SCFTYPE_GVB ||
      data->scftype==MOLFILE_SCFTYPE_MCSCF) {
    /* Try to read second wavefunction, e.g. spin beta */
    n = 2;
  }

  for (i=0; i<n; i++) {
    /* Allocate memory for new wavefunction */
    wavef = add_wavefunction(ts);

    /* Try to read wavefunction and orbital energies */
    if (get_wavefunction(data, ts, wavef) == FALSE) {
      /* Free the last wavefunction again. */
      del_wavefunction(ts);
#ifdef DEBUGGING_ORCA
      printf("mopacplugin) No canonical wavefunction present for timestep %d\n", data->num_frames_read);
#endif
      break;

    } else {
      char action[32];
      char spinstr[32];
      strcpy(spinstr, "");
      if (data->scftype==MOLFILE_SCFTYPE_UHF) {
        if (wavef->spin==SPIN_BETA) {
          strcat(spinstr, "spin  beta, ");
        } else {
          strcat(spinstr, "spin alpha, ");
        }
      }

      /* The last SCF energy is the energy of this electronic state */
      if (ts->scfenergies) {
        wavef->energy = ts->scfenergies[ts->num_scfiter-1];
      } else {
        wavef->energy = 0.f;
      }

      /* Multiplicity */
      wavef->mult = data->multiplicity;


      /* String telling wether wavefunction was added, updated
       * or ignored. */
      strcpy(action, "added");

      /* If there exists a canonical wavefunction of the same spin
       * we'll replace it */
      if (ts->numwave>1 && wavef->type==MOLFILE_WAVE_CANON) {
        int i, found =-1;
        for (i=0; i<ts->numwave-1; i++) {
          if (ts->wave[i].type==wavef->type &&
              ts->wave[i].spin==wavef->spin &&
              ts->wave[i].exci==wavef->exci &&
              !strncmp(ts->wave[i].info, wavef->info, MOLFILE_BUFSIZ)) {
            found = i;
            break;
          }
        }
        if (found>=0) {
          /* If the new wavefunction has more orbitals we
           * replace the old one for this step. */
          if (wavef->num_orbitals >
              ts->wave[found].num_orbitals) {
            /* Replace existing wavefunction for this step */
            replace_wavefunction(ts, found);
            sprintf(action, "%d updated", found);
          } else {
            /* Delete last wavefunction again */
            del_wavefunction(ts);
            sprintf(action, "matching %d ignored", found);
          }
          wavef = &ts->wave[ts->numwave-1];
        }
      }

      printf("orcaplugin) Wavefunction %s (%s):\n", action, wavef->info);
      printf("orcaplugin)   %d orbitals, %sexcitation %d, multiplicity %d\n",
             wavef->num_orbitals, spinstr, wavef->exci, wavef->mult);
    }
  }

  return i;
}


static int get_wavefunction(qmdata_t *data, qm_timestep_t *ts, qm_wavefunction_t *wf) {
  std::vector<float> orbitalEnergies;
  std::vector<int> orbitalOccupancies;
  std::vector<float> wavefunctionCoeffs;
  int num_orbitals = 0;

  char buffer[BUFSIZ];
  char word[8][BUFSIZ];
  long filepos;
  char *line;

  buffer[0] = '\0';
  int i = 0;
  for (i=0; i<8; i++) word[i][0] = '\0';

  if (wf == NULL) {
    PRINTERR;
    return FALSE;
  }

  wf->has_occup = FALSE;
  wf->has_orben = FALSE;
  wf->type = MOLFILE_WAVE_UNKNOWN;
  wf->spin = SPIN_ALPHA;
  wf->exci = 0;
  strncpy(wf->info, "unknown", MOLFILE_BUFSIZ);

  filepos = ftell(data->file);

  do {
    GET_LINE(buffer, data->file);
    line = trimleft(trimright(buffer));
    if(!strcmp(line, "Occupied Orbitals")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
    }
  } while(wf->type == MOLFILE_WAVE_UNKNOWN && strcmp(line, "== MOPAC DONE =="));

  if(wf->type == MOLFILE_WAVE_UNKNOWN) {
    #ifdef DEBUGGING_ORCA
        printf("mopacplugin) get_wavefunction(): No wavefunction found!\n");
    #endif
        fseek(data->file, filepos, SEEK_SET);
        return FALSE;
  } else {
    #ifdef DEBUGGING_ORCA
      printf("mopacplugin) Found wavefunction of type %d.\n", wf->type);
    #endif
  }

  eatline(data->file, 3);

  // number of read values from line;
  int numReadOrbitalIndices = 0;
  int numReadEnergies = 0;
  int numReadOccupancies = 0;
  int numReadCoefficients = 0;
  float numberOfElectrons = 0;
  float occupiedOrbitals = 0;
  std::vector<int> numberContractedBf;
  std::vector<int> wfAngMoment;
  std::vector<std::string> orbitalNames;
  MoCoeff allCoefficients;

  // orbital indices
  int n[8];
  int wavefunctionRead = 0;
  int firstRead = 1;
  int haveAngMom = 0;
  while(!wavefunctionRead) {
    float coeff[8], energies[8];
    float occ[8];
    char dumpName[BUFSIZ];
    char dumpBasisFunc[BUFSIZ];
    filepos = ftell(data->file);

    // reads the orbital indices
    if (firstRead == 1) {
      GET_LINE(buffer, data->file);
      firstRead++;
    }
    numReadOrbitalIndices = sscanf(buffer, "%s %s %d %d %d %d %d %d %d %d", &dumpName, &dumpBasisFunc, &n[0], &n[1], &n[2], &n[3], &n[4], &n[5], &n[6], &n[7]);
    if (!numReadOrbitalIndices || numReadOrbitalIndices == -1) {
      /* If there are no orbital indexes then this must be the
       * end of the wavefunction coefficient table. */
      fseek(data->file, filepos, SEEK_SET);
      wavefunctionRead = 1;
      break;
    }
    eatline(data->file, 1);
    numReadOrbitalIndices -= 2; // because of the "Root No." string in the beginning

    // reads the orbital energies
    GET_LINE(buffer, data->file);
    numReadEnergies = sscanf(buffer, "%f %f %f %f %f %f %f %f", &energies[0], &energies[1], &energies[2],
     &energies[3], &energies[4], &energies[5], &energies[6], &energies[7]);
    if (numReadEnergies != numReadOrbitalIndices) {
      printf("mopacplugin) Molecular Orbital section corrupted! energies.\n");
      break;
    }

    // store the energies in vector
    if (numReadEnergies != -1) {
      for (size_t c = 0; c < numReadEnergies; c++) {
        orbitalEnergies.push_back(energies[c]);
        std::cout << "Energy: " <<energies[c]<< std::endl;
        num_orbitals++;

        orbitalOccupancies.push_back(2);
        // std::cout << "Occupancy: " << occ[c] << std::endl;
	      numberOfElectrons += 2;
	      occupiedOrbitals++;
      }
    }

    // skip --- line
    eatline(data->file, 2);

    std::vector<std::vector<float>> moCoefficients;
    // we expect as many coefficients as numReadOccupancies, numReadEnergies, numReadOrbitalIndices!
    // read them as long as we find new coefficients
    int readingBlock = 1;
    int coefficientNumber = 0;
    int atomIndex = 0;
    int blockNumberOfContracted = 0;
    while(readingBlock) {
      GET_LINE(buffer, data->file);
      numReadCoefficients = sscanf(buffer, "%s %s %d %f %f %f %f %f %f %f %f", &dumpBasisFunc, &dumpName, &atomIndex,
      &coeff[0], &coeff[1], &coeff[2],&coeff[3], &coeff[4], &coeff[5], &coeff[6], &coeff[7]);
      // the coefficient number is the number of read elements minus 3 bc. of the atom and bf name
      coefficientNumber = (numReadCoefficients - 3);
      if (coefficientNumber == numReadOrbitalIndices) {
        // std::cout << "found coeffs: " << dumpName << "," << dumpBasisFunc << std::endl;
        if (firstRead == 2 && !haveAngMom) {
          std::string bfn = dumpBasisFunc;
          std::size_t found = bfn.find_first_not_of("-0123456789 ");
          if (found!=std::string::npos) {
            std::string orbital =  bfn.substr(found);
            orbitalNames.push_back(orbital);
            // std::cout << orbital << std::endl;
          } else {
            printf("mopacplugin) Could not determine orbital description.\n");
            return FALSE;
          }
        }
        // reading coefficients
        std::vector<float> currentMoCoeffs;
        for (size_t cidx = 0; cidx < coefficientNumber; cidx++) {
          currentMoCoeffs.push_back(coeff[cidx]);
          // std::cout << coeff[cidx] << std::endl;
        }
        moCoefficients.push_back(currentMoCoeffs);
      } else {
        // block seems to be finished
        readingBlock = 0;
        haveAngMom = 1;
        eatline(data->file, 1);
        GET_LINE(buffer, data->file);
        blockNumberOfContracted++;
        numberContractedBf.push_back(blockNumberOfContracted);
        blockNumberOfContracted = 0;
      }
      blockNumberOfContracted++;
    }
    allCoefficients.push_back(moCoefficients);
  }

  if ( std::adjacent_find( numberContractedBf.begin(), numberContractedBf.end(), std::not_equal_to<int>() ) != numberContractedBf.end() ) {
    printf("mopacplugin) Molecular orbital section corrupted. Did not read consistent number of contracted basis functions!\n");
    return FALSE;
  }
  printf("%d\n", num_orbitals);

  std::cout << numberContractedBf[0] << std::endl;

  // assign the number of contracted functions to wavefunction size
  data->wavef_size = numberContractedBf[0];
  wf->num_orbitals  = num_orbitals;
  wf->orb_energies = (float *) calloc(num_orbitals, sizeof(float));
  wf->orb_occupancies = (float *) calloc(num_orbitals, sizeof(float));
  wf->wave_coeffs = (float *) calloc(num_orbitals * data->wavef_size, sizeof(float));
  wf->has_occup = TRUE;
  wf->has_orben = TRUE;

  int cnt = 0;
  for (auto en : orbitalEnergies) {
    wf->orb_energies[cnt] = en;
    cnt++;
  }
  cnt = 0;
  for (auto occ : orbitalOccupancies) {
    wf->orb_occupancies[cnt] = occ;
    cnt++;
  }

  int rowIndex = 0, columnIndex = 0;
  int blockIdx = 0;
  int moBlockSize = 0;
  int moBlockIdx = 0;
  for (auto moBlock : allCoefficients) {
    for (auto moRow : moBlock) {
      // std::cout << rowIndex << std::endl;
      for (auto moCo : moRow) {
        if ((columnIndex * data->wavef_size + rowIndex) > num_orbitals * data->wavef_size) {
          std::cout << "something went wrong:" << columnIndex << std::endl;
          std::cout << "something went wrong:" << (columnIndex * data->wavef_size + rowIndex) << " vs. " << num_orbitals * data->wavef_size << std::endl;
          return FALSE;
        }
        // std::cout << orbitalNames[rowIndex] << std::endl;
        wf->wave_coeffs[columnIndex * data->wavef_size + rowIndex] = moCo;
        columnIndex++;
      }
      columnIndex = moBlockSize;
      rowIndex++;
    }
    rowIndex = 0;
    // 0-based!!!
    moBlockSize += moBlock[moBlockIdx].size();
    columnIndex = moBlockSize;
    moBlockIdx++;
    // std::cout << "bs: " << moBlockSize << std::endl;
  }

  // LOGGING for MO coefficients

  // float coeff2 = 0;
  // for (size_t t = 0; t < (num_orbitals * data->wavef_size); t++) {
  //   if (t % data->wavef_size == 0) {
  //     std::cout << "---------- " << t/num_orbitals << " c2: " << coeff2 << std::endl;
  //     coeff2 = 0;
  //   }
  //   coeff2 += wf->wave_coeffs[t]*wf->wave_coeffs[t];
  //   std::cout << wf->wave_coeffs[t] << std::endl;
  // }


  if (data->num_frames_read < 1) {
    data->angular_momentum = (int*)calloc(wfAngMoment.size(), sizeof(int));
  }

  // std::cout << "wfang: " << wfAngMoment.size() <<  " " << 3*data->wavef_size <<std::endl;
  for (size_t ang = 0; ang < wfAngMoment.size(); ang++) {
    data->angular_momentum[ang] = wfAngMoment[ang];
  }

  // TODO: This is just a workaround and might give wrong
  // results when reading unrestricted jobs
  data->num_occupied_A = occupiedOrbitals;
  data->num_occupied_B = occupiedOrbitals;
  data->num_electrons = numberOfElectrons;

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Total number of orbitals: " << num_orbitals << std::endl;
  std::cout << "Number of electrons in read wf: " << numberOfElectrons<< std::endl;
  std::cout << "Number of occupied orbitals: " << occupiedOrbitals << std::endl;
  std::cout << "Number of contracted bf: " << numberContractedBf[0] << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  return TRUE;
}

static int fill_basis_arrays(qmdata_t *data) {
  mopacdata *mopac = (mopacdata *)data->format_specific_data;
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
  std::cout<< "pcount: " << primcount << std::endl;

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
      printf("orcaplugin) WARNING: Couldn't find atomic number for basis set atom %s\n",
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
  printf("orcaplugin) Filled basis arrays.\n");

  return TRUE;
}

static int get_basis(qmdata_t *data) {
  data->basis_set = (basis_atom_t*)calloc(data->numatoms, sizeof(basis_atom_t));
  data->num_basis_atoms = data->numatoms;
  for (size_t atom = 0; atom < data->numatoms; atom++) {
    std::string el(data->atoms[atom].type);
    std::cout << "Atom " << atom << " : " << el << std::endl;
    int ns = 0;
    shell_t* tempShell = get_bas_4element(el, &ns);
    std::cout << "number of s: " << ns << std::endl;
    data->basis_set[atom].shell = tempShell;
    data->basis_set[atom].numshells = ns;
    data->num_shells += ns;
    strncpy(data->basis_set[atom].name, data->atoms[atom].type, 11);
  }
  int pc = 0;
  for(int i=0; i<data->num_basis_atoms; i++) {
    std::cout << "nbas atoms: " << i << std::endl;
    for (int j=0; j<data->basis_set[i].numshells; j++) {
      std::cout << "--- " << j << " " << data->basis_set[i].shell[j].numprims << std::endl;
      pc += data->basis_set[i].shell[j].numprims;
      if (data->basis_set[i].shell[j].numprims > 3) {
        std::cout << data->basis_set[i].name << std::endl;
        std::cout << data->basis_set[i].shell[j].type << std::endl;
        std::cout << data->basis_set[i].shell[j].prim[i].exponent << std::endl;
      }
    }
  }
  std::cout << "mypcount: " << pc << std::endl;
  return fill_basis_arrays(data);
}


/* Read the first trajectory frame. */
static int read_first_frame(qmdata_t *data) {
  /* Try reading the first frame.
   * If there is only one frame then also read the
   * final wavefunction. */
  if (!get_traj_frame(data, data->atoms, data->numatoms)) {
    return FALSE;
  }

  return TRUE;
}


static int parse_static_data(qmdata_t *data, int* natoms) {
  mopacdata *mopac = (mopacdata *)data->format_specific_data;
  if (!get_job_info(data)) return FALSE;

  if (!get_input_structure(data, mopac)) return FALSE;

  if (!get_basis(data)) return FALSE;

  if (!analyze_traj(data, mopac)) {
    printf("mopacplugin) WARNING: Truncated or abnormally terminated file!\n\n");
  }

  *natoms = data->numatoms;

  read_first_frame(data);

  print_input_data(data);
  std::cout << data->runtype << std::endl;

  return TRUE;
}

/* Analyze the trajectory.
 * Read the parameters controlling geometry search and
 * find the end of the trajectory, couinting the frames
 * on the way. Store the filepointer for the beginning of
 * each frame in *filepos_array. */
static int analyze_traj(qmdata_t *data, mopacdata *mopac) {
  char buffer[BUFSIZ], nserch[BUFSIZ];
  char *line;
  long filepos;
  filepos = ftell(data->file);

  data->filepos_array = (long* )calloc(1, sizeof(long ));

  /* currently, only one frame is supported!
   * lines 3130-3348 in gamessplugin.c
   */
  if (data->runtype == MOLFILE_RUNTYPE_ENERGY) {
    /* We have just one frame */
    data->num_frames = 1;
    pass_keyline(data->file, "1SCF WAS SPECIFIED", NULL);
    data->filepos_array[0] = ftell(data->file);

    /* Check wether SCF has converged */
    // if (pass_keyline(data->file,
    //                  "SCF IS UNCONVERGED, TOO MANY ITERATIONS",
    //                  "ENERGY COMPONENTS")==FOUND) {
    //   printf("mopacplugin) SCF IS UNCONVERGED, TOO MANY ITERATIONS\n");
    //   data->status = MOLFILE_QMSTATUS_SCF_NOT_CONV;
    // } else {
    //   data->status = MOLFILE_QMSTATUS_OPT_CONV;
    //   fseek(data->file, data->filepos_array[0], SEEK_SET);
    // }

    pass_keyline(data->file, "== MOPAC DONE ==", NULL);
    data->end_of_traj = ftell(data->file);

    /* Allocate memory for the frame */
    data->qm_timestep = (qm_timestep_t *)calloc(1, sizeof(qm_timestep_t));
    memset(data->qm_timestep, 0, sizeof(qm_timestep_t));

    return TRUE;
  } else if (data->runtype == MOLFILE_RUNTYPE_GRADIENT) {
    int appendedCalculations = 0;
    rewind(data->file);
    pass_keyline(data->file, "FINAL HEAT OF FORMATION", NULL); // TODO: adapt
    data->filepos_array[0] = ftell(data->file);
    data->num_frames = 1;

    while(TRUE) {
      if (!fgets(buffer, sizeof(buffer), data->file)) break;
      line = trimleft(buffer);

      std::string l(line);
      if (l.find("FINAL HEAT OF FORMATION") != std::string::npos && data->runtype==MOLFILE_RUNTYPE_GRADIENT) {
        appendedCalculations++;
        // std::cout << l << std::endl;
        if (data->num_frames > 0) {
          data->filepos_array = (long*)realloc(data->filepos_array, (data->num_frames+1)*sizeof(long));
        }
        data->filepos_array[data->num_frames] = ftell(data->file);
        data->num_frames++;
      }
    }

    if (appendedCalculations) {
      std::cout << "mopacplugin) Found multiple appended gradient calculations: " << data->num_frames << std::endl;
      pass_keyline(data->file, "== MOPAC DONE ==", NULL);
      data->end_of_traj = ftell(data->file);
      fseek(data->file, filepos, SEEK_SET);

      data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames,
                                                  sizeof(qm_timestep_t));
      memset(data->qm_timestep, 0, data->num_frames*sizeof(qm_timestep_t));
    } else {
      data->num_frames = 1;
      pass_keyline(data->file, "== MOPAC DONE ==", NULL);
      data->end_of_traj = ftell(data->file);

      /* Allocate memory for the frame */
      data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames, sizeof(qm_timestep_t));
      memset(data->qm_timestep, 0, sizeof(qm_timestep_t));
    }
    return TRUE;
  }
  else if (data->runtype == MOLFILE_RUNTYPE_OPTIMIZE) {
    std::cout << "mopacplugin) Reading trajectory of optimization not supported at the moment." << std::endl;
    return FALSE;
  }
  else {
    std::cout << "mopacplugin) Jobtype not supported for trajectory reading." << std::endl;
    return FALSE;
  }

  // printf("mopacplugin) Analyzing trajectory...\n");
  data->status = MOLFILE_QMSTATUS_UNKNOWN;

  while (TRUE) {
    if (!fgets(buffer, sizeof(buffer), data->file)) break;
    line = trimleft(buffer);

    std::string l(line);
    if (l.find("GEOMETRY OPTIMIZATION CYCLE") != std::string::npos && data->runtype==MOLFILE_RUNTYPE_OPTIMIZE) {
      // std::cout << l << std::endl;
      if (data->num_frames > 0) {
        data->filepos_array = (long*)realloc(data->filepos_array, (data->num_frames+1)*sizeof(long));
      }
      data->filepos_array[data->num_frames] = ftell(data->file);
      data->num_frames++;
    }
    else if (l.find("THE OPTIMIZATION HAS CONVERGED") != std::string::npos) {
      printf("mopacplugin) ==== End of trajectory (%d frames) ====\n", data->num_frames);
      data->status = MOLFILE_QMSTATUS_OPT_CONV;
    }
    else if (data->status == MOLFILE_QMSTATUS_OPT_CONV) {
      if(l.find("FINAL ENERGY EVALUATION AT THE STATIONARY POINT") != std::string::npos) {
        if (data->num_frames > 0) {
          data->filepos_array = (long*)realloc(data->filepos_array, (data->num_frames+1)*sizeof(long));
        }
        std::cout << "mopacplugin) found equilibrium geometry." << std::endl;
        data->filepos_array[data->num_frames] = ftell(data->file);
        data->num_frames++;
        goto_keyline(data->file, "TOTAL RUN TIME", NULL);
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
    // #ifdef DEBUGGING
    // printf("mopacplugin) atomicnum[%d] = %d\n", i, atom->atomicnumber);
    // #endif

    /* if (data->have_mulliken)
    atom->charge = data->qm_timestep->mulliken_charges[i];
    */
    cur_atom++;
  }

  return MOLFILE_SUCCESS;
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
  std::cout << "READ TIMESTEP" << std::endl;
  qmdata_t *data = (qmdata_t *)mydata;
  qm_timestep_t *cur_ts;
  int offset;
  int i = 0;
  int num_charge_sets = 0;

  if (data->trajectory_done == TRUE) {
    printf("mopacplugin) Trajectory done.\n");
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
  cur_ts = data->qm_timestep+data->num_frames_sent;
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
  if (cur_ts->have_mulliken) {
    offset = num_charge_sets*data->numatoms;
    for (i=0; i<data->numatoms; i++) {
      qm_ts->charges[offset+i] = cur_ts->mulliken_charges[i];
    }
    qm_ts->charge_types[num_charge_sets] = MOLFILE_QMCHARGE_MULLIKEN;
    num_charge_sets++;
  }
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
  /* store the wave function and orbital energies */
  if (cur_ts->wave) {
    std::cout << "mopacplugin) Have wavefunctions: " << cur_ts->numwave << " in frame: " << data->num_frames_sent << std::endl;
    for (i=0; i<cur_ts->numwave; i++) {
      std::cout << "Reading wf" << std::endl;
      qm_wavefunction_t *wave = &cur_ts->wave[i];
      qm_ts->wave[i].type         = wave->type;
      qm_ts->wave[i].spin         = wave->spin;
      qm_ts->wave[i].excitation   = wave->exci;
      qm_ts->wave[i].multiplicity = wave->mult;
      qm_ts->wave[i].energy       = wave->energy;
      strncpy(qm_ts->wave[i].info, wave->info, MOLFILE_BUFSIZ);

      if (wave->wave_coeffs) {
        memcpy(qm_ts->wave[i].wave_coeffs, wave->wave_coeffs, wave->num_orbitals*data->wavef_size*sizeof(float));
      }
      std::cout << "wave_coeffs OK" << std::endl;
      if (wave->orb_energies) {
        memcpy(qm_ts->wave[i].orbital_energies, wave->orb_energies,
               wave->num_orbitals*sizeof(float));
      }
      std::cout << "orb energies OK" << std::endl;
      if (wave->has_occup) {
        memcpy(qm_ts->wave[i].occupancies, wave->orb_occupancies,
               wave->num_orbitals*sizeof(float));
      }
      std::cout << "occup OK" << std::endl;
    }
  }
  //
  if (data->runtype == MOLFILE_RUNTYPE_ENERGY ||
      data->runtype == MOLFILE_RUNTYPE_HESSIAN) {
    /* We have only a single point */
    data->trajectory_done = TRUE;
  }

  data->num_frames_sent++;
  std::cout << "Frames sent: " << data->num_frames_sent << std::endl;

  return MOLFILE_SUCCESS;
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
    printf("mopacplugin) Probing timestep %d\n", data->num_frames_read);

    have = get_traj_frame(data, data->atoms, data->numatoms);
  }

  if (have) {
    // std::cout << "have frame" << std::endl;
    int i;
    qm_timestep_t *cur_ts;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep+data->num_frames_sent;

    std::cout << "numwave: " << cur_ts->numwave << std::endl;

    for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<cur_ts->numwave); i++) {
      meta->num_orbitals_per_wavef[i] = cur_ts->wave[i].num_orbitals;
      meta->has_occup_per_wavef[i]    = cur_ts->wave[i].has_occup;
      meta->has_orben_per_wavef[i]    = cur_ts->wave[i].has_orben;
      // std::cout << "occ: " << cur_ts->wave[i].has_occup << std::endl;
      // std::cout << "energy: " << cur_ts->wave[i].has_orben << std::endl;
    }
    meta->wavef_size      = data->wavef_size;
    meta->num_wavef       = cur_ts->numwave;
    meta->num_scfiter     = cur_ts->num_scfiter;
    meta->num_charge_sets = cur_ts->have_mulliken +
      cur_ts->have_lowdin + cur_ts->have_esp;
    if (cur_ts->gradient) meta->has_gradient = TRUE;

  } else {
    std::cout << "not have frame" << std::endl;
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


#endif


/******************************************************
 *
 * this function extracts the trajectory information
 * from the output file
 *
 * *****************************************************/
static int get_traj_frame(qmdata_t *data, qm_atom_t *atoms,
                          int natoms) {
  mopacdata *mopac = (mopacdata *)data->format_specific_data;
  qm_timestep_t *cur_ts;
  char buffer[BUFSIZ];
  char word[BUFSIZ];
  int units;
  buffer[0] = '\0';
  word[0]   = '\0';

  printf("mopacplugin) Timestep %d:\n", data->num_frames_read);
  printf("mopacplugin) ============\n");

  cur_ts = data->qm_timestep + data->num_frames_read;

  // debugging the trajectory reading file positions
  // printf("nfread: %d \n", data->num_frames_read);
  if (!data->filepos_array) {
    printf("filepos array empty!!!\n");
    return FALSE;
  }

  fseek(data->file, data->filepos_array[data->num_frames_read], SEEK_SET);

  /*
  * distinguish between job types
  * at the moment, only Single Points will work
  * lines 2840 - 3122 in gamessplugin.c
   */

  // reading geometries...

  if ((data->runtype == MOLFILE_RUNTYPE_OPTIMIZE || data->runtype == MOLFILE_RUNTYPE_GRADIENT) && data->num_frames > 1) {
    if (goto_keyline(data->file, "CARTESIAN COORDINATES (ANGSTROEM)", NULL)) {
      GET_LINE(buffer, data->file);
      // thisline(data->file);
      // UNITS ARE ANGSTROEM
      // bohr = 0;
      // sscanf()
    } else {
      printf("mopacplugin) No cartesian coordinates in ANGSTROEM found.\n");
    }
    // skip the ---- line
    eatline(data->file, 1);
    if (!get_coordinates(data->file, &data->atoms, units, &natoms)) {
      printf("mopacplugin) Couldn't find coordinates for timestep %d\n", data->num_frames_read);
    }
  }

  // if (get_scfdata(data, cur_ts) == FALSE) {
  //   printf("mopacplugin) Couldn't find SCF iterations for timestep %d\n",
  //          data->num_frames_read);
  // }

  /* Try reading canonical alpha/beta wavefunction */
  check_add_wavefunctions(data, cur_ts);

  /* Read point charged */
  // if (!cur_ts->have_mulliken && get_population(data, cur_ts)) {
  //   printf("mopacplugin) Mulliken charges found\n");
  // }

  /* Read the energy gradients (= -forces on atoms) */
  // if (get_gradient(data, cur_ts)) {
  //   printf("mopacplugin) Energy gradient found.\n");
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
      std::cout << "mopacplugin) Finished trajectory." << std::endl;
    }

    /* Try to read final wavefunction and orbital energies
     * A preexisting canonical wavefunction for this timestep
     * with the same characteristics (spin, exci, info) will
     * be overwritten by the final wavefuntion if it has more
     * orbitals. */
    // check_add_wavefunctions(data, cur_ts);
  }

  data->num_frames_read++;
  #ifdef DEBUGGING
    std::cout << "mopacplugin) Frames read: " << data->num_frames_read << std::endl;
  #endif

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
  std::cout << plugin.abiversion << std::endl;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "mopac";
  plugin.prettyname = "Mopac";
  plugin.author = "Maximilian Scheurer";
  plugin.majorv = 0;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
  plugin.filename_extension = "mopac";
  plugin.open_file_read = open_mopac_read;
  plugin.read_structure = read_mopac_structure;
  plugin.close_file_read = close_mopac_read;
  //
  plugin.read_qm_metadata = read_mopac_metadata;
  plugin.read_qm_rundata  = read_mopac_rundata;

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


/*****************************************************
 *
 * provide VMD with the sizes of the QM related
 * data structure arrays that need to be made
 * available
 *
 *****************************************************/
static int read_mopac_metadata(void *mydata,
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
static int read_mopac_rundata(void *mydata,
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
    // std::cout << "free data->angular_momentum" << std::endl;
    // free(data->angular_momentum);
  }
#endif

  return MOLFILE_SUCCESS;
}

/**********************************************************
 *
 * close file and free memory
 *
 **********************************************************/
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
  printf(" ATOM      ATOMIC                      COORDINATES (Angstroem)\n");
  printf("           CHARGE         X                   Y                   Z\n");
  for (i=0; i<data->numatoms; i++) {
    printf(" %-8s %6d", data->atoms[i].type, data->atoms[i].atomicnum);

    printf("%17.10f",   data->atoms[i].x);
    printf("%20.10f",   data->atoms[i].y);
    printf("%20.10f\n", data->atoms[i].z);
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
