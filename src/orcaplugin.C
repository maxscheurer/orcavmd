/*
Orca VMD plugin
Authors: Maximilian Scheurer, Marcelo Melo, April 2017
Contributions from: Michael F. Herbst
Inspired from gamessplugin.c
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
#include "Matrix.h"

typedef std::vector<std::vector<std::vector<float>>> MoCoeff;

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

typedef std::vector<std::vector<float>> CoeffRowBlock;

// Matrix for conversion of pure coefficients to cartesian coefficients, d- and f-shell
Matrix *convD = new Matrix("{{ 0, 1, 0, 0, 0},"
                            "{ 0, 0, 1, 0, 0},"
                            "{ 0, 0, 0, 0, 1},"
                            "{-1, 0, 0, 1, 0},"
                            "{-1, 0, 0,-1, 0},"
                            "{ 2, 0, 0, 0, 0}}");

Matrix *convF = new Matrix("{{ 0, 0, 0, 0, 1, 0, 0},"
                            "{ 0, 0,-1, 0, 0, 0, 3},"
                            "{-3, 0, 0, 1, 0, 0, 0},"
                            "{ 0,-1, 0, 0, 0,-3, 0},"
                            "{-3, 0, 0,-1, 0, 0, 0},"
                            "{ 0, 4, 0, 0, 0, 0, 0},"
                            "{ 0, 0, 4, 0, 0, 0, 0},"
                            "{ 0,-1, 0, 0, 0, 1, 0},"
                            "{ 0, 0,-1, 0, 0, 0,-1},"
                            "{ 2, 0, 0, 0, 0, 0, 0}}");


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


// collect information about the jobtype
static int get_job_info(qmdata_t *data);

// analyze trajectory, getting number of frames, file positions etc.
static int analyze_traj(qmdata_t *data, orcadata *orca);

// atom definitions & geometry
static int get_input_structure(qmdata_t *data, orcadata *orca);

// reading coord block
static int get_coordinates(FILE *file, qm_atom_t **atoms, int unit, int *numatoms);

// reading first trajectory frame
static int read_first_frame(qmdata_t *data);

int get_basis(qmdata_t *data);

static int get_wavefunction(qmdata_t *data, qm_timestep_t *ts, qm_wavefunction_t *wf);

// main routine for extracting "step"-dependent info from the output file
static int get_traj_frame(qmdata_t *data, qm_atom_t *atoms, int natoms);

static int get_scfdata(qmdata_t *data, qm_timestep_t *ts);

static int check_add_wavefunctions(qmdata_t *data, qm_timestep_t *ts);

static int fill_basis_arrays(qmdata_t *data);

static int shelltype_int(char type);

static CoeffRowBlock convertPure(CoeffRowBlock pureBlock);

// for VMD
static int read_orca_metadata(void *mydata, molfile_qm_metadata_t *metadata);
static int read_orca_rundata(void *mydata, molfile_qm_t *qm_data);
static int read_qm_timestep_metadata(void *mydata, molfile_qm_timestep_metadata_t *meta);

static int read_timestep(void *mydata, int natoms,
       molfile_timestep_t *ts, molfile_qm_metadata_t *qm_metadata,
			 molfile_qm_timestep_t *qm_ts);

static int read_timestep_metadata(void *mydata, molfile_timestep_metadata_t *meta);

static void print_input_data(qmdata_t *data);

/*
String manipulation
*/
template<typename Out>
void split(const std::string &s, char delim, Out result);

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
                   const std::string& whitespace = " \t")
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

std::vector<std::string> semiempiricals({"MNDO","PM3","AM1"});

static int get_job_info(qmdata_t *data) {
	long filepos;
	char buffer[BUFSIZ];
	filepos = ftell(data->file);
	if (goto_keyline(data->file, "INPUT FILE", NULL)) {
		eatline(data->file, 3);
	}
	std::vector<std::string> inputFile;
	int endOfInput = 0;
	while(!endOfInput) {
		GET_LINE(buffer, data->file);
		std::string test(buffer);
		auto ln = test.find_first_of("123456789");
		auto beforeLine = test.find_first_of(">");
		int lineNumber = stoi(test.substr(ln,beforeLine-ln));
		std::string lContent = trim(test.substr(beforeLine+1));
		inputFile.push_back(lContent);
		if (lineNumber == 1) {
			if (lContent[0] != '!') {
				std::cout << "orcaplugin) File corrupted in"
					" input section. Exiting" << std::endl;
				return FALSE;
			} else {
				std::cout << "Found commands." << std::endl;
			}
			for (auto method : semiempiricals) {
				if (lContent.find(method) != std::string::npos) {
					const char *m = method.c_str();
					strncpy(data->gbasis, m, sizeof(char)*strlen(m));
					strncpy(data->basis_string, "STO-3G", sizeof(char)*strlen("STO-3G"));
					std::cout << "orcaplugin) semiempirical Method used." << std::endl;
					break;
				}
			}
		}
		if (test.find("END OF INPUT") !=std::string::npos) {
			endOfInput = 1;
		}
	}

  data->runtype = MOLFILE_RUNTYPE_ENERGY;
  std::string lower = inputFile[0];
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  if (lower.find("opt") != std::string::npos) {
    data->runtype = MOLFILE_RUNTYPE_OPTIMIZE;
    std::cout << "orcaplugin) Optimization loaded." << std::endl;
  } else if (lower.find("engrad") != std::string::npos) {
    data->runtype = MOLFILE_RUNTYPE_GRADIENT;
  }

  int totalCharge;
  if (goto_keyline(data->file, "Total Charge", NULL)) {
    GET_LINE(buffer,data->file);
		std::string chargeLine(buffer);
    std::vector<std::string> chargeVec = split(reduce(chargeLine), ' ');
    totalCharge = stoi(*(chargeVec.end()-1));
    std::cout << "orcaplugin) Found molecule charge: " << totalCharge << std::endl;
    data->totalcharge = totalCharge;
  } else {
    std::cout << "orcaplugin) No molecule charge found. Exiting" << std::endl;
    return FALSE;
  }

  int multiplicity;
  if (goto_keyline(data->file, "Multiplicity", NULL)) {
    GET_LINE(buffer,data->file);
		std::string multLine(buffer);
    std::vector<std::string> multVec = split(reduce(multLine), ' ');
    multiplicity = stoi(*(multVec.end()-1));
    std::cout << "orcaplugin) Found molecule multiplicity: " << multiplicity << std::endl;
    data->multiplicity = multiplicity;
  } else {
    std::cout << "orcaplugin) No molecule multiplicity found. Exiting" << std::endl;
    return FALSE;
  }

  int nEl;
  if (goto_keyline(data->file, "Number of Electrons", NULL)) {
    GET_LINE(buffer,data->file);
		std::string nElLine(buffer);
    std::vector<std::string> nElVec = split(reduce(nElLine), ' ');
    nEl = stoi(*(nElVec.end()-1));
    std::cout << "orcaplugin) Found number of electrons: " << nEl << std::endl;
    data->num_electrons = nEl;
  } else {
    std::cout << "orcaplugin) Number of electrons not found. Exiting" << std::endl;
    return FALSE;
  }

  rewind(data->file);

	return TRUE;
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
  if (!get_job_info(data)) return FALSE;

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
  int semiempirical = 0;

  if (!strcmp(data->gbasis, "MNDO") ||
      !strcmp(data->gbasis, "AM1")  ||
      !strcmp(data->gbasis, "PM3")) {
    semiempirical = 1;
  }

  /* Search for "ATOMIC BASIS SET" line */
  if (pass_keyline(data->file, "BASIS SET IN INPUT FORMAT", NULL) != FOUND and
		  pass_keyline(data->file, "GAUSSIAN BASIS SET", NULL) != FOUND) {
    printf("orcaplugin) No basis set found!\n");
    return FALSE;
  }

  /* initialize buffers */
  buffer[0] = '\0';
  for (i=0; i<3; i++) word[i][0] = '\0';


  /* skip the next 3 lines */
  // eatline(data->file, 2);

  /* Allocate space for the basis for all atoms */
  data->basis_set = (basis_atom_t*)calloc(data->numatoms, sizeof(basis_atom_t));

  filepos = ftell(data->file);
  i = 0; /* basis atom counter */
  char elementName[11];
  int finished = FALSE;

  basis_atom_t* tempBasis = (basis_atom_t*)calloc(1, sizeof(basis_atom_t));

  prim_t *prim;

  std::vector<std::string> readElements;
  for (size_t atom = 0; atom < data->numatoms; atom++) {


  }

  while (!finished && !semiempirical) {
    // printf("Trying to read bf. \n");
    if (pass_keyline(data->file, "Basis set for element", NULL) == FOUND ) {
      GET_LINE(buffer, data->file);
      numread = sscanf(buffer,"%s %s",&word[0][0], &word[1][0]);
      strcpy(elementName, &word[1][0]);
      printf("New element found: %s\n", &word[1][0]);
      std::string atype(&word[1][0]);
      if (std::find(readElements.begin(), readElements.end(), atype) != readElements.end()) {
        break;
      } else {
        readElements.push_back(atype);
      }
      int elementCompleted = 0;


      shell = (shell_t*)calloc(1, sizeof(shell_t));
      numshells = 0;
      int readingShell = 0;
      int primcounter = 0;

      // this is very sloppy at the moment...
      // for PM3 etc. Orca prints the bf per atom...

      while(!elementCompleted) {
        GET_LINE(buffer, data->file);
        numread = sscanf(buffer,"%s %s %s",&word[0][0], &word[1][0],&word[2][0]);
        // printf("numread: %d -- %s %s %s \n",numread, &word[0][0], &word[1][0],&word[2][0]);
        switch (numread) {
          case 1:
            if (strcmp(trimleft(trimright(&word[0][0])), "end")) {
              // printf("Section ended. \n");
              elementCompleted = 1;
              break;
            }
          case 2:
            shell[numshells].numprims = atoi(trimleft(trimright(&word[1][0])));
            shell[numshells].type = shelltype_int(word[0][0]);
            // printf("orcaplugin) Type: %d NPrims: %d\n", shell[numshells].type, shell[numshells].numprims);
            primcounter = 0;
            prim = (prim_t*)calloc(shell[numshells].numprims, sizeof(prim_t));
            // printf("!! Address of prim is %p\n", (void *)prim);
            shell[numshells].prim = prim;
            numshells++;
            if (numshells) {
              shell = (shell_t*)realloc(shell, (numshells+1)*sizeof(shell_t));
            }
            break;
          case 3:
            prim[primcounter].exponent = atof(&word[1][0]);
            prim[primcounter].contraction_coeff = atof(&word[2][0]);
            // printf("%f - %f\n", prim[primcounter].exponent, prim[primcounter].contraction_coeff);
            primcounter++;
            break;
          default:
            printf("orcaplugin) Unknown line in basis functions.\n");
            elementCompleted = 1;
            break;
        }
      }
      // printf("Number of shells: %d \n", numshells);
      strcpy(tempBasis[i].name, elementName);
      tempBasis[i].numshells = numshells;
      tempBasis[i].shell = shell;
      i++;
      if (i) {
        tempBasis = (basis_atom_t*)realloc(tempBasis, (i+1)*sizeof(basis_atom_t));
      }
      // set prim to nullpointer!
      prim = NULL;
    } else {
      finished = TRUE;
      printf("orcaplugin) Reading basis set finished! \n");
    }
  }

  finished = 0;
  // semiempirical sto-3g
  int ngauss = 3;
  int atomCounter = 0;
  int shellNumber;
  while(!finished && semiempirical) {
  	if(goto_keyline(data->file,"shells",NULL) == FOUND && atomCounter <(data->numatoms-1)) {
  		// thisline(data->file);
  		GET_LINE(buffer, data->file);
  		std::string lineString(buffer);
  		std::vector<std::string> elements = split(reduce(lineString), ' ');
  		shellNumber = stoi(elements[2]);
      // std::cout << "shell number: " << shellNumber << std::endl;
  		shell = (shell_t*) calloc(shellNumber, sizeof(shell_t));
  		data->basis_set[atomCounter].shell = shell;
  		data->basis_set[atomCounter].numshells = shellNumber;
  		for (size_t shellIdx = 0; shellIdx < shellNumber; ++shellIdx) {
  			GET_LINE(buffer, data->file);
  			prim = (prim_t*) calloc(3, sizeof(prim_t));
  			shell[shellIdx].prim = prim;
  			shell[shellIdx].numprims = 3;
        shell[shellIdx].type = shellIdx;
  			for (size_t nbas = 0; nbas < ngauss; ++nbas) {
  				GET_LINE(buffer, data->file);
  				std::string l(buffer);
  				std::vector<std::string> coeff = split(reduce(l), ' ');
  				prim[nbas].exponent = stod(coeff[0]);
  				prim[nbas].contraction_coeff = stod(coeff[1]);
  			}
  		}
  		data->num_shells += shellNumber;
  		data->num_basis_atoms++;
  		data->num_basis_funcs += 3;
      strncpy(data->basis_set[atomCounter].name, data->atoms[atomCounter].type, 11);
  		atomCounter++;
	  } else {
      finished = TRUE;
      prim = NULL;
      std::cout << "orcaplugin) Reading STO-3G basis for semiempirical method finished." << std::endl;
      // free unused stuff
      free(tempBasis);
      // We return here without further ado.
      return fill_basis_arrays(data);
    }
  }

  // As we read GTOs from the Orca output file, we need to
  // loop over all atoms and assign the basis functions

  // printf("orcaplugin) Parsed basis set of %d elements. \n", i);
  // for (size_t j = 0; j < i; j++) {
  //   printf("Element: %s\n", tempBasis[j].name);
  //   printf("- NShells: %d\n", tempBasis[j].numshells);
  //   for (size_t k = 0; k < tempBasis[j].numshells; k++) {
  //     printf("--- NPrims: %d \n", tempBasis[j].shell[k].numprims);
  //     for (size_t o = 0; o < tempBasis[j].shell[k].numprims; o++) {
  //       float expo = tempBasis[j].shell[k].prim[o].exponent;
  //       float cont = tempBasis[j].shell[k].prim[o].contraction_coeff;
  //       printf("----- E= %f , C= %f \n", expo, cont);
  //     }
  //   }
  // }

  // Allocate an array of zeros to store whether the tempBasis
  // of the same index was actually used.
  int* tempBasisUsed = (int*) calloc(i, sizeof(int));

  const char* currentElement;
  basis_atom_t* currentBasis;
  for (size_t n = 0; n < data->numatoms; n++) {
    currentElement = data->atoms[n].type;
    for (size_t j = 0; j < i; j++) {
      if (!strcmp(currentElement, tempBasis[j].name)) {
        // printf("orcaplugin) found basis for element %s\n", currentElement);
        currentBasis = &tempBasis[j];
        tempBasisUsed[j] = 1;
      }
    }
    // printf("orcaplugin) Basis for element %s has %d shells.\n", currentElement, currentBasis->numshells);
    data->basis_set[n].shell = (shell_t *) calloc(currentBasis->numshells, sizeof(shell_t));
    memcpy(data->basis_set[n].shell, currentBasis->shell, currentBasis->numshells * sizeof(shell_t));
    data->basis_set[n].numshells = currentBasis->numshells;
    data->num_shells += currentBasis->numshells;
    for (size_t p = 0; p < currentBasis->numshells; ++p) {
      data->basis_set[n].shell[p].prim = (prim_t *) calloc(currentBasis->shell[p].numprims, sizeof(prim_t));
      memcpy(data->basis_set[n].shell[p].prim, currentBasis->shell[p].prim, currentBasis->shell[p].numprims * sizeof(prim_t));
	    data->num_basis_funcs += currentBasis->shell[p].numprims;
    }
    data->num_basis_atoms++;
    strncpy(data->basis_set[n].name, currentBasis->name, 11);
    currentBasis = NULL;
    currentElement  = NULL;
  }

  // cleaning up memory for tempBasis and its shells
  for(size_t idx = 0; idx < i; ++idx) {
    for (size_t shellIdx = 0; shellIdx < tempBasis[idx].numshells; shellIdx++) {
      shell_t* cshell = &tempBasis[idx].shell[shellIdx];
      free(cshell->prim);
      cshell->prim = NULL;
    }
    free(tempBasis[idx].shell);
    tempBasis[idx].shell = NULL;
  }
  free(tempBasis);
  tempBasis = NULL;
  shell = NULL;
  printf("orcaplugin) Parsed %d uncontracted basis functions.\n", data->num_basis_funcs);
  free(tempBasisUsed);

  // /* allocate and populate flat arrays needed for molfileplugin */
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

  // data->angular_momentum = (int*)calloc(3*data->wavef_size, sizeof(int));

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

  if (data->runtype == MOLFILE_RUNTYPE_OPTIMIZE || data->runtype == MOLFILE_RUNTYPE_GRADIENT && data->num_frames > 1) {
    if (goto_keyline(data->file, "CARTESIAN COORDINATES (ANGSTROEM)", NULL)) {
      GET_LINE(buffer, data->file);
      // thisline(data->file);
      // UNITS ARE ANGSTROEM
      // bohr = 0;
      // sscanf()
    } else {
      printf("orcaplugin) No cartesian coordinates in ANGSTROEM found.\n");
    }
    // skip the ---- line
    eatline(data->file, 1);
    if (!get_coordinates(data->file, &data->atoms, units, &natoms)) {
      printf("orcaplugin) Couldn't find coordinates for timestep %d\n", data->num_frames_read);
    }
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


  /* Read the energy gradients (= -forces on atoms) */
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
      std::cout << "orcaplugin) Finished trajectory." << std::endl;
    }

    /* Try to read final wavefunction and orbital energies
     * A preexisting canonical wavefunction for this timestep
     * with the same characteristics (spin, exci, info) will
     * be overwritten by the final wavefuntion if it has more
     * orbitals. */
    check_add_wavefunctions(data, cur_ts);
  }

  data->num_frames_read++;
  std::cout << "orcaplugin) Frames read: " << data->num_frames_read << std::endl;

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
#ifdef DEBUGGING
      printf("orcaplugin) No canonical wavefunction present for timestep %d\n", data->num_frames_read);
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

std::vector<std::vector<int>> dAngMom{{1,0,1},{0,1,1},{1,1,0},{2,0,0},{0,2,0},{0,0,2}};
std::vector<std::vector<int>> fAngMom{{1,1,1},{2,1,0},{2,0,1},{1,2,0},{0,2,1},{1,0,2},
{0,1,2},{3,0,0},{0,3,0},{0,0,3}};

typedef enum {
  dshell, fshell
} Shelltype;

static CoeffRowBlock convertPure(CoeffRowBlock pureBlock) {
  Shelltype stype;
  CoeffRowBlock resultBlock;
  switch (pureBlock.size()) {
    case 5:
      stype = dshell;
      break;
    case 7:
      stype = fshell;
      break;
    default:
      return resultBlock;
  }

  Matrix *m = new Matrix(pureBlock);
  Matrix *multiplication;
  if (stype == dshell) {
    multiplication = Matrix::multiply(convD, m);
  } else if (stype == fshell) {
    multiplication = Matrix::multiply(convF, m);
  } else {
    return resultBlock;
  }
  // m->printMatrix();
  // multiplication->printMatrix();
  resultBlock = multiplication->toVector();
  delete multiplication;
  delete m;

  // OLD CODE!
  // first get the unchanged d-orbitals xz, yz, xy
  // std::vector<std::string> unchangedOrbitals{"xz", "yz", "xy"};
  // std::vector<int> unchangedIndices{1, 2, 4};
  // for (auto idx : unchangedIndices) {
  //   resultBlock.push_back(pureBlock[idx]);
  // }
  //
  // int alpha = 3;
  // int beta = 0;
  // // now convert the others
  //
  // int rows = pureBlock[0].size(); // should be 5 in this case
  //
  // std::vector<float> x2;
  // for (size_t i = 0; i < pureBlock[0].size(); i++) {
  //   x2.push_back(pureBlock[alpha][i]-pureBlock[beta][i]);
  // }
  //
  // std::vector<float> y2;
  // for (size_t i = 0; i < pureBlock[0].size(); i++) {
  //   y2.push_back(-pureBlock[beta][i]-pureBlock[alpha][i]);
  // }
  //
  // std::vector<float> z2;
  // for (size_t i = 0; i < pureBlock[0].size(); i++) {
  //   z2.push_back(2*pureBlock[beta][i]);
  // }
  //
  // resultBlock.push_back(x2);
  // resultBlock.push_back(y2);
  // resultBlock.push_back(z2);

  return resultBlock;
}

static int get_wavefunction(qmdata_t *data, qm_timestep_t *ts, qm_wavefunction_t *wf) {
  std::vector<float> orbitalEnergies;
  std::vector<int> orbitalOccupancies;
  std::vector<float> wavefunctionCoeffs;
  int num_orbitals = 0;

  char buffer[BUFSIZ];
  char word[6][BUFSIZ];
  long filepos;
  char *line;

  buffer[0] = '\0';
  int i = 0;
  for (i=0; i<6; i++) word[i][0] = '\0';

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
    if(!strcmp(line, "MOLECULAR ORBITALS")) {
      wf->type = MOLFILE_WAVE_CANON;
      strncpy(wf->info, "canonical", MOLFILE_BUFSIZ);
    }
  } while(wf->type == MOLFILE_WAVE_UNKNOWN && strcmp(line, "FINAL SINGLE POINT ENERGY"));

  if(wf->type == MOLFILE_WAVE_UNKNOWN) {
    #ifdef DEBUGGING
        printf("orcaplugin) get_wavefunction(): No wavefunction found!\n");
    #endif
        fseek(data->file, filepos, SEEK_SET);
        return FALSE;
  } else {
    #ifdef DEBUGGING
      printf("orcaplugin) Found wavefunction of type %d.\n", wf->type);
    #endif
  }

  eatline(data->file, 1);

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
  int n[6];
  int wavefunctionRead = 0;
  int firstRead = 1;
  int haveAngMom = 0;
  while(!wavefunctionRead) {
    float coeff[6], energies[6];
    float occ[6];
    char dumpName[BUFSIZ];
    char dumpBasisFunc[BUFSIZ];
    filepos = ftell(data->file);

    // reads the orbital indices
    if (firstRead == 1) {
      GET_LINE(buffer, data->file);
      firstRead++;
    }

    numReadOrbitalIndices = sscanf(buffer, "%d %d %d %d %d %d", &n[0], &n[1], &n[2], &n[3], &n[4], &n[5]);
    if (!numReadOrbitalIndices || numReadOrbitalIndices == -1) {
      /* If there are no orbital indexes then this must be the
       * end of the wavefunction coefficient table. */
      fseek(data->file, filepos, SEEK_SET);
      wavefunctionRead = 1;
      break;
    }

    // reads the orbital energies
    GET_LINE(buffer, data->file);
    numReadEnergies = sscanf(buffer, "%f %f %f %f %f %f", &energies[0], &energies[1], &energies[2],
     &energies[3], &energies[4], &energies[5]);
    if (numReadEnergies != numReadOrbitalIndices) {
      printf("orcaplugin) Molecular Orbital section corrupted! energies.\n");
      break;
    }

    // store the energies in vector
    if (numReadEnergies != -1) {
      for (size_t c = 0; c < numReadEnergies; c++) {
        orbitalEnergies.push_back(energies[c]);
        // std::cout << "Energy: " <<energies[c]<< std::endl;
      }
    }
    // reads the orbital occupancies
    GET_LINE(buffer, data->file);
    numReadOccupancies = sscanf(buffer, "%f %f %f %f %f %f", &occ[0], &occ[1], &occ[2],
      &occ[3], &occ[4], &occ[5]);
    if (numReadOccupancies != numReadOrbitalIndices) {
      printf("orcaplugin) Molecular Orbital section corrupted! %d\n", numReadOccupancies);
      break;
    }

    // stores the occupancies in vector
    if(numReadOccupancies) {
      for (size_t c = 0; c < numReadOccupancies; c++) {
        orbitalOccupancies.push_back((int)occ[c]);
        // std::cout << "Occupancy: " << occ[c] << std::endl;
	numberOfElectrons += occ[c];
	if (occ[c]) {
	  occupiedOrbitals++;
	}
      }
      num_orbitals += numReadOccupancies;
    }

    // skip --- line
    eatline(data->file, 1);

    std::vector<std::vector<float>> moCoefficients;
    // we expect as many coefficients as numReadOccupancies, numReadEnergies, numReadOrbitalIndices!
    // read them as long as we find new coefficients
    int readingBlock = 1;
    int coefficientNumber = 0;
    while(readingBlock) {
      GET_LINE(buffer, data->file);
      numReadCoefficients = sscanf(buffer, "%s %s %f %f %f %f %f %f", &dumpName, &dumpBasisFunc,
      &coeff[0], &coeff[1], &coeff[2],&coeff[3], &coeff[4], &coeff[5]);
      // the coefficient number is the number of read elements minus 2 bc. of the atom and bf name
      coefficientNumber = (numReadCoefficients - 2);
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
            printf("orcaplugin) Could not determine orbital description.\n");
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
      }
    }
    allCoefficients.push_back(moCoefficients);
  }

  if ( std::adjacent_find( numberContractedBf.begin(), numberContractedBf.end(), std::not_equal_to<int>() ) != numberContractedBf.end() ) {
    printf("orcaplugin) Molecular orbital section corrupted. Did not read consistent number of contracted basis functions!\n");
    return FALSE;
  }

  // now loop over MO blocks and convert coefficients if needed
  int readingPureFunction = 0;
  std::vector<std::vector<float>> pureFunction;
  std::vector<std::string> pureFunctionName;
  int expectedNumberOfPureFunctions = 0;
  std::vector<std::string> orbList{"s","px","py","pz","dz2","dxz","dyz","dx2y2","dxy",
  "f0","f+1","f-1","f+2","f-2","f+3","f-3"};
  int orbRowIndex = 0;
  int blockIdx = 0;

  // for (auto name : orbitalNames) {
  //   std::cout << name << std::endl;
  // }

  MoCoeff newAllCoefficients;

  for (auto moBlock : allCoefficients) {
    CoeffRowBlock newRows;
    int blockNumberOfContracted = 0;
    for (auto moRow : moBlock) {
      int angX = 0, angY = 0, angZ = 0;
      std::string orbital = orbitalNames[orbRowIndex];
      std::vector<std::string>::iterator orbIndex = find(std::begin(orbList), std::end(orbList), orbital);
      if (std::end(orbList) != orbIndex) {
        int listIndex = orbIndex-orbList.begin();
        switch (listIndex) {
          case 0:
            break;
          case 1:
            angX = 1;
            break;
          case 2:
            angY = 1;
            break;
          case 3:
            angZ = 1;
            break;
          default:
            break;
        }
        // d-shell found
        // std::cout << "orbital index: " << orbIndex-orbList.begin() << " " << orbital << std::endl;
        if (listIndex > 3) {
          if (!readingPureFunction) {
            pureFunction.clear();
            pureFunctionName.clear();
            if (orbital.compare(0, 1, "d") == 0) {
              expectedNumberOfPureFunctions = 5;
              // std::cout << "Trying to read d-functions." << std::endl;
            } else if (orbital.compare(0, 1, "f") == 0) {
              expectedNumberOfPureFunctions = 7;
              // std::cout << "Trying to read f-functions." << std::endl;
            }
            readingPureFunction = 1;
            pureFunction.push_back(moRow);
            pureFunctionName.push_back(orbital);
          } else {
            pureFunction.push_back(moRow);
            pureFunctionName.push_back(orbital);
            if (pureFunction.size() == expectedNumberOfPureFunctions) {
              // std::cout << "found complete pure function set." << std::endl;
              CoeffRowBlock newBlock = convertPure(pureFunction);
              blockNumberOfContracted+= newBlock.size();
              for (auto r : newBlock) {
                newRows.push_back(r);
              }
              if (!blockIdx) {
                for (size_t i = 0; i < newBlock.size(); i++) {
                  if (expectedNumberOfPureFunctions == 5) {
                    for (auto angMom : dAngMom[i]) {
                      wfAngMoment.push_back(angMom);
                    }
                  } else if (expectedNumberOfPureFunctions == 7) {
                    for (auto angMom : fAngMom[i]) {
                      wfAngMoment.push_back(angMom);
                    }
                  }
                }
              }
              readingPureFunction = 0;
            }
          }
        } else {
          newRows.push_back(moRow);
          blockNumberOfContracted++;
          if (!blockIdx) {
            wfAngMoment.push_back(angX);
            wfAngMoment.push_back(angY);
            wfAngMoment.push_back(angZ);
          }
        }
      } else {
        std::cout << "orcaplugin) ERROR. Only s/p/d/f-shells are supported." << std::endl;
        return FALSE;
      }

      orbRowIndex++;
    }
    numberContractedBf.push_back(blockNumberOfContracted);
    newAllCoefficients.push_back(newRows);
    orbRowIndex = 0;
    blockIdx++;
  }


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
  blockIdx = 0;
  int moBlockSize = 0;
  int moBlockIdx = 0;
  for (auto moBlock : newAllCoefficients) {
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
  /*
  float coeff2 = 0;
  for (size_t t = 0; t < (num_orbitals * data->wavef_size); t++) {
    if (t % data->wavef_size == 0) {
      std::cout << "---------- " << t/num_orbitals << " c2: " << coeff2 << std::endl;
      coeff2 = 0;
    }
    coeff2 += wf->wave_coeffs[t]*wf->wave_coeffs[t];
    std::cout << wf->wave_coeffs[t] << std::endl;
  }
  */

  if (data->num_frames_read < 1) {
    data->angular_momentum = (int*)calloc(wfAngMoment.size(), sizeof(int));
  }

  // std::cout << "wfang: " << wfAngMoment.size() <<  " " << 3*data->wavef_size <<std::endl;
  for (size_t ang = 0; ang < wfAngMoment.size(); ang++) {
    data->angular_momentum[ang] = wfAngMoment[ang];
  }

  // TODO: REMOVE!!!
  // hardcoded for TEST!!!
  data->num_occupied_A = occupiedOrbitals;
  data->num_occupied_B = occupiedOrbitals;
  // data->num_electrons = numberOfElectrons;

  // if (data->num_electrons != numberOfElectrons) {
  //
  // }

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Total number of orbitals: " << num_orbitals << std::endl;
  std::cout << "Number of electrons in read wf: " << numberOfElectrons<< std::endl;
  std::cout << "Number of occupied orbitals: " << occupiedOrbitals << std::endl;
  std::cout << "Number of contracted bf: " << numberContractedBf[0] << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  return TRUE;
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
  if (data->runtype == MOLFILE_RUNTYPE_ENERGY) {
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
  } else if (data->runtype == MOLFILE_RUNTYPE_GRADIENT) {
    int appendedCalculations = 0;
    rewind(data->file);
    pass_keyline(data->file, "Energy+Gradient Calculation", NULL);
    data->filepos_array[0] = ftell(data->file);
    data->num_frames = 1;

    while(TRUE) {
      if (!fgets(buffer, sizeof(buffer), data->file)) break;
      line = trimleft(buffer);

      std::string l(line);
      if (l.find("Energy+Gradient Calculation") != std::string::npos && data->runtype==MOLFILE_RUNTYPE_GRADIENT) {
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
      std::cout << "orcaplugin) Found multiple appended gradient calculations: " << data->num_frames << std::endl;
      pass_keyline(data->file, "FINAL SINGLE POINT ENERGY", NULL);
      data->end_of_traj = ftell(data->file);
      fseek(data->file, filepos, SEEK_SET);

      data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames,
                                                  sizeof(qm_timestep_t));
      memset(data->qm_timestep, 0, data->num_frames*sizeof(qm_timestep_t));
    } else {
      data->num_frames = 1;
      pass_keyline(data->file, "FINAL SINGLE POINT ENERGY", NULL);
      data->end_of_traj = ftell(data->file);

      /* Allocate memory for the frame */
      data->qm_timestep = (qm_timestep_t *)calloc(data->num_frames, sizeof(qm_timestep_t));
      memset(data->qm_timestep, 0, sizeof(qm_timestep_t));
    }
    return TRUE;
  }
  else if (data->runtype == MOLFILE_RUNTYPE_OPTIMIZE) {
    std::cout << "orcaplugin) Reading trajectory of optimization." << std::endl;
    rewind(data->file);
    goto_keyline(data->file, "Geometry Optimization Run", NULL);
  }
  else {
    std::cout << "orcaplugin) Jobtype not supported for trajectory reading." << std::endl;
    return FALSE;
  }

  // printf("orcaplugin) Analyzing trajectory...\n");
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
      printf("orcaplugin) ==== End of trajectory (%d frames) ====\n", data->num_frames);
      data->status = MOLFILE_QMSTATUS_OPT_CONV;
    }
    else if (data->status == MOLFILE_QMSTATUS_OPT_CONV) {
      if(l.find("FINAL ENERGY EVALUATION AT THE STATIONARY POINT") != std::string::npos) {
        if (data->num_frames > 0) {
          data->filepos_array = (long*)realloc(data->filepos_array, (data->num_frames+1)*sizeof(long));
        }
        std::cout << "orcaplugin) found equilibrium geometry." << std::endl;
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
    // std::cout << "have frame" << std::endl;
    int i;
    qm_timestep_t *cur_ts;

    /* get a pointer to the current qm timestep */
    cur_ts = data->qm_timestep+data->num_frames_sent;

    // std::cout << "numwave: " << cur_ts->numwave << std::endl;

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
    // std::cout << "not have frame" << std::endl;
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
  plugin.author = "Maximilian Scheurer, Michael F. Herbst, Marcelo Melo";
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
  /* store the wave function and orbital energies */
  if (cur_ts->wave) {
    std::cout << "orcaplugin) Have wavefunctions: " << cur_ts->numwave << " in frame: " << data->num_frames_sent << std::endl;
    for (i=0; i<cur_ts->numwave; i++) {
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
  //
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
 // move to top later...
static void close_orca_read(void *mydata) {
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
