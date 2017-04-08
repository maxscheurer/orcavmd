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

/*
 * Error reporting macro for the multiple fgets calls in
 * the code
 */
#define GET_LINE(x,y) if (!fgets(x, sizeof(x), y)) return FALSE

#define NOTFOUND 0
#define FOUND    1
#define STOPPED  2

/* ######################################################## */
/*                 orca specific struct                     */
/*
  might be extended in the future...
 */
/* ######################################################## */

typedef struct {
  int version; /* version = 0 : unknown, do not read the file!
                * version = 1 : version 3.0.0, 3.0.1, 3.0.3
                * version = 2 : version 4.0.0
                */
} orcadata;


/* ######################################################## */
/*                    static functions                      */
/* ######################################################## */

// Checks if loaded file is really an Orca output file
static int have_orca(qmdata_t *data, orcadata* orca);




/*************************************************************
 *
 * MAIN ORCA CODE PART
 *
 **************************************************************/

static void* open_orca_read(const char* filename, const char* filetype, int *natoms) {
  FILE* fd;
  qmdata_t *data = NULL;

  #ifdef DEBUGGING
    printf("Open Orca Read called: %s\n", filename);
  #endif

  orcadata* orca;

  // THIS NEEDS TO BE REMOVED LATER!!!
  *natoms = 0;
  // !!!!

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

  // filename in qm data struct
  data->file = fd;

  if(have_orca(data, orca)) {

  } else {
    printf("orcaplugin) This is not an Orca output file!\n");
    return NULL;
  }

  return data;
}

// static int read_orca_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
// {
//   qmdata_t *data = (qmdata_t *)mydata;
//   qm_atom_t *cur_atom;
//   molfile_atom_t *atom;
//   int i = 0;
//   return 0;
// }



static int have_orca(qmdata_t *data, orcadata* orca) {
  int programLine;
  programLine = goto_keyline(data->file, "O   R   C   A", NULL);
  if (programLine != 1) {
    return FALSE;
  }
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
  // plugin.read_structure = read_orca_structure;
  // plugin.close_file_read = close_orca_read;
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
