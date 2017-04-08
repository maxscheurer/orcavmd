/*
Orca VMD plugin
Authors: Maximilian Scheurer, Marcelo Melo, April 2017
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


static void* open_orca_read(const char* filename, const char* filetype, int *natoms) {
  FILE* fd;
  qmdata_t *data = NULL;
  printf("Open Orca Read called; %s\n", filename);
  *natoms = 0;
}

static int read_orca_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
{
  qmdata_t *data = (qmdata_t *)mydata;
  qm_atom_t *cur_atom;
  molfile_atom_t *atom;
  int i = 0;
  return 0;
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
