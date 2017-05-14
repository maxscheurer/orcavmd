/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: main.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.16 $       $Date: 2016/11/28 05:01:54 $
 *
 ***************************************************************************/

/*
 * A general main for testing plugins.
 * Compile using: gcc main.c plugin.c -I../../include -o plugintest
 * Replace plugin.c with the plugin file you want to test.
 * Usage: plugintest <filetype> <file> [<filetype> <file>]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "molfile_plugin.h"

/* Structure and coordintes plugin */
static molfile_plugin_t *splugin = 0;
static molfile_plugin_t *cplugin = 0;
static const char *sfiletype = NULL;
static const char *cfiletype = NULL;

static int register_cb(void *v, vmdplugin_t *p) {
  if (!strcmp(p->type, MOLFILE_PLUGIN_TYPE)) {
    if (!strcmp(p->name, sfiletype))
      splugin = (molfile_plugin_t *)p;
    if (!strcmp(p->name, cfiletype))
      cplugin = (molfile_plugin_t *)p;
  }
  return VMDPLUGIN_SUCCESS;
}

int main(int argc, char *argv[]) {
  const char *sfilename;
  const char *cfilename;
  int rc, natoms;
  void *handle;
  void *chandle;
  molfile_timestep_t timestep;

  if (argc != 3 && argc != 5) {
    fprintf(stderr, "Usage: %s <filetype> <filename> [<filetype> <filename>]\n", argv[0]);
    return 1;
  }

  if (argc == 3) {
    sfiletype = argv[1];
    sfilename = argv[2];
    cfiletype = sfiletype;
    cfilename = sfilename;
  } else {
    sfiletype = argv[1];
    sfilename = argv[2];
    cfiletype = argv[3];
    cfilename = argv[4];
  }

  vmdplugin_init();
  vmdplugin_register(NULL, register_cb);

  if (!splugin) {
    fprintf(stderr, "No plugin for filetype %s was linked in!\n", sfiletype);
    return 1;
  }
  if (!cplugin) {
    fprintf(stderr, "No plugin for filetype %s was linked in!\n", cfiletype);
    return 1;
  }

  /* Read structure */
  if (!splugin->open_file_read) {
    fprintf(stdout, "FAILED: No open_file_read found in structure plugin.\n");
    return 1;
  }
  handle = splugin->open_file_read(sfilename, sfiletype, &natoms);
  if (!handle) {
    fprintf(stderr, "FAILED: open_file_read returned NULL in structure plugin.\n");
    return 1;
  }
  printf("Opened file %s; structure plugin found %d atoms\n", sfilename, natoms);
  if (splugin->read_structure) {
    int optflags;
    molfile_atom_t *atoms;
    atoms = (molfile_atom_t *)malloc(natoms * sizeof(molfile_atom_t));
    rc = splugin->read_structure(handle, &optflags, atoms);
    free(atoms);
    if (rc) {
      fprintf(stderr, "FAILED: read_structure returned %d\n", rc);
      splugin->close_file_read(handle);
      return 1;
    } else {
      printf("Successfully read atom structure information.\n");
    }
    if (splugin->read_bonds) {
      int nbonds, *from, *to, *bondtype, nbondtypes;
      float *bondorder;
      char **bondtypename;
      if ((rc = splugin->read_bonds(handle, &nbonds, &from, &to,
				   &bondorder, &bondtype, &nbondtypes, &bondtypename))) {
        fprintf(stderr, "FAILED: read_bonds returned %d\n", rc);
      } else {
        printf("read_bonds read %d bonds\n", nbonds);
      }
    } else {
      printf("Structure file contains no bond information\n");
    }
  } else {
    fprintf(stderr, "FAILED: File contains no structure information!\n");
    return 1;
  }

  /* Check whether we use one plugin for both structure and coords */
  if (splugin != cplugin) {
    splugin->close_file_read(handle);
    int cnatoms;
    chandle = cplugin->open_file_read(cfilename, cfiletype, &cnatoms);
    printf("Opened coordinates file %s\n", cfilename);
    if (cnatoms != MOLFILE_NUMATOMS_UNKNOWN && cnatoms != natoms) {
      fprintf(stderr, "FAILED: Different number of atoms in structure file (%d) than in coordinates file (%d)!",
	      natoms, cnatoms);
      cplugin->close_file_read(chandle);
      exit(1);
    }
  } else {
    chandle = handle;
  }

  /* Read coordinates */
  if (cplugin->read_timestep) {
    timestep.velocities = NULL;
    int nsteps = 0;
    timestep.coords = (float *)malloc(3*natoms*sizeof(float));
    while (1) {
      molfile_qm_metadata_t *qm_metadata = NULL; // this just a dummy
      molfile_qm_timestep_t qm_timestep;
      molfile_qm_timestep_metadata_t qmmeta;
      memset(&qmmeta, 0, sizeof(molfile_qm_timestep_metadata_t));
      // XXX need to add timestep parameter or other method to specify
      //     which frame this applies to, else keep it the way it is
      //     and rename the plugin function appropriately.
      cplugin->read_qm_timestep_metadata(chandle, &qmmeta);
      memset(&qm_timestep, 0, sizeof(molfile_qm_timestep_t));
      qm_timestep.wave = new molfile_qm_wavefunction_t[qmmeta.num_wavef];
      // printf("%d\n", qmmeta.num_wavef);
      memset(qm_timestep.wave, 0, qmmeta.num_wavef*sizeof(molfile_qm_wavefunction_t));
      // printf("%p\n", qm_timestep.wave[0].wave_coeffs);
      int i;
      for (i=0; (i<MOLFILE_MAXWAVEPERTS && i<qmmeta.num_wavef); i++) {
        qm_timestep.wave[i].wave_coeffs =
          new float[qmmeta.num_orbitals_per_wavef[i]*qmmeta.wavef_size];
        if (qmmeta.has_orben_per_wavef[i]) {
          qm_timestep.wave[i].orbital_energies =
            new float[qmmeta.num_orbitals_per_wavef[i]];
        }
        if (qmmeta.has_occup_per_wavef[i]) {
          qm_timestep.wave[i].occupancies =
            new float[qmmeta.num_orbitals_per_wavef[i]];
        }
      }
      rc = cplugin->read_timestep(chandle, natoms, &timestep, qm_metadata, &qm_timestep);
      if (rc != MOLFILE_SUCCESS) {
        break;
      }
      nsteps++;

      delete [] qm_timestep.scfenergies;
      if (qm_timestep.gradient) delete [] qm_timestep.gradient;
      if (qm_timestep.charges)  delete [] qm_timestep.charges;
      if (qm_timestep.charge_types) delete [] qm_timestep.charge_types;

      for (i=0; i<qmmeta.num_wavef; i++) {
        delete [] qm_timestep.wave[i].wave_coeffs;
        delete [] qm_timestep.wave[i].orbital_energies;
        delete [] qm_timestep.wave[i].occupancies;
      }
      delete [] qm_timestep.wave;
    }

    free(timestep.coords);

    if (rc != MOLFILE_SUCCESS) {
      fprintf(stderr, "FAILED: read_next_timestep returned %d\n", rc);
    } else {
      printf("Successfully read %d timesteps\n", nsteps);
    }
  }

  /* Close plugin(s) */
  cplugin->close_file_read(chandle);

  vmdplugin_fini();
  printf("Tests finished.\n");
  return 0;
}
