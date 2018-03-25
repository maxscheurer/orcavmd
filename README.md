# orcavmd

Molfile plugin for VMD to read Orca output files.

For compilation, the VMD plugin include directory is needed! (i.e. ../../include)

## Usage
For orbital visualization, add the following keywords to your Orca input file:
```
! [...other input on first line...] Printbasis
%scf print[p_mos] 1 end
```

In VMD, load an Orca file with `mol load orca your_output_file.out `.

## Features
- Load Orca output files (version 3.0.3 and above)
- Jobtypes:
  * Single Point Calculations (ground state, tested for DFT, HF, MNDO/PM3/AM1)
  * Geometry Optimizations
  * Energy+Gradient Calculations (output files can be concatenated in order to visualize QM/MM results)

## Planned Features
- Jobtypes:
   * Excited State Calculations
   * Born-Oppenheimer MD (will be hard, because Orca does not provide sufficient information per timestep)
  

# mopacvmd

Molfile plugin for VMD to read Mopac output files.

__very experimental at the moment! use with caution!__

## Usage
For orbital visualization, add the following keywords in the Mopac input file:
```
EIGEN VECTORS ALLVEC
```

In VMD, load a Mopac file with `mol load mopac your_output_file.out `.

## Features
- Load Mopac output files
- Jobtypes:
  * Single Point Calculations
  * Energy+Gradient Calculations (output files can be concatenated in order to visualize QM/MM results)
 
- Supported Elements:
  * H, O, C, N, S, Cl, F, Br, F
 
