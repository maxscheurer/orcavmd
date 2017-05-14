# orcavmd

Molfile plugin for VMD to read Orca output files.

For compilation, the VMD plugin include directory is needed! (i.e. ../../include)

# Features
- Load Orca output files (version 3.0.3 and above)
- Jobtypes:
  * Single Point Calculations (ground state, tested for DFT, HF. semiempirical methods)
  * Geometry Optimizations
  * Energy+Gradient Calculations (output files can be appended in order to visualize QM/MM results)

# Planned Features
- Jobtypes:
  * Excited State Calculations
  * Born-Oppenheimer MD (will be hard, because Orca does not provide sufficient information per timestep)
