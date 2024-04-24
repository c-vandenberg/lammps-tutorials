## Level 1 - Lennard-Jones Fluid `improved-input.min.lammps` & `improved-input.min.lammps` Scripts

This tutorial separates out the energy minimization and molecular dynamics stages, making the process more modular and allowing us to restart the simulation from a previously saved configuration (in this case, the system state after energy minimization). 

Additionally, it imposes specific initial positions to the atoms, instead generating their positions randomly. The initial position of the type 2 atoms is within a cylindrical region within the simulation box, and the initial position of the type 1 atoms is outside of this cylindrical region.
