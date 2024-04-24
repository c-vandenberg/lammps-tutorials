## Level 1 - Lennard-Jones Fluid `improved-input.min.lammps` & `improved-input.min.lammps` Scripts

This tutorial separates out the energy minimization and molecular dynamics stages, making the process more modular and allowing us to restart the simulation from a previously saved configuration (in this case, the system state after energy minimization). 

Additionally, it imposes specific initial positions to the atoms, instead generating their positions randomly like we did with the `first-input.lammps` input script. This allows us to visualise the evolution of the atom populations during mixing of the binary system.

Finally, it extracts the number of atoms in each region as a function of time, and computes the average per atom coordination number per atom between the two atom types.
