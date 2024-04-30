## Level 1 - Lennard-Jones Fluid `improved-input.min.lammps` & `improved-input.min.lammps` Scripts

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/db4cc9a4-0c3a-4368-ae02-580f5bf4c890" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/f98f2826-41f1-4b97-b23e-ae92d54b9e05" width="300" /> 
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/64c8f927-f7c7-4336-919a-2098f40b0d8d" width="300" />
</p>

This tutorial separates out the energy minimization and molecular dynamics stages, making the process more modular and allowing us to restart the simulation from a previously saved configuration (in this case, the system state after energy minimization). 

Additionally, it imposes specific initial positions to the atoms, instead generating their positions randomly like we did with the `first-input.lammps` input script. This allows us to visualise the evolution of the atom populations during mixing of the binary system.

Finally, it extracts the number of atoms in each region as a function of time, and computes the average per atom coordination number per atom between the two atom types.
