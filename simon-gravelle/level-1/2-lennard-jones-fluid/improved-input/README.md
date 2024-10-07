## 2.2 Level 1 - Lennard-Jones Fluid `improved-input.min.lammps` & `improved-input.min.lammps` Scripts

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/3861bbcf-e3af-40d8-a078-55b5b89d27f1" alt="improved-input-start md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/84e2ca5b-d80d-476f-a99c-170372d4ae3c" alt="improved-input-mid md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/10ca7534-edd6-4a93-82cf-d2009784eaba" alt="improved-input-end md" width="300" />
</p>

This tutorial separates out the energy minimization and molecular dynamics stages, making the process more modular and allowing us to restart the simulation from a previously saved configuration (in this case, the system state after energy minimization). 

Additionally, it imposes specific initial positions to the atoms, instead generating their positions randomly like we did with the `first-input.lammps` input script. This allows us to visualise the evolution of the atom populations during mixing of the binary system.

Finally, it extracts the number of atoms in each region as a function of time, and computes the average per atom coordination number per atom between the two atom types.
