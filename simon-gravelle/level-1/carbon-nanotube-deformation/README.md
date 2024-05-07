## Level 1 - Lennard-Jones Fluid: Pulling on a Carbon Nanotube

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/6179516c-da6a-4a5c-9b89-7c63428e76f4" alt="cnt-unbreakable-bonds start md" width="500" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/20c9c899-290b-4ea0-9247-dd5721701d88" alt="cnt-unbreakable-bonds-mid md" width="500" />
</p>


The objective of this tutorial is to deform a carbon nanotube (CNT) using LAMMPS. A small CNT will be simulated within an empty simulation box, an external force will be exerted on the CNT, and its deformation will be measured over time.

This tutorial will illustrate the difference between classical and reactive force fields:
* With a classical force field (in this case OPLS-AA Force Field), the bond between the atoms are unbreakable
* With a reactive force field (in this case AIREBO Force Field), the breaking of chemical bonds is possible when the deformation is strong enough
