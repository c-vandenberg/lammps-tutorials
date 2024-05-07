## Level 1 - Lennard-Jones Fluid: Pulling on a Carbon Nanotube

The objective of this tutorial is to deform a carbon nanotube (CNT) using LAMMPS. A small CNT will be simulated within an empty simulation box, an external force will be exerted on the CNT, and its deformation will be measured over time.

This tutorial will illustrate the difference between classical and reactive force fields:
* With a classical force field (in this case OPLS-AA Force Field), the bond between the atoms are unbreakable
* With a reactive force field (in this case AIREBO Force Field), the breaking of chemical bonds is possible when the deformation is strong enough
