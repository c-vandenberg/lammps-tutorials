## Level 1 - Carbon Nanotube Deformation

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/ea6277cb-2d37-4ccf-8c48-1b6ba812d904" alt="cnt-membrane" width="" />
</p>

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/54201568-22cb-4021-96c0-e3fba125e62b" alt="cnt-unbreakable-bonds end md" width="500" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/1a38a1ca-9d88-47aa-82c9-3bb6225924c5" alt="cnt-breakable-bonds-end md" width="500" />
</p>

The objective of this tutorial is to deform a carbon nanotube (CNT) using LAMMPS. A small CNT will be simulated within an empty simulation box, an external force will be exerted on the CNT, and its deformation will be measured over time.

This tutorial will illustrate the difference between classical and reactive force fields:
* With a classical force field (in this case OPLS-AA Force Field), the bond between the atoms are unbreakable
* With a reactive force field (in this case AIREBO Force Field), the breaking of chemical bonds is possible when the deformation is strong enough
