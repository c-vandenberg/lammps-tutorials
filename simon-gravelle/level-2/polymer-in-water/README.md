## Level 2 - Polymer in Water

The objective of this tutorial is to use LAMMPS to solvate the small hydrophobic polymer PolyEthylene Glycol (PEG) in a reservoir of water.

An all-atom description is used for both PEG and water, and the long range Coulomb interactions are solved using the [PPPM method](https://docs.lammps.org/Howto_dispersion.html). Once the PEF and water system are properly equilibrated at the desired temperature & pressure, a constant stretching force will be applied to both ends of the PEG, and the changes in its length will be measured as a function of time.

This tutorial was inspired by a publication by Liese *et al* <sup>1</sup>, in which molecular dynamics simulations are compared with single-molecule force spectroscopy.

### References
[1] Liese, S. *et al*. (2016) ‘Hydration effects turn a highly stretched polymer from an entropic into an energetic spring’, *ACS Nano*, 11(1), pp. 702–712.