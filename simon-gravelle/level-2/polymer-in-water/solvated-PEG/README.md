# Solvating the PEG in Water

https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/58a02fc5-91b0-4af8-b660-782eb3208f98

## Introduction

Now that we have created & equilibrated both our water reservoir and single PEG polymer, we can now combine the two systems to solvate the PEG in the water reservoir.

## Input Script Command Syntax

Breaking down the new commands we encounter in this input script:

```
# 1) Initialization
units real
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style harmonic
pair_style lj/cut/coul/long 12
kspace_style pppm 1e-5
special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes
```
* `special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes dihedral yes`:
  * `angle yes` - Specifies that the 1-4 interactions weighting factor will be ignored if the atoms are not listed as the first and last atoms in any dihedral defined in the simulation

```
# 2) System Definition
## Import the two previously generated data file H2O.data and PEG.data. Again we are using the `extra/x/per/atom` commands for memory allocation
## When using the `read_data` command more than once. we need to use the `add append` arguments
read_data ../pure-H2O/H2O.data &
    extra/bond/per/atom 3 &
    extra/angle/per/atom 6 &
    extra/dihedral/per/atom 10 &
    extra/special/per/atom 14
read_data ../single-PEG/PEG.data add append
```
* `extra/x/per/atom` - Allows us to build a composite system combining both the water reservoir and PEG molecule, with sufficient memory allocation to handle any additional interactions that may occur during the simulation. This allows for a flexible & dynamic system where new bonds & angles can form if necessary
* `read_data ../single-PEG/PEG.data add append` - When using the `read_data` command more than once. we need to use the `add append` arguments to append any subsequent data to the existing system
