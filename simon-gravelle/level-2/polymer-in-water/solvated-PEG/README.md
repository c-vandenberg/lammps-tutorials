# Solvating the PEG in Water

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