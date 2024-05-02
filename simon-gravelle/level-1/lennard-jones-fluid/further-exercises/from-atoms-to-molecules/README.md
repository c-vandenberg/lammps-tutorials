# Level 1 - Lennard-Jones Fluid Further Exercises: From Atoms to Molecules

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/66a60bcc-6d64-4be1-9df2-9975ffd02821" width="300"/>
</p>

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/b616a218-bef6-4604-acd9-7babcbad2eb5" alt="polymer-start md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/b795916e-5353-45cd-b830-8be2697f2522" alt="polymer-mid md" width="300" /> 
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/82054788-16f9-4188-ad53-1c019fb87a1b" alt="polymer-end md" width="300" />
</p>

## Problem
Add a bond between type 2 atoms in order to create dumbbell molecules instead of single atoms. Additionally, create a small polymer/long chain of atoms linked by bonds and defined by angles

**Hint**


Use a **molecule template** to easily insert as many atoms connected by bonds (i.e. molecules) as you want. A molecule template typically begins as follows:
```
2 atoms
1 bonds

Coords

1 0.5 0 0
2 -0.5 0 0

(...)
```
A bond section also needs to be added.

## Solution

### Dumbbell Molecule

The first change to the input script is to choose an `atom_style` that allows for the atoms to be connected by bonds (i.e. `molecular`). It is also necessary tp specify the `bond_style` (i.e. the **type of potential that will keep the bonds together**). Here it is `harmonic`:
```
atom_style molecular
bond_style harmonic
```

Additionally, when creating the simulation region we need to allocate space in memory for the bond:
```
create_box 2 simulation_box bond/types 1 extra/bond/per/atom 1
```

We also need to create a **molecule template** and then import & use it when creating the atoms:
```
molecule dumbell dumbell.mol
create_atoms 1 random 500 341341 simulation_box
create_atoms 0 random 5 678865 simulation_box mol dumbell 8754
```

Lastly, we need to specify some parameters for the bond, namely its **rigidity (5)** and **equilibrium length (2.5)**
```
bond_coeff 1 5 2.5
```

### Polymer Molecule
