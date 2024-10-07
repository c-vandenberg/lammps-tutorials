# 2.5 Further Exercises: From Atoms to Molecules

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/272a447c-9639-4343-a5ce-61f009c66142", alt="dumbbell-molecule-min" width="300"/>
</p>

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/f5e9d1f1-5f44-4eab-bffd-2f3dba1278d2" alt="polymer-start md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/efe0ae66-5742-432c-b00d-ca2cc52260c7" alt="polymer-mid md" width="300" /> 
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/29a578fb-97e3-4676-b100-3c9180fd0cfa" alt="polymer-end md" width="300" />
</p>

## 2.5.1 Problem
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

## 2.5.2 Solution

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
