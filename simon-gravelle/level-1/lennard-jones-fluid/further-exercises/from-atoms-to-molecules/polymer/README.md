# Polymer Molecule - Molecule Template & Input Script Command Syntax

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/b616a218-bef6-4604-acd9-7babcbad2eb5" alt="polymer-start md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/b795916e-5353-45cd-b830-8be2697f2522" alt="polymer-mid md" width="300" /> 
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/82054788-16f9-4188-ad53-1c019fb87a1b" alt="polymer-end md" width="300" />
</p>

## Molecule Template - `polymer.mol`

Breaking down the `polymer.mol` file:

**1. Header & General Information Section**
```
# Polymer molecule
9 atoms
8 bonds
7 angles
```
* `9 atoms` - Specifies that the molecule consists of nine atoms
* `8 bond` - Specifies that there are eight bonds within this molecule
* `7 angles` - Specifies that there are seven angles within this molecule

**2. Coordinates Section**
```
Coords
1 -10.0 0 0
2 -7.5 0 0
3 -5.0 0 0
4 -2.5 0 0
5 0.0 0 0
6 2.5 0 0
7 5.0 0 0
8 7.5 0 0
9 10.0 0 0

```
* `Coords` - Demarcates the coordinates section where the positions of the atoms are defined
* `1 -10.0 0 0` - Defines the first atom's ID as `1`, followed by its coordinates `x=-10.0`, `y=0`, `z=0`
* Etc

**3. Types Section**
```
Types
1 2
2 2
3 2
4 2
5 2
6 2
7 2
8 2
9 2
```
* `Types` - Demarcates the types section where the atom types are defined
* `1 2` to `9 2` - All nine atoms are given an ID of `1` to `9` and all are of atom type `2`

**4. Bonds Section**
```
Bonds
1 1 1 2
2 1 2 3
3 1 3 4
4 1 4 5
5 1 5 6
6 1 6 7
7 1 7 8
8 1 8 9
```
* `Bonds` - Demarcates the bonds section where the connections between the atoms are detailed
* `1 1 1 2` to `8 1 8 9` - All 8 bonds are given an ID of `1` to `8`, a bond type of `1` and are connected in a linear fashion

**5. Angles Section**
```
Angles
1 1 1 2 3
2 1 2 3 4
3 1 3 4 5
4 1 4 5 6
5 1 5 6 7
6 1 6 7 8
7 1 7 8 9
```
* `Angles` - Demarcates the angles section where the angles between sets of three connected atoms are defined
* `1 1 1 2 3` - The first `1` is the angle ID, the second `1` is the angle type, and the `1 2 3` specifies the IDs of the three atoms that form the angle

## Energy Minimization Input Script - `polymer.min.lammps`

Breaking down the new commands we encounter in this minimization script:

```
# 3) Simulation settings
mass 1 1
mass 2 1
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
bond_coeff 1 5 2.5
angle_coeff 1 2 180
```
* `angle_coeff 1 2 180` - The `1` refers to the **angle ID** which we defined in the molecule template. The `2` represents the **angle force constant**, and `180` represents the **equilibrium angle in degrees**
