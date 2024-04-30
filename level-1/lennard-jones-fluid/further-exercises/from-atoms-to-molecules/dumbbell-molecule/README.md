# Dumbbell Molecule - Molecule Template & Input Script Command Syntax

## Molecule Template - `dumbbell.mol`
A `.mol` file in LAMMPS defines the structure of a molecule, including its atoms, their positions/coordinates, types and the bonds between them. 

Breaking down the `dumbbell.mol` file:

**1. Header & General Information Section**
```
# Dumbell molecule
2 atoms
1 bonds
```
* `2 atoms` - Specifies that the molecule consists of two atoms
* `1 bond` - Specifies that there is one bond within this molecule

**2. Coordinates Section**
```
Coords
1 1 0 0
2 -1 0 0
```
* `Coords` - Demarcates the coordinates section where the positions of the atoms are defined
* `1 1 0 0` - Defines the first atom's ID as `1`, followed by its coordinates `x=1`, `y=0`, `z=0`. The atom is placed one unit distance along the x-axis
* `2 -1 0 0` - Defines the second atom's ID as `2`, followed by its coordinates `x=-1`, `y=0`, `z=0`. This atom is placed on the opposite side of the origin along the x-axis, making this molecule symmetric about the origin

**3. Types Section**
```
Types
1 2
2 2
```
* `Types` - Demarcates the types section where the atom types are defined
* `1 2` - Atom `1` is given type `2`
* `2 2` - Atom `2` is given type `2`

**4. Bonds Section**
```
Bonds
1 1 1 2
```
* `Types` - Demarcates the bonds section where the connections between the atoms are detailed
* `1 1 1 2` - The first `1` is the bond ID, the second `1` is the bond type and `1 2` are the IDs of the atoms connected by this bond. This specifies a single bond between atoms `1` and `2`

## Energy Minimization Input Script - `dumbbell-molecule.min.lammps`

Breaking down the new commands we encounter in this minimization script:

```
# 1) Initialization
units lj
dimension 3
atom_style molecular
bond_style harmonic
pair_style lj/cut 2.5
boundary p p p
```
* `atom_style molecular` - Specifies that the atoms can be connected by bonds
* `bond_style harmonic` - Specifies that the bonds are to behave according to a harmonic potential

```
# 2) System definition
# Cubic Box
region simulation_box block -20 20 -20 20 -20 20

# Small Elongated Box
# region simulation_box block -20 20 -6 6 -6 6

create_box 2 simulation_box bond/types 1 extra/bond/per/atom 1

# Import and use molecule template to populate box with molecules
molecule dumbbell ../dumbbell.mol
create_atoms 1 random 1500 341341 simulation_box
create_atoms 0 random 20 678865 simulation_box mol dumbbell 8751
```
* `create_box 2 simulation_box bond/types 1 extra/bond/per/atom 1`:
  * The `bond/types 1` parameter specifies that there will be one type of bond in the system. 
  * The `extra/bond/per/atom 1` parameter specifies that there should be memory allocation for one additional bond per atom. 
  * This is useful in systems where bonds will need to be dynamically created or destroyed during the simulation, or when the maximum number of bonds any given atom may have at the start of the simulation isn't known. 
  * It ensures LAMMPS allocates sufficient memory to handle these bonds without needing to reallocate memory during runtime

```
# 3) Simulation settings
mass 1 1
mass 2 1
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
bond_coeff 1 5 2.0
```
* `bond_coeff 1 5 2.0` - The `1` refers to the **bond ID** which we defined in the molecule template. The `5` is the **energy constant**, which for harmonic potential, would be the **spring constant, k** (affects bond stiffness). The `2.0` is the **equilibrium distance**, the ideal length of the bond at which there is no force exerted by the bond on the bonded atoms
