# Deformation of Carbon Nanotube with Unbreakable Bonds

## Introduction

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/2a77460a-021c-4732-a78c-a08935f35ef9" alt="cnt-unbreakable-bonds start md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/707aa6f3-5138-4259-88fc-754b4c901c32" alt="cnt-unbreakable-bonds-mid md" width="300" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/bca5260b-a3b9-4d1d-8e1b-d3f7fe2617d7" alt="cnt-unbreakable-bonds-end md" width="300" />
</p>

In classical molecular dynamics force fields (also known as non-reactive force fields) such as OPLS-AA Force Field, the chemical bonds between the atoms are set at the start of the simulation. As the simulation is running, regardless of the forces applied to the atoms, the bonds will remain intact. Therefore, they are designed to model the interaction within molecules and materials where the bonding structure does not change throughout the simulation. 

The bonds between neighbouring atoms are typically modelled as springs with a given equilibrium distance r<sub>0</sub> and a spring constant k<sub>b</sub>: U<sub>b</sub> = k<sub>b</sub>(r - r<sub>0</sub>)<sup>2</sup>.

Additionally, angular and dihedral constraints are usually applied to maintain the relative orientations of neighbour atoms

## Input Script Command Syntax

Breaking down the new commands we encounter in this input script:

### 1) Initialization
```
## Define simulation temperature, T=300K
variable T equal 300

units real
atom_style molecular
boundary f f f
pair_style lj/cut 14

bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style harmonic

special_bonds lj 0.0 0.0 0.5
```
* `variable T equal 300` - Defines the simulation temperature, T=300K
* `units real` - Specifies the simulation is using `real` units instead of the unitless `lj` units we have used previously. In the `real` units system:
  * Distances = **Å**
  * Energy = **kcal/mol**
  * Force = **kcal/mol/Å**
  * Velocity = **Å/fs**
  * Mass = **atomic mass units (amu)**
  * Time = **picoseconds (fs)**
  * Temperature = **K**
  * Pressure = **atm**
* `bond_style harmonic` - Bonds are modelled as harmonic springs in OPLS-AA
* `angle_style harmonic`- Bond angles are modelled as harmonic springs in OPLS-AA
* `dihedral_style opls` - Dihedral/torsion angles can exhibit more complex behaviour than simple harmonic motion because they involve the interaction of four atoms and can have multiple energy minima. The OPLS style is therefore used as it uses a Fourier series, allowing it to cross the minima barriers and capture the multiple minima
* `improper_style harmonic` - Improper torsions are modelled as harmonic springs in OPLS-AA
* `special_bonds lj 0.0 0.0 0.5` - Atoms connected by a bond typically do not interact through Lennard-Jones interation. Therefore, they must be excluded from the Lennard-Jones potential calculation. This command scales:
  * **1-2 Lennard-Jones interactions** (i.e. atoms directly bonded) by **0.0**, ignoring these interactions
  * **1-3 Lennard-Jones interactions** (i.e. atoms two bonds apart) by **0.0**, ignoring these interactions
  * **1-4 Lennard-Jones interactions** (i.e. atoms three bonds apart) by **0.5**, halving the calculated interaction value

```
# 3) Simulation settings
## Import OPLS-AA Force Field Bonded & Non-Bonded Potentials
include opls-aa.lammps

## Calculate the negative of the x, y & z centre of mass (cm) components of the carbon_atoms group
## Recenter the carbon nanotube at the origin (0, 0, 0) by displacing them by the calculated center of mass variables
group carbon_atoms type 1
variable carbon_xcm equal -1*xcm(carbon_atoms,x)
variable carbon_ycm equal -1*xcm(carbon_atoms,y)
variable carbon_zcm equal -1*xcm(carbon_atoms,z)
displace_atoms carbon_atoms move ${carbon_xcm} ${carbon_ycm} ${carbon_zcm}

## Center the simulation box about the origin (0, 0, 0) by altering dimensions to be from -40 to 40 in x, y & z directions
change_box all x final -40 40 y final -40 40 z final -40 40

## Define top edge region, bottom edge region, & middle region of carbon nanotube
variable zmax equal bound(carbon_atoms,zmax)-0.5
variable zmin equal bound(carbon_atoms,zmin)+0.5
region region_top block INF INF INF INF ${zmax} INF
region region_bottom block INF INF INF INF INF ${zmin}
region region_middle block INF INF INF INF ${zmin} ${zmax}

## Group atoms according to these regions
group carbon_top region region_top
group carbon_bottom region region_bottom
group carbon_middle region region_middle
```
* `variable carbon_xcm equal -1*xcm(carbon_atoms,x)` - Defines variable `carbon_xcm` and sets it equal to the negative of the center of mass x-coordinates of the `carbon_atoms` group. The same is done for the y and z coordinates. These calculations are used to determine how much each atom needs to be moved to center the nanotube at the origin
* `displace_atoms carbon_atoms move ${carbon_xcm} ${carbon_ycm} ${carbon_zcm}` - Move all atoms in the `carbon_atoms` group based on these calculated negative center of mass coordinates
* `change_box all x final -40 40 y final -40 40 z final -40 40` - Alters the dimensions of the simulation box for all three axes to be -40 Å to 40 Å. This ensures the box is symmetrical about the origin
* `variable zmax equal bound(carbon_atoms,zmax)-0.5` - Defines variable `zmax` and sets it equal to `bound(carbon_atoms,zmax)-0.5`
  * `bound(carbon_atoms,zmax)-0.5` - The `bound` function has two arguments; `group` and `boundary` (`bound(group, boundary)`), and it returns the coordinate boundary of a specified group of atoms along the x, y or z dimension (either the maximum or minimum boundary). This boundary can be reduced to shrink the zone of interest (e.g. by `-0.5 Å` in this case)
  * Here, this boundary function calculates the coordinate of the maximum boundary in the z-direction, -0.5 Å
  * The next boundary function calculates the coordinate of the minimum boundary in the z-direction, +0.5 Å
* `region region_top block INF INF INF INF ${zmax} INF` - Defines region `region_top` of shape `block`. The parameters `INF INF INF INF ${zmax} INF` define the spatial bounds of the block:
  * `INF INF` - The first two parameters specify the -x and +x bound respectively. Here they both extend infinitely
  * `INF INF` - The second two parameters specify the -y and +y bound respectively. Here they both extend infinitely
  * `${zmax} INF` - The third two parameters specify the -z and +z bound respectively. Here the lower bond is set as the previously calculated `${zmax}`, whereas the upper bound extends infinitely. This defines a region from `${zmax}` upwards infinitely in the z-direction
  * The same is then done for a `region_bottom` region `INF ${zmin}` and a middle region inbetween `${zmin}` and `${zmax}`

```
# 4) Visualization
## Randomly delete non-edge carbon atoms to enhance visualization and highlight stretching and compression of carbon nanotube
variable zmax_delete equal ${zmax}-2
variable zmin_delete equal ${zmin}+2
region region_delete block INF INF INF INF ${zmin_delete} ${zmax_delete}
group carbon_delete region region_delete
delete_atoms random fraction 0.02 no carbon_delete NULL 482793 bond yes
```
* `variable zmax_delete equal ${zmax}-2` - Defines variable `zmax_delete` and sets it equal to `${zmax}` - 2 Å. The same is done with `${zmin}` + 2 Å. This narrows the atom deletion region to preserve the edge atoms in the top and bottom regions (as these are what will be applying the force later)
* `region region_delete block INF INF INF INF ${zmin_delete} ${zmax_delete}` - Defines region `region_delete` of shape `block` and sets it as between` ${zmin_delete`} and `${zmax_delete}` in the z-direction, and extending infinitely in the x and y directions
* `group carbon_delete region region_delete` - Defines group `carbon_delete` that includes all atoms located within `region_delete`
* `delete_atoms random fraction 0.02 no carbon_delete NULL 482793 bond yes`:
  * `random` - Specifies atoms should be deleted randomly
  * `fraction 0.02` - Indicates that 2% of the atoms in the specified group (`carbon_delete`) should be deleted
  * `no` - Ensures that deletion does not use a type-based criterion
  * `carbon_delete` - Group from which atoms will be randomly deleted 
  * `NULL` - Typically used to specify a type argument, but here it indicates no specific type is targeted. 
  * `482793` - Seed number for the random number generator
  * `bond yes` - Removes bonds involving deleted atoms. This maintains the consistency of the molecular topology.

```
# 5) Run
## Re-order atoms by their IDs since atoms have deleted
reset_atoms id sort yes

## Initialize velocities of atoms in the carbon_middle group at T (300 K) with no overall translational momentum (`mom yes``) or rotational momentum (`rot yes`)
## As a result of `mom yes` and `rot yes`, the systems net movement and rotation is null
velocity carbon_middle create ${T} 48455 mom yes rot yes

## Apply the NVE ensemble to the system
fix nve_ensemble all nve

## Calculate the temperature of the carbon_middle group
compute carbon_middle_temp carbon_middle temp

## Apply a Berendsen thermostat to the carbon_middle group to regulate its temperature to ${T} with a relaxation time of 100 time units
fix berendsen_thermostat carbon_middle temp/berendsen ${T} ${T} 100
```
* `fix berendsen_thermostat carbon_middle temp/berendsen ${T} ${T} 100`:
  * Defines fix `berendsen_thermostat` that applies to the carbons in group `carbon_middle`
  * Fix is of type `temp/brendsen`, and specifies both the initial & target temperature as `T=300K`
  * The `100` parameter is the time constant in femtoseconds for the thermostat (i.e. how quickly the system's temperature is adjusted towards the target temperature. Smaller value=higher adjustment)
  * The Berendsen thermostat algorithm is designed to weakly couple the system to a heat bath with a specific temperature by adjusting the velocities of the atoms in the group to scale their kinetic energy towards the target temperature