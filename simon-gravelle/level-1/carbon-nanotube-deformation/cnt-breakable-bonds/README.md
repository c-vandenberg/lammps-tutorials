# Deformation of Carbon Nanotube with Breakable Bonds

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/96af3932-6af6-498f-b4cb-2cddc549ef3a" alt="cnt-breakable-bonds-start md" width="325" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/66c6692a-5ed2-44db-815a-037e76780562" alt="cnt-breakable-bonds-mid md" width="325" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/69d17e3e-b0a0-4462-bead-279f00e56b46" alt="cnt-breakable-bonds-end md" width="325" />
</p>

## Introduction

Reactive molecular dynamics force fields such as AIREBO Force Field extend the capabilities of classical force fields by allowing the simulation to include the breaking & formation of chemical bonds. Therefore, they are used in dynamic systems where chemical reactions are expected to occur (e.g. synthetic & enzymatic reactions).

Another key difference compared to classical force fields is higher complexity of the pre-calculated parameters is higher for reactive force fields. This is because they need to accurately model the potential energy surfaces of reacting systems.

## Differences in Topology File

With OPLS-AA, because it is a classical/non-reactive force field, it requires an explicit list of bonds, angles, dihedral/torsion angles and improper torsions to calculate the system's initial potential energy.

With reactive force fields like AIREBO however, atom connectivity, angles, dihedral/torsion angles and improper torsions (along with other properties required to calculate the system's potential energy) are calculated **dynamically during the simulation**.

Therefore, all bond, angle, dihedral and impropers information must be removed from the topology file.

## Input Script Command Syntax

Most of the commands in `breakable-bonds-input.lammps` input script are the same as `unbreakable-bonds-input.lammps`, however there are some new commands we can break down:

```
# 1) Initialization
## Define simulation temperature, T=300K
variable T equal 300

units metal
atom_style atomic
boundary f f f
pair_style airebo 2.5 1 1
```
* `units metal` - The force field parameters of the AIREBO force field use `metal` units. In the `metal` units system:
  * Distances = **Å**
  * Energy = **eV**
  * Force = **eV/Å**
  * Velocity = **Å/ps**
  * Mass = **atomic mass units (amu)**
  * Time = **picoseconds (ps)**
  * Temperature = **K**
  * Pressure = **bar**

```
# 2) System definition
## Import VMD generated initial carbon nanotube topology. Topology from cnt_molecular.data has been altered to be compatible with the AIREBO potential
read_data cnt_atomic.data

# 3) Simulation settings
## Import parameters/coefficients for the AIREBO potential. This is the file from the original 'Stuart' paper
pair_coeff * * CH.airebo C
```
* `pair_coeff * * CH.airebo C`:
  * Command: `pair_coeff` - Specifies the parameters for the pair potential used to calculate interactions between atoms. In the context of AIREBO, it defines how atoms interact through bonded & non-bonded potentials dynamically during the simulation
  * Parameters:
    * `* *` - Represent the interaction between two given atoms. Here they are wildcards that indicate the specified potential energy parameters apply to all combinations of atoms defined in the system
    * `CH.airebo` - The file containing the potential functions & parameters necessary for the AIREBO force field to calculate all properties of the system. Written by the authors of AIREBO force field (*Stuart et al*)
    * `C` - Specifies the atom type(s) in the `cnt_atomic.data` file. Since there are only carbons present in the carbon nanotube, the single `C` designation suffices
