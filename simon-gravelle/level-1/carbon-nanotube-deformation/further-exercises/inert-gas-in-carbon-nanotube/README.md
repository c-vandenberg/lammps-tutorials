# Level 1 - Carbon Nanotube Deformation Further Exercises: Inert Gas (Argon) in The Carbon Nanotube

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/47e5a5fb-dfa4-4951-a655-13f7fd775a96" alt="cnt-breakable-bonds-start md" width="325" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/bf52d51f-4afe-415e-b240-b40b7783ae46" alt="cnt-breakable-bonds-mid md" width="325" />
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/39e0d3db-ef43-461f-aa89-b0b7fc673978" alt="cnt-breakable-bonds-end md" width="325" />
</p>

## Problem

In this exercise we need to modify the input from the unbreakable CNT to add atoms of argon within the CNT. The argon atoms are to be given the following `pair_coeff`:
```
pair_coeff 2 2 0.232 3.3952
```

## Solution

We must modify the `cnt_molecular.data` molecular topography file to specify `2 atom types` and `2 39.948 # Ar`. We also need to add the `pair_coeff` specified in earlier to the `opls-aa.lammps` file:

```
## Bonded & Non-Bonded Potentials for C and Ar. Parameters Taken from OPLS-AA (Optimised Potentials for Liquid Simulations-All-Atom) force field
pair_coeff 1 1 0.066 3.4
pair_coeff 2 2 0.232 3.3952
bond_coeff 1 469 1.4
angle_coeff 1 63 120
dihedral_coeff 1 0 7.25 0 0
improper_coeff 1 5 180
```

We then create a region called `inside_CNT` and populate it randomly with argon atoms via:
* `create_atoms 2 random 40 323485 inside_CNT overlap 1.8 maxtry 50`:
  * `overlap 1.8` - Ensures that new atoms are not placed too close to existing atoms to avoid overlaps. Here, a distance of 1.8 Å is the minimum allowable distance between newly created atoms and existing atoms
  * `maxtry 50` - Specifies the maximum number of attempts to find a non-overlapping position for each new atom. Here, if no suitable position is found after 50 attempts, a new atom is created elsewhere and the process starts again

We then apply the Berendsen thermostat to the CNT C atoms and the Ar atoms separately to avoid having a large temperature difference between the two atom types. This is because:
* C atoms and Ar atoms have different masses and heat capacities and will also likely have different initial temperatures in the simulation
* As a result, if we apply a global thermostat to both simultaneously, this could lead to uneven energy distribution and large temperature differences between the C atoms and Ar atoms

Finally, we anchor the CNT in place/keep the CNT atoms near their starting positions using command `fix myspr carbon_atoms spring/self 5`:
* `spring/self` - Specifies that the fix will apply **harmonic springs** to each atom in the group. These springs will keep each atom near its initial position using a **restoring force**
* `5` - The **spring constant, k** which determines the strength of the harmonic spring
