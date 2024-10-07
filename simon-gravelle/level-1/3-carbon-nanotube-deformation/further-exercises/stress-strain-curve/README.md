# 3.3 Further Exercises: Plot the Carbon Nanotube Stress-Strain Curves

## 3.3.1 Problem

In this exercise we will calculate stress-strain curve for both breakable bonds CNT and unbreakable bonds CNT.

## 3.3.2 Solution

* The stress is calculated as the total force applied to the CNT from the deformation, divided by the surface area of the CNT
* We also need to include constants and conversion factors in our scripts, taking note of the difference in units between the OPLS-AA & AIREBO Force Fields

## 3.3.3 Data Analysis

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/9a6cd7c6-8432-4dbf-9265-c1d2ae5b10b3" alt="stress_vs_strain" width="" />
</p>

The stress vs strain trend for the unbreakable CNT molecule is similar to that seen for the [CNT length vs time graphs](https://github.com/c-vandenberg/lammps-tutorials/blob/master/simon-gravelle/level-1/carbon-nanotube-deformation/cnt-unbreakable-bonds/README.md#carbon-nanotube-length-during-deformation). We see a linear relationship between stress and strain that is continuous due to the chemical bonds being unbreakable.

The breakable CNT molecule on the other hand shows a non-linear relationship between stress and strain that suddenly drops once stress reaches ~ 500 GPa. This indicates that the stress threshold for CNT chemical bond stretching is ~ 500 GPa, after which the bonds will break.
