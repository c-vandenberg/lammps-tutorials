# Deforming/Stretching the Water Solvated PEG Polymer Molecule

https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/aa14bf6b-29a2-475f-8f1a-cdfd3ae44f42

## Introduction

We will now apply a constant force to the oxygen atoms on both ends of the PEG polymer molecule until it stretches. This will allow us to probe the deformation via the changing intra-molecular parameters (see `PEG-deformation-evaluation` further exercise). 

The force magnitude was chosen via trial & error to be large enough to overcome the thermal agitation & the entropic contribution from both water and PEG molecules.

## Data Analysis
<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/11110a72-6d6a-448a-ac33-b6df8c9b7294" alt ="peg length vs time" width="" />
</div>

The PEG end-to-end length is defined as the Euclidean distance between the center of mass coordinates for the two oxygen atoms at each end of the PEG molecule. Euclidean distance between any two points in a 3D Euclidean space is given by the equation:

<br>
<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?\color{white}d%20=%20\sqrt{(x_2%20-%20x_1)^2%20+%20(y_2%20-%20y_1)^2%20+%20(z_2%20-%20z_1)^2}" alt ="euclidean distance formula" width="" />
</div>
<br>

As you can see, the force is applied to the oxygen atoms at each end of the PEG molecule at ~ *t* = 31.5 ps. Maximum PEG length is reached at ~ *t* = 54 ps, indicating that the molecule is fully stretched.


## Input Script Command Syntax

Breaking down the new commands we encounter in this input script:

```
# 4) Run
variable oxygen_6_x_coordinate equal xcm(pull_oxygen_6,x)
variable oxygen_6_y_coordinate equal xcm(pull_oxygen_6,y)
variable oxygen_6_z_coordinate equal xcm(pull_oxygen_6,z)
variable oxygen_7_x_coordinate equal xcm(pull_oxygen_7,x)
variable oxygen_7_y_coordinate equal xcm(pull_oxygen_7,y)
variable oxygen_7_z_coordinate equal xcm(pull_oxygen_7,z)
variable oxygen_6_7_distance equal sqrt((v_oxygen_6_x_coordinate-v_oxygen_7_x_coordinate)^2&
+(v_oxygen_6_y_coordinate-v_oxygen_7_y_coordinate)^2+(v_oxygen_7_x_coordinate-v_oxygen_7_z_coordinate)^2)
```
* Here, the end-to-end (Euclidean) distance is calculated using the center of mass coordinates for the x, y & z coordinates for each group
* Euclidean distance between any two points in a 3D Euclidean space is given by the equation `sqrt((x2 − x1)^2 + (y2 − y1)^2 + (z2 − z1)^2)`

