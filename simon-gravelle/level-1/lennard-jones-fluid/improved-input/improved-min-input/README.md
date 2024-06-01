# Level 1 - Lennard-Jones Fluid `improved_input.min.lammps` Script

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/3be49a8e-ac39-4757-96f9-6e841e181149" alt="improved-input min" width="500">
</p>

This `improved-input.min.lammps` input script creates a cubic 3-dimensional simulation box, with periodic boundary conditions, and the same atomic configuration and unit type as `first-input.lammps`. Unlike `first-input.lammps` however, it also two additional regions for placing the atoms:
* A cylindrical region - The initial position of the type 2 atoms
* A region outside this cylindrical region - The initial position of the type 1 atoms

Another novelty compared to `first-input.lammps` is that the `improved_input.min.lammps` script outputs the final state of the system following energy minimization to a `minimized_coordinate.data` file which can be later used to restart the simulation from this final state.

## Input Script Command Syntax

Most of the commands are present in the `first-input.lammps` script and have already been explained. There are some new ones however:

### 2) System definition
* `region region_inside_cylinder cylinder z 0 0 10 INF INF side in`: Defines a cylindrical region along the z-axis centered coordinates x=0 & y=0 (i.e. the origin), with a radius of 10 (unitless) and infinite height. The `side in` argument specifies that the region is the space inside the cylinder
* `region region_outside_cylinder cylinder z 0 0 10 INF INF side out`: Defines a region outside the aforementioned cylinder
* `create_atoms 1 random 1000 341341 region_outside_cylinder`: Randomly places 1000 atoms of type 1 within the region outside the cylinder, using integer 341341 as a random seed
* `create_atoms 2 random 150 127569 region_inside_cylinder`: Randomly places 150 atoms of type 1 within the region inside the cylinder, using integer 127569 as a random seed

### 5) Run
* `write_data minimized_coordinate.data`: Outputs the final state of the system following energy minimization to file `minimized_coordinate.data` which can be used as starting point for subsequent simulations or analysis
