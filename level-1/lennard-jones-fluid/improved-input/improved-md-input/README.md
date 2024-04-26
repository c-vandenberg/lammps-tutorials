# Level 1 - Lennard-Jones Fluid `improved_input.md.lammps` Script

This `improved-input.md.lammps` input script starts an MD simulation directly from the previously saved energy minimized system configuration in `improved-input.min.lammps`.

Once the energy minimized system state has been loaded, the script then removes any atoms that have migrated from one region to the other (i.e. from inside the cylinder to outside the cylinder, or vice versa) during minimization. 

During the MD simulation, as well as dumping all atoms per 1000 timesteps (similar to  it outputs the average number of atoms of each type in the cylinder as a function of time. It also outputs the coordination number per atom between atoms of type 1 and 2, i.e. the average number of atoms of type 2 in the vicinity of the atoms of type 1. Coordination number can be used as an indicator of the degree of mixing of the binary mixture. 

The conditions of the MD simulation (e.g. the NVE ensemble and Langevin thermostat) are the same as in the `first-input.lammps` input script.

## Input Script Command Syntax

The commands in the `5) Run` section are all the same as those in the `first-input.lammps` script and have already been explained. There are new commands in the `4) Visualization` section however:
* `variable number_atom_type_1_inside_cylinder equal count(group_atom_type_1,region_inside_cylinder)`: Sets variable `number_atom_type_1_inside_cylinder` equal to the number of type 1 atoms inside the cylinder region as a function of time. The next command does the same for the type 2 atoms. Both of these counts are useful for tracking how the population of different atom types changes within specific spatial boundaries during the simulation
* `fix fix_count_atom_type_1_inside_cylinder all ave/time 10 200 2000 v_number_atom_type_1_inside_cylinder &
    file atom-population/output_atom_type_1_population_inside_cylinder_vs_time.dat`: 
  * `fix fix_count_atom_type_1_inside_cylinder`: Defines a `fix` command with ID `fix_count_atom_type_1_inside_cylinder`
  * `all ave/time 10 200 2000 v_number_atom_type_1_inside_cylinder`: Sets the logic of the previous `compute` command. It calculates a time-averaged value for the previously defined `number_atom_type_1_inside_cylinder` variable. The `v_` prefix before `number_atom_type_1_inside_cylinder` refers to data that comes from that `variable` command. Every 10 timesteps the number of type 1 atoms in the cylinder is recorded, the data is averaged 200 times, and these averages are then written to the output file every 2000 timesteps
  * `file atom-population/output_atom_type_1_population_inside_cylinder_vs_time.dat`: This `file` command specifies where these time averaged values are outputted to
  * The next command does the same for the type 2 atoms
* `compute compute_atom_type_1_coordination_number group_atom_type_1 coord/atom cutoff 2.0 group group_atom_type_2`:
  * `compute compute_atom_type_1_coordination_number`: Defines a `compute` command with ID ``compute compute_atom_type_1_coordination_number``
  * `group_atom_type_1 coord/atom cutoff 2.0 group group_atom_type_2`: Calculates the coordination number of each type 1 atom with respect to type 2 atoms with a cutoff distance of 2 (unitless)
* `compute compute_average_atom_type_1_coordination_number all reduce ave c_compute_atom_type_1_coordination_number`:
  * `compute compute_average_atom_type_1_coordination_number`: Defines a `compute` command with ID `compute_average_atom_type_1_coordination_number`
  * `all reduce ave c_compute_atom_type_1_coordination_number`: Sets the logic of the previous `compute` command. It averages the coordination numbers computed in `compute_atom_type_1_coordination_number` over all atoms. This averaging by `compute ave` is necessary as the `compute coord/atom` command returns an array where each value corresponds to two given atoms i-j, and an array can’t be printed by `fix ave/time`
* `fix fix_average_atom_type_1_coordination_number all ave/time 10 200 2000 &
    c_compute_average_atom_type_1_coordination_number file atom-coordination-number/output_average_atom_type_1_coordination_number.dat`:
  * `fix fix_average_atom_type_1_coordination_number`: Defines a `fix` command with ID `fix_average_atom_type_1_coordination_number`
  * `all ave/time 10 200 2000 & c_compute_average_atom_type_1_coordination_number`: Sets the logic of the previous `fix` command. It calculates a time-averaged value for the previously calculated average coordination numbers in `average_atom_type_1_coordination_number`. The `c_` prefix before `average_atom_type_1_coordination_number` refers to data that comes from that `compute` command. Every 10 timesteps the average coordination number is recorded, the data is averaged 200 times, and these averages are then written to the output file every 2000 timesteps
  * `file atom-coordination-number/output_average_atom_type_1_coordination_number.dat`:This `file` command specifies where these time averaged values are outputted to