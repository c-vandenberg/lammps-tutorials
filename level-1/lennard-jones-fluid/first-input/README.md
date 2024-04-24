# Level 1 - Lennard-Jones Fluid first-input.lammps Script

This `first-input.lammps` input script creates a simple cubic 3-dimensional simulation box, with periodic boundary conditions (i.e. the simulation box is conceptually replicated infinitely in 3-dimensions), that uses Lennard-Jones units.

It also populates this simulation box with two types of atomic style atoms with the same mass but different energy and distance parameters. It also specifies a cut-off distance for any Lennard-Jones potential interactions.

Finally, it performs an energy minimization of the system, outputting basic thermodynamic properties, before performing an MD simulation.

## Input Script Syntax

### PART A - ENERGY MINIMIZATION
This section of the script is concerned with minimizing the energy of the system to find a stable or lower energy configuration before beginning dynamics simulations

### 1) Initialization
* `units lj`: Specifies that the simulation uses Lennard-Jones (LJ) units. LJ units are dimensionless and are based on Lennard-Jones potential parameters:
  * Energies are expressed in units of **ϵ**. This is the well depth and is a measure of the strength of the particle-particle interaction
  * Distances are expressed in units of **σ**. This is the distance at which the particle-particle potential energy is zero
  * Masses are expressed in units of **m**. This is the atomic mass
* `dimension 3`: Sets the simulation in 3-dimensions
* `atom_style atomic`: Sets the atom style to atomic. This means that each atom is represented solely by its type and position, without additional features like charges or bonds. Therefore, each atom is treated as just a dot with a mass
* `pair_style lj/cut 2.5`: Defines the interaction of the particles as Lennard-Jones potential with a cut-off distance equal to r<sub>c</sub> = 2.5 (unitless). Therefore, for r < r<sub>c</sub>:


* `boundary p p p`: Sets periodic boundary conditions along all three directions of space (the 3 p stand for x, y, and z respectively)

### 2) System Definition
* `region simulation_box block -20 20 -20 20 -20 20`: Defines a cubic simulation region called `simulation_box` that extends from -20 to 20 (unitless) along all 3 directions of space
* `create_box 2 simulation_box`: Creates the pre-defined `simulation_box` block that contains two types of atoms
* `create_atoms 1 random 1500 341341 simulation_box`: Creates 1500 atoms of type 1 randomly within the region `simulation_box`. The integer 341341 is a random seed that can be changed in order to create different initial conditions for the simulation
* `create_atoms 2 random 100 127569 simulation_box`: Creates 100 atoms of type 1 randomly within the region `simulation_box` with random seed 127569

### 3) Simulation settings
* `mass 1 1`: Attribute mass, **m** = 1 (unitless) to atoms of type 1
* `mass 2 1`: Attribute mass, **m** = 1 (unitless) to atoms of type 2
* `pair_coeff 1 1 1.0 1.0`: Defines the Lennard-Jones interaction parameters/coefficients for the interactions between atoms of type 1 (respectively the energy parameter **ϵ<sub>11</sub>** = 1.0 and the distance parameter **σ<sub>11</sub>** = 1.0)
* `pair_coeff 2 2 0.5 3.0`: Defines the Lennard-Jones interaction parameters/coefficients for the interactions between atoms of type 2 (respectively the energy parameter **ϵ<sub>22</sub>** = 0.5 and the distance parameter **σ<sub>22</sub>** = 3.0)

### 4) Visualization
* `thermo 10`: Specifies that LAMMPS print thermodynamic output (e.g. temperature, energy) in the terminal every 10 steps
* `thermo_style custom step temp pe ke etotal press`: Configures custom thermodynamic output to simulation step number(`step`), temperature (`temp`), potential energy (`pe`), kinetic energy(`ke`), total energy(`etotal`), and pressure (`press`)

### 5) Run
* `minimize 1.0e-4 1.0e-6 1000 10000`: Specifies LAMMPS to perform an energy minimization. By default, LAMMPS uses the conjugate gradient (CG) algorithm. This involves adjusting the coordinates of the atoms that are too close to each other until one of the four stopping criteria is reached. The four stopping criteria are:
  * **Energy tolerance** - The change in energy between the two iterations. Here it is < 1.0e-4
  * **Force tolerance** - The maximum force between two atoms in the system. Here it is < 1.0e-6
  * **Maximum Number of Iterations** - The maximum iterations in the energy minimization algorithm. Here it equals 1000
  * **Maximum Number of Force Evaluations** - The maximum number of times the force and energy are evaluated. Here it equals 10,000

### PART B - MOLECULAR DYNAMICS
This section of the script is concerned with performing the actual molecular dynamics simulation that will start from the final state of the energy minimization.

Molecular dynamics (MD) is based on the numerical solution of the Newtonian equations of motion for every atom `i`:
* ∑<sub>ij</sub> F<sub>ij</sub> = m<sub>i</sub> * a<sub>i</sub>
* Where:
  * ∑ = Summation of all atoms other than `i`
  * F<sub>ij</sub> = Force between the atom pairs `i-j`
  * m<sub>i</sub> = Mass of atom `i`
  * a<sub>i</sub> = Acceleration of atom `i`

The Newtonian equations are solved every timestep to predict the evolution of the positions and velocities of atoms and molecules over time.

At every timestep of an MD simulation, the following operations occur:
1. The forces between the atoms are calculated from the potential energy (here the Lennard-Jones potential)
2. The acceleration of each atom is evaluated from the Newtonian equation
3. The velocity and position of each atom is updated according to the calculated acceleration (typically using the Verlet algorithm, or similar)

### 1) Visualization
* `thermo 50`: Specifies that LAMMPS print thermodynamic output (e.g. temperature, energy) in the terminal every 10 steps
* `dump dump_all_atoms_per_100_timestep all atom 100 dump_all_atoms_per_100_timestep.lammpstrj`: A `dump` command configures the output of atomic positions and other specified data to a file for visualization for further analysis. Breaking this full command down:
  * `dump_all_atoms_per_100_timestep`: Custom identifier for this particular `dump` command
  * `all`: Specifies which atoms are to be included in the dump file. In this case all atoms in the simulation are to be included
  * `atom`: Defines the style of the dump, in this case "atom". This style outputs properties of individual atoms, typically their positions, velocities & other attributes depending on the atom style used in the simulation
  * `100`: The frequency (in timesteps) at which the data is written to the file. In this case, every 100 timesteps
  * `dump_all_atoms_per_100_timestep.lammpstrj`: The name of the output file where the data will be written. The `.lammpstrj` extension is a file that can be read by MD simulation visualization tools such as VMD or OVITO

### 2) Run
* `fix fix_nve_ensemble all nve`: A `fix` command is any operation that is applied to the system during timestepping or minimization. Fixes perform their operations at different stages of the timestep, in the order they are specified in the input script. Examples include updating atom positions & velocities via time integration, controlling temperature, applying constraints forces to atoms etc. Breaking this full command down:
  * `fix_nve_ensemble`: Custom identifier for this particular `fix` command
  * `all`: Specifies that the fix is applied to all atoms in the simulation
  * `nve`: Applies the NVE ensemble to the simulation, which stands for "constant Number of particles (**N**), Volume (**V**), and Energy(**E**)". It simulates the natural evolution of an isolated system consistent with the microcanonical ensemble where no energy or particles are exchanged with the environment. Therefore, all microstates whose energy falls within a given range have equal probability, and those outside that range are given a probability of zero
* `fix fix_langevin_thermostat all langevin 1.0 1.0 0.1 1530917`: Breaking this full command down:
  * `fix_langevin_thermostat`: Custom identifier for this particular `fix` command
  * `all`: Specifies that the fix is applied to all atoms in the simulation
  * `langevin`: Applies a Langevin thermostat to the simulation. A Langevin thermostat maintains the temperature of the simulation through a modification of Newton's equations of motion
  * `1.0 1.0`: The target temperatures (unitless) at the start and end of the simulation respectively
  * `0.1`: A dampening factor, which controls the rate of energy exchange between the system ad thermostat. The smaller the value the faster the energy exchange
  * `1530917`: Integer 1530917 is a random seed for the stochastic (random probability distribution), ensuring reproducibility of the simulation's randomness
* `timestep 0.005`: Sets the simulation timestep size, 0.005 (unitless)
* `run 10000`: Run MD simulation for specified number of timesteps, 10,000 timesteps in this case

## Input Script Syntax
