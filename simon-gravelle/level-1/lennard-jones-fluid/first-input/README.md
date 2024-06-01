# Level 1 - Lennard-Jones Fluid `first-input.lammps` Script

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/e65aa729-dc73-40f7-85d8-2d8445dbd537" alt="first-input" width="500">
</p>

## Exercise

This `first-input.lammps` input script creates a simple cubic 3-dimensional simulation box, with periodic boundary conditions (i.e. the simulation box is conceptually replicated infinitely in 3-dimensions), that uses Lennard-Jones units.

It also populates this simulation box with two types of atomic style atoms with the same mass but different energy and distance parameters. It also specifies a cut-off distance for any Lennard-Jones potential interactions.

Finally, it performs an energy minimization of the system, outputting basic thermodynamic properties, before performing an MD simulation.

## Data Analysis

### Lennard-Jones Potential as a Function of Interatomic Distance

In the `3) Simulation settings` section of the LAMMPS input script, the line `pair_coeff 1 1 1.0 1.0` defines the Lennard-Jones interaction parameters/coefficients for interactions between atoms of type 1. Here, the Lennard-Jones energy parameter **ϵ<sub>11</sub>** = 1.0 and the Lennard-Jones distance parameter **σ<sub>11</sub>** = 1.0.

Similarly, the line `pair_coeff 2 2 0.5 3.0` defines the Lennard-Jones interaction parameters/coefficients for interactions between atoms of type 2. Here, the Lennard-Jones energy parameter **ϵ<sub>11</sub>** = 0.5 and the Lennard-Jones distance parameter **σ<sub>11</sub>** = 3.0.

For interactions between atoms of type 1 and atoms of type 2, by default LAMMPS calculates the cross coefficients  **ϵ<sub>12</sub>** and **σ<sub>12</sub>** using the geometric average:

<br>
<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cepsilon_%7Bij%7D%20%3D%20%5Csqrt%7B%5Cepsilon_%7Bii%7D%5Cepsilon_%7Bjj%7D%7D", alt='epsilon-geometric-ave-equation'/>
</div>
<br>
<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Csigma_%7Bij%7D%20%3D%20%5Csqrt%7B%5Csigma_%7Bii%7D%5Csigma_%7Bjj%7D%7D", alt='sigma-geometric-ave-equation'/>
</div>
<br>

Therefore, LAMMPS will have calculated **ϵ<sub>12</sub>** and **σ<sub>12</sub>** for us as:

<br>
<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Cepsilon_%7B0.5%7D%20%3D%20%5Csqrt%7B%5Cepsilon_%7B1%7D%5Cepsilon_%7B0.5%7D%7D%20%3D%200.707", 
    alt='epsilon-atom-1-atom-2-geometric-ave-equation'/>
</div>
<br>
<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%5Csigma_%7B1.0%7D%20%3D%20%5Csqrt%7B%5Csigma_%7B1.0%7D%5Csigma_%7B3.0%7D%7D%20%3D%201.732", 
    alt='epsilon-atom-1-atom-2-geometric-ave-equation'/>
</div>
<br>

These Lennard-Jones parameters/coefficients can then be used to calculate the Lennard-Jones potential **E<sub>ij</sub>** for all three interactions (**E<sub>11</sub>**, **E<sub>22</sub>** and **E<sub>12</sub>**) via the equation:

<br>
<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?\color{white}V(r)%20=%204\epsilon%20\left[%20\left(\frac{\sigma}{r}\right)^{12}%20-%20\left(\frac{\sigma}{r}\right)^{6}%20\right]" 
    alt="lennard-jones-potential-equation"/>
</div>

where:
- *V(r)* is the potential energy as a function of interatomic distance *r*,
- *ε* = the depth of the potential well (the deeper the well, the stronger the interatomic interaction)
- *σ* = the finite distance at which the interatomic potential is zero,
- *r* = the interatomic distance,
- ![sigma_r_12](https://latex.codecogs.com/svg.latex?\color{white}\left(\frac{\sigma}{r}\right)^{12}) = the replusive force
- ![sigma_r_6](https://latex.codecogs.com/svg.latex?\color{white}\left(\frac{\sigma}{r}\right)^{6}) =  the attractive force

These calculated Lennard-Jones potentials as a function of interatomic distance are shown below:

<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/134a20a6-d02a-4185-9a2a-20b8da8c6e06" alt="lennard_jones_potential" width="">
</div>

The graph shows that as the interatomic distance decreases, the Lennard-Jones/bonding potential energy decreases and so the probability of interaction increases. As the two atoms come closer they eventually reach a distance region where they become bound; their Lennard-Jones potential/bonding potential energy decreases from zero to a negative quantity. 

While the atoms are bound, the interatomic distance continues to decrease until the particles reach an equilibrium, specified by the separation distance at which the minimum potential energy is reached (the minimum/depth of the potential well *ε*). 

The more negative the minimum potential energy/the deeper the well depth, the greater the interaction between the two atoms. As you can see, the self-interaction between atom type 1 atoms is the strongest, followed by the self-interaction between atom type 2 atoms, with atom type 1-atom type 2 interactions being the weakest.

However the two bound atoms are further pressed together, past their equilibrium distance, repulsion begins to occur. This is a situtation where the atoms are so close together that their electrons are forced to occupy each other’s orbitals. This electron-electron replusion causes the Lennard-Jones/bonding potential energy to increase rapidly as the distance of separation decreases (highly unfavourable). This illustrates the high exponent of the repulsive force in the Lennard-Jones potential equation:

<div align="center">
  <img src="https://latex.codecogs.com/svg.latex?\color{white}\left(\frac{\sigma}{r}\right)^{12}" 
    alt="lennard-jones-potential-equation-repulsive-component"/>
</div>

## Input Script Command Syntax

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
  * `0.1`: A dampening factor, which controls the rate of energy exchange between the system and the thermostat. The smaller the value the faster the energy exchange
  * `1530917`: Integer 1530917 is a random seed for the stochastic (random probability distribution), ensuring reproducibility of the simulation's randomness if the same random seed is used
* `timestep 0.005`: Sets the simulation timestep size, 0.005 (unitless)
* `run 10000`: Run MD simulation for specified number of timesteps, 10,000 timesteps in this case
