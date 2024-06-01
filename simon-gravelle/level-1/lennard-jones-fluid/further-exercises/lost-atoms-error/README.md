# Level 1 - Lennard-Jones Fluid Further Exercises: Lost Atoms Error

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/bcdf929d-c941-45fd-9ded-7d1945a0e104" alt="lost-atoms-error" width="500" />
</p>

## Problem
For this exercise, the following input script is provided:

```
# 1) Initialization
units lj
dimension 3
atom_style atomic
pair_style lj/cut 2.5
boundary p p p

# 2) System definition
region simulation_box block -20 20 -20 20 -20 20
create_box 1 simulation_box
create_atoms 1 random 1000 341841 simulation_box

# 3) Simulation settings
mass 1 1
pair_coeff 1 1 1.0 1.0

# 4) Visualization
dump dump_all_atoms_per_100_timestep all atom 100 dump_all_atoms_per_100_timestep.lammpstrj
thermo 100
thermo_style custom step temp pe ke etotal press

# 5) Run
fix fix_nve_ensemble all nve
fix fix_langevin_thermostat all langevin 1.0 1.0 0.1 1530917
timestep 0.005

run 10000
```

In its current state, this input script returns the following error:

```
ERROR: Lost atoms: original 1000 current 984
```

This is one of the most common errors encountered when using LAMMPS. The goal of this exercise is to fix this `Lost atoms` error without using any other commands than the ones already present. The values of the parameters can be altered and/or each command can be replicated as many times as needed.

**Hint**
* This script is failing because atoms are created randomly in space, some of them are likely overlapping, and no energy minimization is performed prior to start the molecular dynamics simulation.

## Solution
The `Lost atoms` error typically indicates that some of the atoms are moving outside the simulation box. This occurs because, in the absence of energy minimization, some of the randomly created atoms can be highly overlapping/too close to each other, and are immediately subjected to strong repulsive forces. As a result, these atoms have unreasonably high initial velocities that lead to them being displaced beyond the boundary the simulation domain.

The solution is to **reduce the initial timestep value** as well as **reducing the temperature imposed by the Langevin thermostat**. Additionally, in order to make sure that the temperature of the system quickly reaches a reasonable value, the **dampening parameter of the Langevin thermostat needs to be reduced**:

```
# 5) Run
fix fix_nve_ensemble all nve
fix fix_langevin_thermostat all langevin 0.001 0.001 0.001 1530917
timestep 0.0001
run 10000
```

The reasons why this is a valid solution are explained below.

### Reducing the Initial Temperature
* The `fix langevin` command is used to maintain the system at a target temperature by adding random thermal noise and a frictional drag force
* Reducing this temperature from `1.0` to `0.001` greatly reduces the initial kinetic energy of the atoms
* This reduction in kinetic energy results in the atoms having much smaller initial velocities
* This reduces the likelihood of their displacement beyond the boundary of the simulation domain when subjected to the high repulsive forces that result from any initial overlapping positions

### Reducing the Timestep
* Reducing the timestep to `0.0001` is a crucial adjustment for maintaining stability, especially when the starting conditions/interactions may lead to high accelerations
* A smaller timestep allows the time integrator (in this case `fix fix_nve_ensemble`) to integrate Newton's equations of motion more frequently, allowing it to more accurately capture the rapidly changing dynamics of the system. 
* With a smaller timestep, the integration is more precise and less likely to result in the errors that lead to lost atoms

### Changing the Dampening Parameter
* The damping parameter in the Langevin thermostat controls the rate of energy exchange between the system and the thermostat. It therefore determines how quickly the system's kinetic energy is adjusted to match the target temperature
* A smaller dampening parameter (from `0.1` to `0.001`) means the thermostat will act more aggressively to bring the system temperature to the desired value
* This will counter any drastic fluctuations in kinetic energy resulting from initial overlaps/high repulsive forces and thus dampens out any excessive velocities



