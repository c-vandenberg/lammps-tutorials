# Level 1 - Carbon Nanotube Deformation Further Exercises: Flying Ice Cube Artifact

<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/179c4d67-5d97-4c65-a655-f4f29386cb25" alt="cnt-flying-ice-cube-artifact" width="300" />
</p>

## Introduction
The "flying ice cube" effect is one of the most famous artifacts (unintended behaviours)/phenomena of molecular simulations. It is where the simulated system, typically a box containing molecules or atoms, begins to exhibit a collective motion in one direction. The effect is named metaphorically; just as an ice cube might slide across a surface, the entire group of atoms/molecules in the simulation "fly" together, moving as a whole through the simulation box.

### What causes the Flying Ice Cube Effect?
1. **Issues with Conservation of Momentum**
    * According to classical physics, in a perfectly isolated system, total momentum should be conserved (i.e. remains constant)
    * In MD simulations however, if the initial conditions or subsequent manipulations (e.g. removing or modifying atoms/molecules) don't correctly account for overall momentum, it can lead to the whole system gaining net momentum in a particular direction
2. **Thermostat Artifacts**
   * Some methods of temperature control (thermostats) used in MD simulations can inadvertently impart a drift to the center of mass of the system
   * This occurs when the thermostat algorithm might apply corrections uniformly to atom/molecule velocities without necessarily zeroing out the total momentum of the system. 
   * If not carefully implemented, these adjustments can accumulate, leading to a net momentum
3. **Numerical Instabilities & Errors**
   * Accumulation of small numerical errors over long simulation times can also contribute to this effect, especially in large, complex systems or simulations run over very long timescales

### Why is the Flying Ice Cube Effect a Problem?
1. **Non-Physical Behaviour**
   * The flying ice cube effect is generally considered undesirable because it introduces a non-physical artifact into the simulation
   * In reality, there should be no external force imparting a net momentum to the entire system, so this effect can distort the dynamics and properties studied
2. **Interference with Boundary Conditions**
   * In simulations using periodic boundary conditions (where atoms/molecules exiting one side of the simulation box re-enter from the opposite side), the flying ice cube effect can lead to incorrect interactions and artifacts as the whole system cyclically exits and enters the boundaries

### How can the Flying Ice Cube Effect Be Fixed?
1. **Zeroing Total Momentum**
   * A common approach to prevent or correct for this effect is to periodically adjust the velocities of all atoms/molecules to zero out the total momentum of the system
   * This adjustment ensures that no net translational movement affects the system's dynamics
   * An example a command that would zero momentum is `velocity all zero linear`:
     * `velocity` - Command used to set or modify the velocities of a group of atoms
     * `all` - Specifies that the command should apply to all atoms/molecules in the simulations. Can be replaced to target specific atom/molecule groups if needed
     * `zero liner` - Zeroes out the linear (translational) momentum of the system, ensuring that the sum of the velocities, weighted by their masses, is zero in the x, y & z directions
   * Another example of a command would be `fix zero_momentum all momentum 100 linear 1 1 1`:
     * `fix zero_momentum` - Defines fix `zero_momentum`
     * `all` - Specifies that the command should apply to all atoms/molecules in the simulations. Can be replaced to target specific atom/molecule groups if needed
     * `momentum` - Specifies the fix style used to control the momentum of the system
     * `100` - The frequency (in timesteps) at which the fix operation is to be applied. Therefore, it zeroes the system's momentum every 100 timesteps
     * `linear` - Specifies that the linear (translational) momentum should is to be adjusted/zeroed
2. **Careful Use of Thermostats**
   * Choosing a thermostat that conserves or does not affect the total momentum, or configuring them to properly adjust for any changes in particle velocities that they apply can also prevent or correct the effect
3. **Monitoring Simulation Parameters**
   * Regular checks on the kinetic energy, temperature & center of mass movement of the system can alert you to the development of any non-physical behaviours early in the simulation process

## Problem
In this exercise, we are given a new carbon nanotube molecular data file (`flying_ice_cube_cnt_molecular.data`) and a seemingly simple input script: 

```
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

read_data cnt_molecular.data
include ../../unbreakable-bonds/parm.lammps

group carbon_atoms type 1

fix mynve all nve
fix myber all temp/berendsen ${T} ${T} 100
fix myrct all recenter 0 0 0

thermo 1000
dump mydmp all atom 1000 dump.lammpstrj

timestep 1.0
run 10000
```

When the simulation is run, we see that the temperature is very close to 300 K, as expected. However, if we look at the system using VMD, we can see that the atoms are not moving.

The task is to identify the origin of the issue and fix the input.

## Solution
If we look in the `flying_ice_cube_cnt_molecular.data` file we can see that we have large velocities in the x-direction, which accumulates and causes a net momentum in the x-direction:
```
Velocities

24 0.007983439029626362 -6.613056392124822e-05 7.867644943646289e-05
1 0.007906200203484036 3.252025147011299e-05 -4.4209216231039336e-05
25 0.007861090484107148 9.95045322688365e-06 -0.00014277147407215768
(...
```

If the system has a high net momentum as it does here, the kinetic energy derived from this momentum makes it appear as if the system has a higher temperature. This is because the kinetic energy used to calculate the temperature includes both the relative motion of the atoms within the system, and the motion of the system as a whole. This causes the Berendsen thermostat to fail.

The Berendsen thermostat, attempting to control the temperature, rescales the velocities based on what it calculates as temperature (which as we said, is influenced heavily by the systems high net momentum in the x-direction). This temperature calculation is artificially high due to the system's high net momentum, and so the thermostat improperly scales down the velocities, attempting to reduce what it perceives as excess thermal energy.

As a result, despite the high momentum, the internal velocities of the atoms are scaled down too much, making the system appear frozen or much less dynamic than expected. This is because the thermostat addresses the wrong aspect of the system's kinetic energy (the part due to net momentum, not internal thermal energy)

There are several ways to fix this:
1. **Zeroing Out Net Momentum** - Using a command like `fix zero_momentum all momentum 100 linear 1 1 1` (used in the script) or `velocity all zero linear` after reading the initial velocities but before temperature control will adjust the velocities so that the net momentum in the x-direction does not accumulate
2. **Use a More Robust Thermostat** - Switching to a thermostat like **Nosé-Hoover** is also a solution. Nosé-Hoover thermostat is better at handling systems with net momentum because it adjusts velocities in a manner that more accurately maintains the canonical ensemble without being misled by the total kinetic energy
