# Level 1 - Carbon Nanotube Deformation Further Exercises: Carbon Nanotube Membrane

<p align="center">
  <img src="" alt="cnt-unbreakable-bonds start md" width="300" />
  <img src="" alt="cnt-unbreakable-bonds-mid md" width="300" />
  <img src="" alt="cnt-unbreakable-bonds-end md" width="300" />
</p>

## Problem

In this exercise we are to replicate the carbon nanotube (CNT) in the x and y directions and equilibrate the system to create an infinite membrane of multiple CNTs. We are also to apply a shear stress/sheer deformation to the CNT membrane.

## Solution

Firstly, to simulate the effect of an infinite system in the x and y directions, we must change the boundary conditions to periodic in the x and y directions using `boundary p p p`. These conditions make the simulation box behave as if it is a single unit cell in an infinite system.

To replicate the CNT we must use the `replicate x y z` command. Here we want to replicate the CNTs by 3 in x-direction and 3 in y-direction to give 9 in unit cell.

However, before we replicate we should adjust the size of the simulation box/unit cell so that the individual CNTs are packed closely by changing `change_box all x final -40 40 y final -40 40 z final -40 40` to `change_box all x final -7 7 y final -7 7` (this is mainly for visualization purposes).

In the `# 5) Run` section:
* We apply the NPT (constant **N**umber of particles, **P**ressure & **T**emperature) ensemble via command `fix npt_ensemble all npt temp ${T} ${T} 100 x 1 1 1000 y 1 1 1000`:
  * `fix npt_ensemble`: Defines a fix called `npt_ensemble`
  * `all`: Applies fix to all atoms
  * `npt`: Specifies that the fix applies an NPT ensemble
  * `temp ${T} ${T} 100`: Specifies the initial temperature (`${T}`), target temperature (`${T}`) and a thermostat coupling constant of `100` (time units over which the temperature is adjusted)
  * `x 1 1 1000`: Controls pressure in the x-direction:
    * `1 1`: The target pressure range (minimum, maximum) in the x-direction
    * `1000`: A barostat coupling constant (time units over which pressure is adjusted)
  * `y 1 1 1000`: Controls pressure in the y-direction
* We apply the NPT ensemble in order to equilibrate both the temperature & pressure of the system while allowing the simulation box to adjust in size:
  * The NPT ensemble allows the simulation box dimensions to change dynamically in response to the applied pressure
  * The simulation box therefore adjusts/scales along the x & y directions to accommodate the desired pressure, allowing the system to reach an equilibrium state
* We then stop the NPT ensemble via command `unfix npt_ensemble` and switch to the NVT (constant **N**umber of particles, **V**olume & **T**emperature) ensemble via `fix nvt_ensemble all nvt temp ${T} ${T} 100`
  * This is to ensure we maintain the temperature of the system while keeping the volume fixed after the system has equilibrated to the desired pressure via the NPT ensemble
  * Switching to the NVT ensemble ensures that the simulation box dimensions remain constant while we apply the deformation
* Finally, we apply a deformation to all atoms in the simulation via command `fix deformation all deform 1 xy erate 5e-5`:
  * `fix deformation`: Defines a fix called `deformation`
  * `all`: Applies fix to all atoms
  * `deform`: Specifies that the fix applies a deformation
  * `1 xy`: The fix is applied every 1 timesteps in both the x & y directions
  * `erate 5e-5`: Sets an **engineering rate** (`erate`) 5 * 10^-5 per timestep for the deformation
