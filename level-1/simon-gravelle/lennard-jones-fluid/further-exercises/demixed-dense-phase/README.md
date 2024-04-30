# Level 1 - Lennard-Jones Fluid Further Exercises: Create a De-mixed Dense Phase

## Problem
Using one of the input scripts from either `first-input` or `improved-input`, fine-tune parameters such as atom numbers and atomic interaction to create a simulation with the following properties:
* System with high atomic density
* Both type 1 and type 2 atoms must be the same size
* The type 1 and type 2 atoms must de-mix

**Hint**
* An easy way to create a dense phase is to adjust the simulation box dimensions. This can be achieved by using `fix nve` instead of `fix nph`, which adjusts simulation box dimensions to control pressure

## Solution

### De-mixing
The key to creating de-mixing phase is to **adjust the Lennard-Jones parameters**:

```
# 3) Simulation settings
mass 1 1
mass 2 1
pair_coeff 1 1 5.0 1.0
pair_coeff 2 2 5.0 1.0
pair_coeff 1 2 0.05 1.0
```

* Atoms of both type 1 and type 2 have the same distance parameter (`σ`) value of 1.0, meaning they both have the same diameter
* Additionally, atoms of both type 1 and type 2 have a large energy parameter (`ϵ`) value of 5.0, meaning self-interaction is highly favourable
* Finally, the small energy parameter (`ϵ`) value of 0.05 for interaction between atoms of type 1 and atoms of type 2 means that interaction between the two is highly disfavoured

### High Density
As stated before, to create a dense phase we can use `fix nph`, which will adjust the simulation box dimensions in order to control pressure:
* `fix fix_nph_ensemble all nph iso 1.0 1.0 0.5`:
  * `fix fix_nph_ensemble`: Defines a `fix` command with ID `fix_nph_ensemble`
  * `all nph`: Applies the **isenthalpic (nph) ensemble** to all atoms in the simulation. The NPH ensemble is used to perform MD simulations under a constant **N**umber of particles, **P**ressure and enthalpy (**H**)
  * `iso`: Applies **isotropic** pressure control, meaning that the pressure is adjusted uniformly in all directions
  * `1.0 1.0 0.5`: Sets the starting pressure to `1.0`, the end pressure to `1.0` and a dampening factor of `0.5` (dampening factor controls the rate at which the barostat adjusts the system pressure to the target pressure)