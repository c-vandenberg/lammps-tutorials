# Level 1 - Lennard-Jones Fluid first-input.lammps Script Explained

## PART A - ENERGY MINIMIZATION
This section of the script is concerned with minimizing the energy of the system to find a stable or lower energy configuration before beginning dynamics simulations

### 1) Initialization
* `units lj`: Specifies that the simulation uses Lennard-Jones (LJ) units. LJ units are dimensionless and are based on Lennard-Jones potential parameters:
  * Energies are expressed in units of **ϵ**. This is the well depth and is a measure of the strength of the particle-particle interaction
  * Distances are expressed in units of **σ**. This is the distance at which the particle-particle potential energy is zero
  * Masses are expressed in units of **m**. This is the atomic mass
* `dimension 3`: Sets the simulation in 3-dimensions
* `atom_style atomic`: Sets the atom style to atomic. This means that each atom is represented solely by its type and position, without additional features like charges or bonds. Therefore, each atom is treated as just a dot with a mass
* `pair_style lj/cut 2.5`: Defines the interaction of the particles as Lennard-Jones potential with a cut-off distance equal to r<sub>c</sub> = 2.5 (unitless). Therefore, for r < r<sub>c</sub>:

$$
V(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]
$$

* `boundary p p p`: 