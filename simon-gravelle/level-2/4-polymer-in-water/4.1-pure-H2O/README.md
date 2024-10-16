# 4.1 Preparing The Water Reservoir

<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/dcbb5fe3-9e4e-4de5-a64a-931e866c2523" alt="pure H2O" width="" />
</div>

## 4.1.1 Introduction

First, a rectangular water reservoir will be created & equilibrated in the absence of the PEG chain. The water model we will use is the Single Point Charge/Flexible Water (SPC/Fw) model<sup>1</sup>.

The SPC/Fw model is a refinement of the simple SPC model & Extended SPC (SPC/E) model where the water molecule is represented with three interaction sites, corresponding to two hydrogen atoms and the oxygen atom. 

In the simple SPC model, the intramolecular degrees of freedom are usually frozen to give rigid models, while the intermolecular interactions are described by Lennard-Jones & Coulombic potentials between sites with fixed-point charges. The main limitations of the simple SPC model are:
1. **The self-diffusion constant is poorly defined**:
   * The self-diffusion constant is a measure of how fast water molecules move through the bulk liquid.
   * The simple SPC model tends to underestimate the self-diffusion of water and therefore, water molecules tend to move slower than they do in reality in simulations that use this model.
3. **The dielectric constant is poorly defined**:
   * The dielectric constant is a measure of the amount of electric potential energy (in the form of induced polarization) is stored in a material.
   * A dielectric is an insulating material, therefore the dielectric constant of an insulator measures the ability of the insulator to store electric energy in an electric field.
   * In general, the lower a materials dielectric constant, the less polarizable it is.
   * The simple SPC model does not account for the polarizability of water as it considers the charge distribution to be fixed. This fixed charge distribution therefore does not allow for the charge redistribution that results from polarization.

The Extended SPC (SPC/E) model introduced a self-polarization energy correction term to account for the polarizability of the water molecules, and thus better define the dielectric constant. This was achieved by slightly increasing the atomic partial charges on the hydrogen & oxygen sites, but still does not account for intramolecular degrees of freedom/flexibility. However, because polarization effects are highly environment dependent, the accuracy of this simple approach is significantly limited in heterogeneous systems where the dielectric properties of water are influenced by that of another dielectric medium (e.g. membranes, proteins, ion channels etc.).

*Wu et al.* improved on this by investigating the effect of equilibrium bond length on the self-diffusion constant, and the effect of equilibrium bond angles on the dielectric constant. They found the sensitivity of both parameters to these molecular geometries was very high. Therefore, they derived a new flexible simple point-charge water model (the SPC/Fw model) by slightly perturbing the equilibrium bond length & angle to optimize the bulk water self-diffusion & dielectric constants towards their experimental values.

## 4.1.2 Data Analysis

### Pure H<sub>2</sub>O Temperature as a Function of Time

<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/60efcd37-386e-4d96-94cb-9d815ef2a92f" alt="pure H2O temp vs time" width="" />
</div>

During equilibration of pure H<sub>2</sub> system, system temperature rapidly rises and plateaus at 300K as defined in the LAMMPS input script. Equilibrium temperature is reached in ~ 0.8 ps.

### Pure H<sub>2</sub>O Volume as a Function of Time

<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/c9917ab1-e1e5-43bc-8dc1-68109147a59b" alt="pure H2O vol vs time" width="" />
</div>

The NPT ensemble applies isotropic pressure control meaning that the simulation box expands or contracts equally in all directions to control pressure. Once its applied during equilibration of pure H<sub>2</sub> system, system volume decreases gradually, reaching an equilibrium volume of ~ 34 nm<sup>3</sup> in ~ 11 ps.

### Pure H<sub>2</sub>O Density as a Function of Time

<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/be837b5f-c771-4cff-8848-ba12baf20d17" alt="pure H2O density vs time" width="" />
</div>

Due to the inverse relationship with volume, the system density increases gradually as system volume decreases gradually decreases. The equilibrium density of ~ 31 nm<sup>-3</sup> is reached in ~ 13 ps.

## 4.1.3 Input Script Command Syntax

Breaking down the new commands we encounter in this input script:

```
# 1) Initialization
units real
atom_style full
bond_style harmonic
angle_style harmonic
dihedral_style harmonic
pair_style lj/cut/coul/long 12
kspace_style pppm 1e-5
special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes
```
* `pair_style lj/cut/coul/long 12` - Specifies that atoms interact via both Lennard-Jones potential & Coulombic interactions, with a cut-off of 12 Å
  * The cutoff of 12 Å applies to both Lennard-Jones and Coulombic interactions, slightly differently for each interaction
  * For Lennard-Jones interactions (`lj/cut` here), atoms interact with each other only if they are separated by a distance smaller than the cutoff
  * For Coulombic interactions (`coul/long` here), interactions between atoms closer than the cutoff are computed directly, and interactions between atoms outside that cutoff are computed in the reciprocal space
* `kspace_style pppm 1e-5` - Defines the **K-space long range solver** used to compute long-range Coulombic interactions as PPPM (Particle-Particle Particle-Mesh) with an accuracy of 1e-5
  * Derived by *Luty et al.*<sup>2</sup>, the PPPM method is based on separating the total interaction between particles into the sum of short-range interactions (computed by direct particle-particle summation), and long-range interactions (calculated by solving Poisson's equation using periodic boundary conditions)
* `special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 1.0 angle yes` - Modifies the non-bonded interactions (Lennard-Jones & Coulombic) between atoms connected by bonds. This command scales:
  * **1-2 Lennard-Jones interactions** (i.e. atoms directly bonded) by **0.0**, ignoring these interactions
  * **1-3 Lennard-Jones interactions** (i.e. atoms two bonds apart/angle interactions) by **0.0**, ignoring these interactions
  * **1-4 Lennard-Jones interactions** (i.e. atoms three bonds apart/dihedral interactions) by **0.5**, halving the calculated interaction value
  * **1-2 Coulombic interactions** (i.e. atoms directly bonded) by **0.0**, ignoring these interactions
  * **1-3 Coulombic interactions** (i.e. atoms two bonds apart/angle interactions) by **0.0**, ignoring these interactions
  * **1-4 Coulombic interactions** (i.e. atoms three bonds apart/dihedral interactions) by **1.0**, therefore are not altered/scaled
  * `angle yes` - Specifies that the 1-3 interactions weighting factor will be ignored if the atoms are not listed as the first and last atoms in any angle defined in the simulation

```
# 2) System Definition
region box block -15 15 -15 15 -15 15
create_box 9 box &
bond/types 7 &
angle/types 8 &
dihedral/types 4 &
extra/bond/per/atom 3 &
extra/angle/per/atom 6 &
extra/dihedral/per/atom 10 &
extra/special/per/atom 14
```
* `extra/bond/per/atom 3` - Used to allocated additional memory for bonds on a per-atom basis, so that enough is allocated for bonds that are created during the simulation. Here, additional memory is allocated for `3` additional bonds for each atom
* `extra/angle/per/atom 6` - Allocates additional memory for `6` angles for each atom
* `extra/dihedral/per/atom 20` - Allocates additional memory for `10` dihedrals for each atom
* `extra/special/per/atom 14` - Allocates additional memory for `14` special pairwise interactions for each atom

## 4.1.4 References
[1] Wu, Y., Tepper, H.L. and Voth, G.A. (2006) ‘Flexible simple point-charge water model with improved liquid-state properties’, *The Journal of Chemical Physics*, 124(2). <br>
[2] Luty, B.A. and van Gunsteren, W.F. (1996) ‘Calculating electrostatic interactions using the particle−particle particle−mesh method with nonperiodic long-range interactions’, *The Journal of Physical Chemistry*, 100(7), pp. 2581–2587.
