# Preparing The Water Reservoir

## Introduction

First, a rectangular water reservoir will be created & equilibrated in the absence of the PEG chain. The water model we will use is the Single Point Charge/Flexible Water (SPC/Fw) model is used <sup>1</sup>.

The SPC/Fw model is a refinement of the simple SPC model & Extended SPC (SPC/E) model where the water molecule is represented with three interaction sites, corresponding to two hydrogen atoms and the oxygen atom. In the simple SPC model, the intramolecular degrees of freedom are usually frozen to give rigid models, while the intermolecular interactions are described by Lennard-Jones & Coulombic potentials between sites with fixed-point charges. The main limitations of the simple SPC model are:
1. **The self-diffusion constant is poorly defined** - The self-diffusion constant is a measure of how fast water molecules move through the bulk liquid. The simple SPC model tends to underestimate the self-diffusion of water and therefore, water molecules tend to move slower than they do in reality in simulations that use this model
2. **The dielectric constant is poorly defined** - The dielectric constant is a measure of the amount of electric potential energy (in the form of induced polarization) is stored in a material. A dielectric is an insulating material, therefore the dielectric constant of an insulator measures the ability of the insulator to store electric energy in an electric field. In general, the lower a materials dielectric constant, the less polarizable it is. The simple SPC model does not account for the polarizability of water as it considers the charge distribution to be fixed. This fixed charge distribution therefore does not allow for the charge redistribution that results from polarization

The Extended SPC (SPC/E) model introduced a self-polarization energy correction term to account for the polarizability of the water molecules, and thus better define the dielectric constant. This was achieved by slightly increasing the atomic partial charges on the hydrogen & oxygen sites, but still does not account for intramolecular degrees of freedom/flexibility. However, because polarization effects are highly environment dependent, this simple approach significantly reduces the accuracy of the model in heterogeneous systems where the dielectric properties of water are influenced by that of another dielectric medium (e.g. membranes, proteins, ion channels etc.).

*Wu et al.* improved on this by investigating the effect of equilibrium bond length on the self-diffusion constant, and the effect of equilibrium bond angles on the dielectric constant. They found the sensitivity of both parameters to these molecular geometries was very high. Therefore, they derived a new flexible simple point-charge water model (the SPC/Fw model) by slightly perturbing the equilibrium bond length & angle to optimize the bulk water self-diffusion & dielectric constants towards their experimental values.

## Input Script Command Syntax



## References
[1] Wu, Y., Tepper, H.L. and Voth, G.A. (2006) ‘Flexible simple point-charge water model with improved liquid-state properties’, *The Journal of Chemical Physics*, 124(2).
