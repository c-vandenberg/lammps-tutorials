# Radial Distribution Function of Water Solvated PEG Polymer Molecule

## Introduction
The **radial distribution function (RDF)** is the most useful measure of the "structure" of a fluid at molecular length scales. Fluid here means any dense, disordered system which has local variation in the position of its constituent particles but is macroscopically isotropic (e.g. liquids & solutions, biomolecular systems, amorphous materials, crystalline solids etc.).

The RDF gives a statistical description of the local packing and particle density of the system, by describing the average particle distribution of particles around a central reference particle. Put simply, it is defined as the ratio of the local density of particles at a distance *r* from a reference particle, to the average particle density of the system.

The RDF is typically computed via these (simplified) steps:
1. **Select a Reference Particle** - Choose a particle as the reference point
2. **Measure Distances** - Measure the distances from the reference particle to all other particles within a cutoff radius
3. **Bin the Distances Into a Histogram** - Divide the distance range into smaller intervals (**bins**) and count the number of particles that fall within each bin. These bins are then plotted on a histogram (a graphical representation of the RDF. Usually an RDF-distance graph)
4. **Normalise** - Adjust the raw counts of particles in each bin dividing the total number of particles in each bin by (number of reference particles * volume of the spherical shell/bin * average number density of all particles in the system). This is normalization in the context of the RDF.

An RDF-distance graph/histogram will give a series of peaks that correspond to the first, second, third etc. coordination/solvation shell of the reference particle. Generally, the larger the RDF value for a given peak, the higher the particle density around the reference particle. Additionally, we can calculate the coordination number/average number of neighbours surrounding the reference particle if we integrate the peak/calculate the area under it.

In this exercise we will extract the RDF (or *g(r)*) between the oxygen atom of the water molecules and two oxygen atoms of the PEG molecule. We will compare the RDF before and after the force is applied to the PEG.

## Input Script Command Syntax

Breaking down the new commands we encounter in this input script:
```
# 4) Run
compute PEG_H2O_RDF all rdf 200 1 8 2 8 cutoff 10
fix ave_PEG_H2O_RDF all ave/time 10 4000 50000 c_PEG_H2O_RDF[*] file ave-time-output-data/ave_PEG_H2O_RDF_initial.dat mode vector
```
* `compute PEG_H2O_RDF all rdf 200 1 8 2 8 cutoff 10`:
  * `compute PEG_H2O_RDF all` - Specifies a `compute` command called `PEG_H2O_RDF` that applies to all atoms in the system
  * `200` - The number of bins to divide the distance range into when generating the histogram. This defines the resolution of the RDF
  * `1 8` - Specifies that the RDF should be calculated between atom type 1 (PEG oxygen) and atom type 8 (H<sub>2</sub>O oxygen)
  * `2 8` - Specifies that the RDF should be calculated between atom type 2 (PEG oxygen) and atom type 8 (H<sub>2</sub>O oxygen)
  * `cutoff 10` - Sets the maximum distance for which the RDF is calculated as 10 Å
* `fix ave_PEG_H2O_RDF all ave/time 10 4000 50000 c_PEG_H2O_RDF[*] file ave-time-output-data/ave_PEG_H2O_RDF_initial.dat mode vector`:
  * `fix ave_PEG_H2O_RDF all` - Specifies a `fix` command called `ave_PEG_H2O_RDF` that applies to all atoms in the system
  * `ave/time` - Specifies that this fix will perform time-averaging of a specified quantity
  * `10` - Number of timesteps between each sampling 
  * `4000` - Number of samples to average over before resetting
  * `50000` - Total number of timesteps over which the averaging will be performed
  * `c_PEG_H2O_RDF[*]` - Refers to the output of the previously defined compute named PEG_H2O_RDF. The [*] indicates that all components of this compute are included.
  * `file ave-time-output-data/ave_PEG_H2O_RDF_initial.dat` - The file where the averaged RDF data will be output
  * `mode vector` - Specifies that the output is a vector.