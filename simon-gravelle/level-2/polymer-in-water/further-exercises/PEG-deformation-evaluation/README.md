# Evaluate The Deformation of The PEG Polymer Molecule

## Introduction
Once the PEG molecule is fully stretched, its structure differs from when it is not stretched. We can probe this deformation by extracting typical intra-molecular parameters such as dihedral angles.

In this exercise we will extract the histograms of the angular distribution of the PEG dihedral angles in before and after stretching.

## Data Analysis
<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/d9a6a492-4261-421b-bf31-bdfa2cdbf9e9" alt ="PEG deformation probability vs dihedrals" width="" />
</div>

A dihedral angle *ϕ* = 60° & *ϕ* = 180° corresponds to the gauche conformation, whereas the a dihedral angle *ϕ* = 60° represents an eclipsed conformation

<div align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/599d314a-33d4-4a07-b6ae-28524c36246a" alt ="n-butane dihedral angle rotational energy profile" width="500" />
  <p> <b>Fig 2</b> Dihedral angle rotational energy profile of n-Butane</p>
</div>

From <b>Fig 2</b>, we can see that the gauche/anti-periplanar conformations represent energy minimas in the dihedral angle rotational energy profile of n-Butane, whereas the eclipsed conformations represent energy maximas. 

We can see that both the unstretched & stretched PEG molecules have the highest dihedral angle probability density at the gauche conformations (*ϕ* = 60° & *ϕ* = 180°) and the lowest dihedral angle probability density at the eclipsed conformations (*ϕ* = 0° & *ϕ* = 120°). This is because they are the lowest and highest energy conformations respectively.

The stretched PEG molecule has higher probability densities at *ϕ* = 60° & *ϕ* = 180° and lower probability densities at ϕ* = 0° & *ϕ* = 120°. This is consistent with assumption that the stretched PEG molecule has a more ordered/less conformationally free structure, reducing the range of observed dihedral angles.
