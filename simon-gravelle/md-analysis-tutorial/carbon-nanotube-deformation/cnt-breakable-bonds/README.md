# MDAnalysis Tutorials - Carbon Nanotube (CNT) Breakable Bonds

In this tutorial, we will use the topology file and subsequent trajectory file generated in the deformation of a CNT with breakable bonds exercise.

For analysing and processing data via Python, this exercise utilises both a procedural approach (code within the Jupyter Lab notebook), and an object-oriented approach (code in `src` directory Python classes).

Before the main exercises there is some basic MDAnalysis manipulations in the Jupyter Lab notebook (counting atoms & timestep frames, accessing bonded atom indices, and extracting atomic positions)

## Evolution of CNT Average Bond Length & Bond Number as a Function of Time
<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/bd03ad37-523d-41d9-8dc5-cd9d3c060891" alt="cnt_average_bond_length_bond_number_vs_timestep" width="" />
</p>

In order to measure the evolution of the number of bonds over time, we loop over the simulation trajectory and manually extract the inter-atomic distance over time.

Within this simulation trajectory loop, there is a nested loop to iterate over the indices of the atoms that are detected as bonded & calculate the distance between the two atoms. The mean bond lengths and the total number of bonds will then be stored in two separate lists.

This data is then plotted on two scatter graph subplots via Matplotlib Pyplot. The subplots display 'bond length vs timestep frame' and 'bond number vs timestep frame'.

## Bond Length Distributions
<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/f061e4e7-9d28-4fb7-b29b-ca5f485940a6" alt="cnt_bond_length_distributions" width="" />
</p>

Using similar logic as above, this exercise extracts bond length distribution at the beginning of the simulation (frames 1 - 20), and during maximum deformation (frames 200 - 220).

This is achieved by carrying out a histogram calculation of all bond lengths in 50 bins ranging from 1.3 to 1.65 Å. These histogram counts are then normalized to obtain the probability density of each bond length.

As you would expect, prior to deformation the CNT bond length distributions show a sharp & high probability density between 1.41 to 1.43 Å, indicating that the majority of CNT bond lengths lie in this range during the simulation prior to deformation. Once maximum deformation is reached however, the bond length probabiltiy density is shows a much broader range with a higher probability density for bond lengths up to 1.65 Å. This reflects the stretching of the CNT bonds as a result of the deformation.
