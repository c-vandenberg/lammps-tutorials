# MDAnalysis Tutorials - Carbon Nanotube (CNT) Breakable Bonds

In this tutorial, we will use the topology file and subsequent trajectory file generated in the deformation of a CNT with breakable bonds exercise.

For analysing and processing data via Python, this exercise utilises both a procedural approach (code within the Jupyter Lab notebook), and an object-oriented approach (code in `src` directory Python classes).

Before the main exercises there is some basic MDAnalysis manipulations in the Jupyter Lab notebook (counting atoms & timestep frames, accessing bonded atom indices, and extracting atomic positions)

## Evolution of CNT Average Bond Length & Bond Number as a Function of Time

In order to measure the evolution of the number of bonds over time, we loop over the simulation trajectory and manually extract the inter-atomic distance over time.

Within this simulation trajectory loop, there is a nested loop to iterate over the indices of the atoms that are detected as bonded & calculate the distance between the two atoms. The mean bond lengths and the total number of bonds will then be stored in two separate lists.

This data is then plotted on two scatter graph subplots via Matplotlib Pyplot. The subplots display 'bond length vs timestep frame' and 'bond number vs timestep frame'.

## Bond Length Distributions

Using similar logic as above, this exercise extracts bond length distribution at the beginning of the simulation (frames 1 - 20), and during maximum deformation (frames 200 - 220).

This is achieved by carrying out a histogram calculation of all bond lengths in 50 bins ranging from 1.3 to 1.65 Å. These histogram counts are then normalized to obtain the probability density of each bond length.
