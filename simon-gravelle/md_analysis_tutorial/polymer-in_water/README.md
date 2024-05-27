# MDAnalysis Tutorials - Polymer in Water

In this tutorial, we will use the topology file and subsequent trajectory file generated in the water solvated PEG exercise.

For analysing and processing data via Python, this exercise utilises both a procedural approach (code within the Jupyter Lab notebook), and an object-oriented approach (code in `src` directory Python classes).

## Extract Temporal Evolution of Hydrogen Type 4 Atom (First Atom in PEG Group)

In this exercise, the position of the first atom of the peg group (i.e. the hydrogen type 4 atom) is extracted over all 300 frames, and its coordinates in each timestep frame stored in a list.

Matplotlib Pyplot is then used to visualise the x & y coordinates occupied by the hydrogen type 4 atom during the equilibration simulation:
<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/087ccf0d-8065-4f0e-9edb-0916d983b65b" alt="two_dimensional_H_atom_four_temporal_evolution" width="" />
</p>

Matplotlib Pyplot is also used to visualise the x, y & z coordinates occupied by the hydrogen type 4 atom during the equilibration simulation, along with the timestep frames:
<p align="center">
  <img src="https://github.com/c-vandenberg/lammps-tutorials/assets/60201356/bf54f3c7-7dfe-4e76-adaa-56baea2ff3ba" alt="three_dimensional_H_atom_four_temporal_evolution" width="" />
</p>
