#!/usr/bin/env python3

import os
import MDAnalysis as MDAnalysis
import numpy
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D
from ..constants import md_analysis_constants


class SolvatedPEGTrajectoryPlot:
    def __init__(self):
        self._md_universe = None

    def md_universe(self, topology_file_path: str, trajectory_file_path: str):
        topology_file_format: str = os.path.splitext(topology_file_path)[1]
        trajectory_file_format: str = os.path.splitext(trajectory_file_path)[1]

        supported_topology_formats = MDAnalysis.topology.topology.guess_format.__globals__['_READERS']
        supported_trajectory_formats = MDAnalysis.coordinates.core.reader.__globals__['_READERS']

        if topology_file_format not in supported_topology_formats:
            raise ValueError(f"Unsupported topology file format: {topology_file_format}")

        if trajectory_file_format not in supported_trajectory_formats:
            raise ValueError(f"Unsupported trajectory file format: {trajectory_file_format}")

        self._md_universe: MDAnalysis.Universe = MDAnalysis.Universe(
            topology_file_path,
            trajectory_file_path,
            topology_format=topology_file_format,
            format=trajectory_file_format
        )



# Instantiate MD Universe object with `./solvated_PEG.data` molecular topology data
# & `./solvated_PEG_dump.lammpstrj` simulation trajectory coordinates
md_universe: MDAnalysis.Universe = MDAnalysis.Universe(
    "./solvated_PEG.data",
    "./solvated_PEG_dump.lammpstrj",
    topology_format="data",
    format="lammpsdump"
)

peg_molecule = md_universe.select_atoms("type 1 2 3 4 5 6 7")
h2o_molecule = md_universe.select_atoms("type 8 9")

print("Atoms in PEG molecule:", peg_molecule.atoms.n_atoms)
print("Atoms in H2O molecule:", h2o_molecule.atoms.n_atoms)

for atom in peg_molecule[:6]:
    atom_id: int = atom.id
    atom_type: str = atom.type
    atomic_mass: float = atom.mass
    atomic_charge: float = numpy.round(atom.charge, 2)
    print("Atom ID:", atom_id, "|",
          "Atom Type:", atom_type, "|",
          "Atomic Mass:", atomic_mass, "g/mol | Atomic Charge:", atomic_charge, "e")

hydrogen_atom_4: md_analysis.AtomGroup = peg_molecule[0]
position_vs_time: list = []
scatter_plot_size: list = []

for timestep in md_universe.trajectory:
    x, y, z = hydrogen_atom_4.position
    position_vs_time.append([timestep.frame, x, y, z])  # ts.frame is the ID of the timestep frame

# Extract x and y coordinates
timestep_frames: list = [position[0] for position in position_vs_time]
x_coordinates: list = [position[1] for position in position_vs_time]
y_coordinates: list = [position[2] for position in position_vs_time]
z_coordinates: list = [position[3] for position in position_vs_time]

# Create graph
two_d_figure, two_d_axes = pyplot.subplots(figsize=(10, 6))

# Stylise & plot the x and y coordinates using Matplotlib
two_d_figure.patch.set_facecolor('black')
two_d_axes.set_facecolor('black')
two_d_axes.scatter(x_coordinates, y_coordinates, c='cyan', marker='o', alpha=0.6)
two_d_axes.set_ylim(12, two_d_axes.get_ylim()[1])

pyplot.title('X and Y Coordinates of Hydrogen Atom 4 During Equilibration', color='white')
pyplot.xlabel('x (Å)', color='white')
pyplot.ylabel('y (Å)', color='white')
two_d_axes.tick_params(colors='white', which='both')

for spine in two_d_axes.spines.values():
    spine.set_edgecolor('white')

pyplot.show()

three_d_figure: pyplot.Figure = pyplot.figure(figsize=(10, 6))
three_d_axes: pyplot.Axes = three_d_figure.add_subplot(111, projection='3d')

three_dimension_scatter = three_d_axes.scatter(
    x_coordinates,
    y_coordinates,
    z_coordinates,
    c=timestep_frames,
    cmap='viridis_r',
    marker='o'
)

# Add color bar to show frame mapping
cb = pyplot.colorbar(three_dimension_scatter, ax=three_d_axes, shrink=0.5)
cb.set_label('Timestep Frame')

# Set titles and labels
three_d_axes.set_title('3D Coordinates of Hydrogen Atom 4 During Equilibration')
three_d_axes.set_xlabel('x (Å)')
three_d_axes.set_ylabel('y (Å)')
three_d_axes.set_zlabel('z (Å)')

# Show the plot
pyplot.show()
