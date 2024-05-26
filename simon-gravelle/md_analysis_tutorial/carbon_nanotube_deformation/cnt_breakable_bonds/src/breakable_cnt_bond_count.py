#!/usr/bin/env python3

import MDAnalysis as md_analysis
import numpy
import matplotlib.pyplot as pyplot
from numpy import ndarray

# Instantiate MD Universe object with `./cnt_atomic.data` molecular topology data
# & `./cnt_breakable_bonds_dump.lammpstrj` simulation trajectory coordinates
md_universe: md_analysis.Universe = md_analysis.Universe(
    "../../data/raw/cnt_atomic.data",
    "../../data/raw/cnt_breakable_bonds_dump.lammpstrj",
    topology_format="data",
    format="lammpsdump",
    atom_style='id type x y z',
    guess_bonds=True,
    vdwradii={'1': 1.7}
)

# Instantiate carbon atoms (atom type 1) AtomGroup object, and perform basic MD analysis
cnt = md_universe.select_atoms("type 1")

# Extract Mean Bond Lengths & Number of Bonds
bond_length_vs_timestep_frame: list = []
bond_number_vs_timestep_frame: list = []

for timestep in md_universe.trajectory:
    frame = timestep.frame
    current_timestep_bond_lengths: list = []
    for index_1, index_2 in cnt.atoms.bonds.indices:
        position_1 = md_universe.atoms.positions[md_universe.atoms.indices == index_1]
        position_2 = md_universe.atoms.positions[md_universe.atoms.indices == index_2]

        bond_length: float = numpy.sqrt(numpy.sum((position_1 - position_2) ** 2))

        if bond_length < 1.8:
            current_timestep_bond_lengths.append(bond_length)

    mean_bond_length = numpy.mean(current_timestep_bond_lengths)
    number_of_bonds: float = round(len(current_timestep_bond_lengths) / 2)  # Divide by two to avoid counting twice

    bond_length_vs_timestep_frame.append([frame, mean_bond_length])
    bond_number_vs_timestep_frame.append([frame, number_of_bonds])

# Save 'bond length vs timestep frame' and 'bond number vs timestep frame' data to file
numpy.savetxt("../../data/processed/bond_length_vs_timestep_frame.dat", bond_length_vs_timestep_frame)
numpy.savetxt("../../data/processed/bond_number_vs_timestep_frame.dat", bond_number_vs_timestep_frame)

# Load data
bond_length_vs_timestep_frame: ndarray = numpy.loadtxt("../../data/processed/bond_length_vs_timestep_frame.dat")
bond_number_vs_timestep_frame: ndarray = numpy.loadtxt("../../data/processed/bond_number_vs_timestep_frame.dat")

# Create two line graphs/subplots
line_graphs, (bond_length_vs_timestep_frame_axes, bond_number_vs_timestep_frame_axes) = pyplot.subplots(
    2,
    1,
    figsize=(10, 6)
)
# Set background colour to black and grid coloir to white for both subplots
line_graphs.patch.set_facecolor('black')
bond_length_vs_timestep_frame_axes.set_facecolor('black')
bond_number_vs_timestep_frame_axes.set_facecolor('black')

# Set colour of ticks and spines on botha axes to white
for axes in [bond_length_vs_timestep_frame_axes, bond_number_vs_timestep_frame_axes]:
    axes.tick_params(colors='white', which='both')
    for spine in axes.spines.values():
        spine.set_edgecolor('white')


# Extract bond length vs timestep frame data and plot line graph subplot
bond_length_vs_timestep_frame_axes.plot(bond_length_vs_timestep_frame[:, 0], bond_length_vs_timestep_frame[:, 1],
                                        color='cyan')
bond_length_vs_timestep_frame_axes.set_title('CNT Bond Length vs Timestep Frame (a)', color='white')
bond_length_vs_timestep_frame_axes.set_xlabel('t (ps)', color='white')
bond_length_vs_timestep_frame_axes.set_ylabel('Bond Length (Å)', color='white')
bond_length_vs_timestep_frame_axes.set_ylim([1.35, 1.65])

# Extract number of bonds vs timestep frame data and plot line graph subplot
bond_number_vs_timestep_frame_axes.plot(bond_number_vs_timestep_frame[:, 0], bond_number_vs_timestep_frame[:, 1],
                                        color='cyan')
bond_number_vs_timestep_frame_axes.set_title('CNT Number of Bonds vs Timestep Frame (b)', color='white')
bond_number_vs_timestep_frame_axes.set_xlabel('t (ps)', color='white')
bond_number_vs_timestep_frame_axes.set_ylabel('Number of Bonds', color='white')
bond_number_vs_timestep_frame_axes.set_ylim([500, 520])

# Set x-axis limits for both plots (same for both)
bond_length_vs_timestep_frame_axes.set_xlim([0, 300])
bond_number_vs_timestep_frame_axes.set_xlim([0, 300])

# Adjust the spacing between subplots
pyplot.subplots_adjust(hspace=0.5)

line_graphs.text(
    0.5,
    0.0005,
    r'$\bf{Fig\ 1}$ Evolution of carbon nanotube (CNT) average bond length (a) '
    r'and bond number (b) as a function of time.',
    ha='center',
    va='center',
    color='white',
    fontsize=12
)

# Adjust layout and show line graphs
pyplot.tight_layout()
pyplot.show()
