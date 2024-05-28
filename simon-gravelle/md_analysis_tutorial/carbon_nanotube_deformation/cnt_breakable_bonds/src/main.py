#!/usr/bin/env python3

import numpy
from numpy import ndarray
from typing import List, Union, Tuple
from visualisation.breakable_cnt_bonds_plot import BreakableCNTBondsPlot
from MDAnalysis import AtomGroup, Universe


def main():
    # Instantiate MD Universe object with `../../data/raw/cnt_atomic.data` molecular topology data
    # & `../../data/raw/cnt_breakable_bonds_dump.lammpstrj` simulation trajectory coordinates
    md_universe: Universe = Universe(
        "../data/raw/cnt_atomic.data",
        "../data/raw/cnt_breakable_bonds_dump.lammpstrj",
        topology_format="data",
        format="lammpsdump",
        atom_style='id type x y z',
        guess_bonds=True,
        vdwradii={'1': 1.7}
    )

    # Instantiate dedicated class for plotting CNT breakable bonds data
    breakable_cnt_bonds_plot: BreakableCNTBondsPlot = BreakableCNTBondsPlot()

    # Instantiate carbon atoms (atom type 1) AtomGroup object
    cnt_atom_group: AtomGroup = md_universe.select_atoms('type 1')

    # Extract 'bond length vs timestep frame' and 'bond number vs timestep frame' data
    breakable_cnt_bonds_plot.extract_mean_bond_lengths_bond_numbers(
        md_universe=md_universe,
        cnt_atom_group=cnt_atom_group,
        bond_length_vs_timestep_path='../data/processed/bond_length_vs_timestep_frame.dat',
        bond_number_vs_timestep_path='../data/processed/bond_number_vs_timestep_frame.dat'
    )

    # Load 'bond length vs timestep frame' and 'bond number vs timestep frame' data
    bond_length_vs_timestep_frame = numpy.loadtxt("../data/processed/bond_length_vs_timestep_frame.dat")
    bond_number_vs_timestep_frame = numpy.loadtxt("../data/processed/bond_number_vs_timestep_frame.dat")

    # Define subplot configurations
    subplots_data_arrays: List[ndarray] = [bond_length_vs_timestep_frame, bond_number_vs_timestep_frame]

    # Create 'bond length vs timestep frame' and 'bond number vs timestep frame' subplots
    breakable_cnt_bonds_plot.line_graph_subplots(
        data_arrays=subplots_data_arrays,
        subplot_titles=['CNT Bond Length vs Timestep Frame (a)', 'CNT Number of Bonds vs Timestep Frame (b)'],
        x_labels=['t (ps)', 't (ps)'],
        y_labels=['Bond Length (Å)', 'Number of Bonds'],
        y_lims=[(1.35, 1.65), (500, 520)],
        x_lims=(0, 300),
        graph_title='CNT Average Bond Length & Bond Number vs Timestep Frame',
        figure_text=(
            r'$\bf{Fig\ 1}$ Evolution of carbon nanotube (CNT) average bond length (a) and bond number (b) as a '
            r'function of time.'
        )
    )

    # Extract bond length distribution data
    breakable_cnt_bonds_plot.extract_bond_length_distributions(
        md_universe=md_universe,
        cnt_atom_group=cnt_atom_group,
        first_bond_length_distribution_path='../data/processed/starting_bond_length_distribution.dat',
        second_bond_length_distribution_path='../data/processed/maximum_deformation_bond_length_distribution.dat',
        bond_length_cutoff=1.8,
        number_of_bins=50,
        bond_length_range=(1.3, 1.65),
        first_frame_range=(0, 20),
        second_frame_range=(200, 220)
    )

    # Load bond length distribution data and plot on custom line graph
    starting_bond_length_distributions_data = numpy.loadtxt('../data/processed/starting_bond_length_distribution.dat').T
    maximum_deformation_bond_length_distributions_data = (
        numpy.loadtxt('../data/processed/maximum_deformation_bond_length_distribution.dat').T)

    # Define line graph configurations
    bond_length_distributions_data = [
        starting_bond_length_distributions_data,
        maximum_deformation_bond_length_distributions_data
    ]

    breakable_cnt_bonds_plot.single_line_graph(
        data_arrays=bond_length_distributions_data,
        figure_size=(10, 6),
        line_labels=['At Start (Frames 1 - 20)', 'During Maximum Deformation (Frames 200 - 220)'],
        line_colours=['cyan', 'orange'],
        x_label='Bond Length (Å)',
        y_label='Probability',
        y_lim=(0.00, 0.13),
        x_lim=(1.30, 1.65),
        graph_title='Carbon Nanotube Bond Length Probability Densities',
        figure_text=(r'$\bf{Fig\ 2}$ Carbon nanotube (CNT) bond length distribution at start of simulation & at '
                     r'maximum deformation.'),
        font_size=12,
        label_size=10,
        line_width=1.5
    )


if __name__ == '__main__':
    main()
