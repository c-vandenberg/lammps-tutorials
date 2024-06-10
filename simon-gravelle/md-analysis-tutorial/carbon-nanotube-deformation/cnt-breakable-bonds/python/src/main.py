#!/usr/bin/env python3
import os
import numpy
from numpy import ndarray
from typing import List
from visualisation.breakable_cnt_bonds_plot import BreakableCNTBondsPlot
from MDAnalysis import AtomGroup, Universe


def main():
    base_dir: str = os.getcwd()

    # Instantiate MD Universe object with `../../data/raw/cnt_atomic.data` molecular topology data
    # & `../../data/raw/cnt_breakable_bonds_dump.lammpstrj` simulation trajectory coordinates
    md_universe: Universe = Universe(
        os.path.join(base_dir, "../../data/raw/topology/cnt_atomic.data"),
        os.path.join(base_dir, "../../data/raw/trajectory/cnt_breakable_bonds_dump.lammpstrj"),
        topology_format="data",
        format="lammpsdump",
        atom_style='id type x y z',
        guess_bonds=True,
        vdwradii={'1': 1.7}
    )

    # Instantiate carbon atoms (atom type 1) AtomGroup object
    cnt_atom_group: AtomGroup = md_universe.select_atoms('type 1')

    # Extract 'bond length vs timestep frame' and 'bond number vs timestep frame' data
    BreakableCNTBondsPlot.extract_mean_bond_lengths_bond_numbers(
        md_universe=md_universe,
        cnt_atom_group=cnt_atom_group,
        bond_length_vs_timestep_path=os.path.join(
            base_dir,
            '../../data/processed/bond-length-vs-time/bond_length_vs_timestep_frame.dat'
        ),
        bond_number_vs_timestep_path=os.path.join(
            base_dir,
            '../../data/processed/bond-number-vs-time/bond_number_vs_timestep_frame.dat'
        )
    )

    # Load 'bond length vs timestep frame' and 'bond number vs timestep frame' data
    bond_length_vs_timestep_frame: ndarray = numpy.loadtxt(
        os.path.join(
            base_dir,
            "../../data/processed/bond-length-vs-time/bond_length_vs_timestep_frame.dat"
        )
    )
    bond_number_vs_timestep_frame: ndarray = numpy.loadtxt(
        os.path.join(
            base_dir,
            "../../data/processed/bond-number-vs-time/bond_number_vs_timestep_frame.dat"
        )
    )

    # Define subplot configurations
    subplots_data_arrays: List[ndarray] = [bond_length_vs_timestep_frame, bond_number_vs_timestep_frame]

    # Create 'bond length vs timestep frame' and 'bond number vs timestep frame' subplots
    BreakableCNTBondsPlot.line_graph_subplots(
        data_arrays=subplots_data_arrays,
        subplot_titles=['CNT Bond Length vs Timestep Frame (a)', 'CNT Number of Bonds vs Timestep Frame (b)'],
        x_labels=['t (ps)', 't (ps)'],
        y_labels=['Bond Length (Å)', 'Number of Bonds'],
        x_lim=(0, 300),
        y_lims=[(1.35, 1.65), (500, 520)],
        graph_title=r'$\bf{CNT\ Average\ Bond\ Length\ &\ Bond\ Number\ vs\ Timestep\ Frame}$',
        figure_text=(
            r'$\bf{Fig\ 1}$ Evolution of carbon nanotube (CNT) average bond length (a) and bond number (b) as a '
            r'function of time.'
        )
    )

    # Extract bond length distribution data
    BreakableCNTBondsPlot.extract_bond_length_distributions(
        md_universe=md_universe,
        cnt_atom_group=cnt_atom_group,
        first_bond_length_distribution_path=(
            os.path.join(
                base_dir,
                '../../data/processed/bond-length-distribution/starting_bond_length_distribution.dat'
            )
        ),
        second_bond_length_distribution_path=(
            os.path.join(
                base_dir,
                '../../data/processed/bond-length-distribution/maximum_deformation_bond_length_distribution.dat'
            )
        ),
        bond_length_cutoff=1.8,
        number_of_bins=50,
        bond_length_range=(1.3, 1.65),
        first_frame_range=(0, 20),
        second_frame_range=(200, 220)
    )

    # Load bond length distribution data and plot on custom line graph
    starting_bond_length_distributions_data = numpy.loadtxt(
        os.path.join(
            base_dir,
            '../../data/processed/bond-length-distribution/starting_bond_length_distribution.dat'
        )
    ).T
    maximum_deformation_bond_length_distributions_data = (
        numpy.loadtxt(
            os.path.join(
                base_dir,
                '../../data/processed/bond-length-distribution/maximum_deformation_bond_length_distribution.dat'
            )
        ).T
    )

    # Create line graph to plot bond length probability densities
    bond_length_distributions_data = [
        starting_bond_length_distributions_data,
        maximum_deformation_bond_length_distributions_data
    ]

    BreakableCNTBondsPlot.single_line_graph(
        data_arrays=bond_length_distributions_data,
        figure_size=(10, 6),
        line_colours=['cyan', 'orange'],
        line_labels=['At Start (Frames 1 - 20)', 'During Maximum Deformation (Frames 200 - 220)'],
        x_label='Bond Length (Å)',
        y_label='Probability',
        y_lim=(0.00, 0.13),
        x_lim=(1.30, 1.65),
        graph_title=r'$\bf{Carbon\ Nanotube\ Bond\ Length\ Probability\ Densities}$',
        figure_text=(r'$\bf{Fig\ 2}$ Carbon nanotube (CNT) bond length distribution at start of simulation & at '
                     r'maximum deformation.'),
        figure_text_font_size=12,
        figure_text_x_coord=0.5,
        figure_text_y_coord=0.005,
        font_size=12,
        tick_label_size=10,
        line_width=1.5
    )


if __name__ == '__main__':
    main()
