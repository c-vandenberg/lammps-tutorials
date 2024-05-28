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
        md_universe,
        cnt_atom_group,
        '../data/processed/bond_length_vs_timestep_frame.dat',
        '../data/processed/bond_number_vs_timestep_frame.dat'
    )

    # Load 'bond length vs timestep frame' and 'bond number vs timestep frame' data
    bond_length_vs_timestep_frame = numpy.loadtxt("../data/processed/bond_length_vs_timestep_frame.dat")
    bond_number_vs_timestep_frame = numpy.loadtxt("../data/processed/bond_number_vs_timestep_frame.dat")

    # Define subplot configurations
    subplots_data_arrays: List[ndarray] = [bond_length_vs_timestep_frame, bond_number_vs_timestep_frame]

    # Create 'bond length vs timestep frame' and 'bond number vs timestep frame' subplots
    breakable_cnt_bonds_plot.line_graph_subplots(
        subplots_data_arrays,
        ['CNT Bond Length vs Timestep Frame (a)', 'CNT Number of Bonds vs Timestep Frame (b)'],
        ['t (ps)', 't (ps)'],
        ['Bond Length (Å)', 'Number of Bonds'],
        [(1.35, 1.65), (500, 520)],
        (0, 300),
        'CNT Average Bond Length & Bond Number vs Timestep Frame',
        (r'$\bf{Fig\ 1}$ Evolution of carbon nanotube (CNT) average bond length (a) and bond number (b) as a '
         r'function of time.')
    )

    # Extract bond length distribution data
    breakable_cnt_bonds_plot.extract_bond_length_distributions(
        md_universe,
        cnt_atom_group,
        '../data/processed/starting_bond_length_distribution.dat',
        '../data/processed/maximum_deformation_bond_length_distribution.dat',
        1.8,
        50,
        (1.3, 1.65),
        (0, 20),
        (200, 220)
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
        bond_length_distributions_data,
        (10, 6),
        ['At Start (Frames 1 - 20)', 'During Maximum Deformation (Frames 200 - 220)'],
        ['cyan', 'orange'],
        'Bond Length (Å)',
        'Probability',
        (0.00, 0.13),
        (1.30, 1.65),
        'Carbon Nanotube Bond Length Probability Densities',
        (r'$\bf{Fig\ 2}$ Carbon nanotube (CNT) bond length distribution at start of simulation & at '
         r'maximum deformation.'),
        12,
        10,
        1.5
    )


if __name__ == '__main__':
    main()
