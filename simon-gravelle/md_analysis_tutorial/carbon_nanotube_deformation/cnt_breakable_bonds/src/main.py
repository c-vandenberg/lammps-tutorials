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



    # Extract bond length distribution data
    breakable_cnt_bonds_plot.extract_bond_length_distributions(
        md_universe,
        cnt_atom_group,
        '../data/processed/starting_bond_length_distribution.dat',
        '../data/processed/maximum_deformation_bond_length_distribution.dat'
    )

    starting_bond_length_distributions_data = numpy.loadtxt('../data/processed/starting_bond_length_distribution.dat').T
    maximum_deformation_bond_length_distributions_data = (
        numpy.loadtxt('../data/processed/maximum_deformation_bond_length_distribution.dat').T)

    bond_length_distributions_data = [
        starting_bond_length_distributions_data,
        maximum_deformation_bond_length_distributions_data
    ]

    bond_length_plot_title: List = [
        'CNT Bond Length Distributions'
    ]
    bond_length_plot_x_labels: List = ['Bond Length (Å)']
    bond_length_plot_y_labels: List = ['Probability']
    bond_length_plot_y_lims: Union[Tuple, List[Tuple]] = (0.00, 1.17)
    bond_length_plot_x_lim: Union[Tuple, List[Tuple]] = (1.25, 1.65)
    bond_length_plot_figure_title: str = (r'$\bf{Fig\ 2}$ Bond length distribution carbon nanotube (CNT) at start of '
                                          r'simulation & at maximum deformation.')

    breakable_cnt_bonds_plot.plot_bond_length_distributions(starting_bond_length_distributions_data,
                                                            maximum_deformation_bond_length_distributions_data)


if __name__ == '__main__':
    main()
