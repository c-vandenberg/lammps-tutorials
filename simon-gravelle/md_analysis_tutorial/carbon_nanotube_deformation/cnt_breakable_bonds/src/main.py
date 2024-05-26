import numpy
from numpy import ndarray
from visualisation.breakable_cnt_bonds_plot import BreakableCNTBondsPlot
from MDAnalysis import AtomGroup, Universe
from constants.breakable_cnt_bonds_constants import BreakableCNTBondsConstants


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

    breakable_cnt_bonds_plot: BreakableCNTBondsPlot = BreakableCNTBondsPlot()

    # Instantiate carbon atoms (atom type 1) AtomGroup object
    cnt = md_universe.select_atoms('type 1')

    breakable_cnt_bonds_plot.extract_bond_lengths_bond_numbers(
        md_universe,
        cnt,
        '../data/processed/bond_length_vs_timestep_frame.dat',
        '../data/processed/bond_number_vs_timestep_frame.dat'
    )

    # Load data
    bond_length_vs_timestep_frame = numpy.loadtxt("../data/processed/bond_length_vs_timestep_frame.dat")
    bond_number_vs_timestep_frame = numpy.loadtxt("../data/processed/bond_number_vs_timestep_frame.dat")

    # Define plot configurations
    data_arrays = [bond_length_vs_timestep_frame, bond_number_vs_timestep_frame]
    subplot_titles = [
        'CNT Bond Length vs Timestep Frame (a)',
        'CNT Number of Bonds vs Timestep Frame (b)'
    ]
    x_labels = ['t (ps)', 't (ps)']
    y_labels = ['Bond Length (Å)', 'Number of Bonds']
    y_lims = [(1.35, 1.65), (500, 520)]
    x_lim = (0, 300)
    figure_title = (r'$\bf{Fig\ 1}$ Evolution of carbon nanotube (CNT) average bond length (a) '
                    r'and bond number (b) as a function of time.')

    # Create dynamic subplots
    breakable_cnt_bonds_plot.line_graph_subplots(data_arrays, subplot_titles, x_labels, y_labels, y_lims, x_lim, figure_title)


if __name__ == '__main__':
    main()
