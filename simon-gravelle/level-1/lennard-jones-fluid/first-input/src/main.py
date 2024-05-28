#!/usr/bin/env python3
import numpy
from numpy import ndarray
from typing import List
from visualisation.lennard_jones_fluid_interactions_plot import LennardJonesFluidInteractionsPlot


def main():
    # Define intermolecular range values and Lennard-Jones parameters
    r: ndarray = numpy.arange(0.5, 10, 0.001)
    lj_parameters: List = [
        (1.0, 1.0),  # sigma_11, epsilon_11
        (3.0, 0.5),  # sigma_22, epsilon_22
        (numpy.sqrt(1.0 * 3.0), numpy.sqrt(1.0 * 0.5))  # sigma_12, epsilon_12
    ]

    # Calculate Lennard-Jones potentials
    lennard_jones_fluid_interactions_plot: LennardJonesFluidInteractionsPlot = LennardJonesFluidInteractionsPlot()
    lj_potentials: List[ndarray] = lennard_jones_fluid_interactions_plot.calculate_lennard_jones_potentials(
        r,
        lj_parameters
    )

    lennard_jones_fluid_interactions_plot.single_line_graph(
        data_arrays=lj_potentials,
        figure_size=(18, 6),
        line_labels=[r'$E_{12}$', r'$E_{22}$', r'$E_{11}$'],
        line_colours=['gray', 'cyan', 'orange'],
        x_label=r'$r$',
        y_label=r'$E(r)$',
        y_lim=(-1.055, 1.055),
        x_lim=(0, 6.1),
        graph_title='Lennard-Jones Potential vs Intermolecular Distance',
        figure_text=r'$\bf{Fig\ 1}$ Lennard-Jones potential $E_{ij}(r)$ as a function of the intermolecular distance',
        font_size=20,
        label_size=20,
        line_width=4
    )


if __name__ == '__main__':
    main()