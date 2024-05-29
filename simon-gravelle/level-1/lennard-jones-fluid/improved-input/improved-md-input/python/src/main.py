#!/usr/bin/env python3
import sys
import os
import numpy
from numpy import ndarray
from typing import List

sys.path.append(
    os.getenv('LAMMPS_MD_ANALYSIS_BASE_DIRECTORY')
)

from modules.line_graph import LineGraph


def main():
    time, atom_type_1_pop_inside_cylinder = numpy.loadtxt(
        '../../data/raw/atom-population/output_atom_type_1_population_inside_cylinder_vs_time.dat'
    ).T
    time -= time[0]
    time *= 0.005

    _, atom_type_2_pop_inside_cylinder = numpy.loadtxt(
        '../../data/raw/atom-population/output_atom_type_2_population_inside_cylinder_vs_time.dat'
    ).T

    atom_population_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, atom_type_1_pop_inside_cylinder)),
        numpy.vstack((time, atom_type_2_pop_inside_cylinder))
    ]

    improved_input_md_line_graph: LineGraph = LineGraph()

    improved_input_md_line_graph.single_line_graph(
        data_arrays=atom_population_vs_time_data_array,
        figure_size=(18, 10),
        line_labels=['Atom Type 1 Population Inside Cylinder', 'Atom Type 2 Population Inside Cylinder'],
        line_colours=['orange', 'cyan'],
        x_label=r'$t$',
        y_label=r'$N(inside)$',
        y_lim=(0, 225),
        x_lim=(0, 1500),
        graph_title='Atom Type 1 & Atom Type 2 Population Inside Cylinder vs Time',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of number of atoms within the $region_cylinder_in$ region as a, '
                    r'function of time',
        figure_text_font_size=15,
        font_size=15,
        label_size=20,
        line_width=3.5
    )

    _, atom_type_1_coordination_number = numpy.loadtxt(
        '../../data/raw/atom-coordination-number/output_average_atom_type_1_coordination_number.dat'
    ).T

    atom_coordination_number_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, atom_type_1_coordination_number))
    ]

    improved_input_md_line_graph.single_line_graph(
        data_arrays=atom_coordination_number_vs_time_data_array,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$',
        y_label=r'$Coordination\ Number$',
        y_lim=(0.00, 0.05),
        x_lim=(0, 1500),
        graph_title='Atom Type 1 Coordination Number vs Time',
        figure_text=r'$\bf{Fig\ 2}$ Evolution of atom type 1 coordination number within the $region_cylinder_in$ '
                    r'region as a function of time',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
