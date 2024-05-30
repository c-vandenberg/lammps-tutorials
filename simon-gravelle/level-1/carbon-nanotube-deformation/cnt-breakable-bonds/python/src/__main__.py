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
    # File contains columns of data, the .T transpose operation switches them to rows, making it easier to
    # unpack into separate variables
    time, cnt_length = numpy.loadtxt(
        '../data/raw/length-vs-time/output_cnt_length.dat'
    ).T

    # Convert time from femtoseconds to picoseconds
    time /= 1000

    # Combine time, atom type 1 & atom type 2 data into two separate ndarrays
    atom_population_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, cnt_length)),
    ]

    # Instantiate line graph object and create 'atom population vs time' line graph
    unbreakable_cnt_bonds_line_graph: LineGraph = LineGraph()
    unbreakable_cnt_bonds_line_graph.single_line_graph(
        data_arrays=atom_population_vs_time_data_array,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$',
        y_label=r'$Length(Å)$',
        y_lim=(55, 67),
        x_lim=(0, 16),
        graph_title=r'$\bf{Carbon\ Nanotube\ (CNT)\ Length\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of CNT length as a function of time',
        figure_text_font_size=15,
        font_size=15,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()