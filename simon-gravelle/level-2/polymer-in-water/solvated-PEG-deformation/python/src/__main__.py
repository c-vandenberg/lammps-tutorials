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
    # Unpack time and PEG length data
    time, peg_length = numpy.loadtxt(
        '../../data/raw/length-vs-time/ave_PEG_end_to_end_distance.dat'
    ).T

    # Convert time from femtoseconds to picoseconds
    time /= 1000

    # Convert PEG length data from Å (10e-10) to nm (10e-9)
    peg_length /= 10

    # Combine time and atom coordination number data into ndarray
    peg_length_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, peg_length))
    ]

    # Instantiate LineGraph object and create 'atom type 1 coordination number vs time' line graph
    solvated_peg_line_graph: LineGraph = LineGraph()
    solvated_peg_line_graph.single_line_graph(
        data_arrays=peg_length_vs_time_data_array,
        figure_size=(18, 10),
        line_colours='cyan',
        x_label=r'$t$ (ps)',
        y_label=r'$PEG\ Length\ (End to End)$ (nm)',
        y_lim=(0, 6),
        x_lim=(0, 60),
        graph_title=r'$\bf{PEG\ Length\ (End to End)\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of end-to-end PEG length as a function of time. The force is applied at '
                    r'~ $t$ = 31.5',
        figure_text_font_size=15,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
