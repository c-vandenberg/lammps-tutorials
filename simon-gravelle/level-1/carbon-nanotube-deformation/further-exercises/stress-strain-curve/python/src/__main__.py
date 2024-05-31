#!/usr/bin/env python3
import sys
import os
import lammps_logfile
import numpy
from numpy import ndarray
from typing import List

sys.path.append(
    os.getenv('LAMMPS_MD_ANALYSIS_BASE_DIRECTORY')
)

from modules.line_graph import LineGraph


def main():
    cnt_unbreakable_bonds_path: str = '../../cnt-unbreakable-bonds/data/raw/stress-strain/'
    cnt_stress_strain_data_array: List[ndarray] = []

    # File contains columns of data, the .T transpose operation switches them to rows, making it easier to
    # unpack into separate variables. Ignore time with `_` placeholder
    # Extract CNT unbreakable bonds strain
    _, cnt_unbreakable_bonds_strain = numpy.loadtxt(
        cnt_unbreakable_bonds_path + 'output_cnt_strain.dat'
    ).T

    # Extract CNT unbreakable bonds stress
    _, cnt_unbreakable_bonds_stress = numpy.loadtxt(
        cnt_unbreakable_bonds_path + 'output_cnt_stress.dat'
    ).T

    # Convert stress from Pa to GPa
    cnt_unbreakable_bonds_stress /= 1e9

    # Combine time, and CNT length into ndarray
    cnt_stress_strain_data_array.append(
        numpy.vstack((cnt_unbreakable_bonds_strain, cnt_unbreakable_bonds_stress))
    )

    cnt_breakable_bonds_path: str = '../../cnt-breakable-bonds/data/raw/stress-strain/'

    # Extract CNT breakable bonds strain
    _, cnt_breakable_bonds_strain = numpy.loadtxt(
        cnt_breakable_bonds_path + 'output_cnt_strain.dat'
    ).T

    # Extract CNT breakable bonds stress
    _, cnt_breakable_bonds_stress = numpy.loadtxt(
        cnt_breakable_bonds_path + 'output_cnt_stress.dat'
    ).T

    # Convert strain from Pa to GPa
    cnt_breakable_bonds_stress /= 1e9

    # Combine time, and CNT length into ndarray
    cnt_stress_strain_data_array.append(
        numpy.vstack((cnt_breakable_bonds_strain, cnt_breakable_bonds_stress))
    )

    # Instantiate line graph object and create 'CNT stress v strain' line graph
    unbreakable_cnt_bonds_line_graph: LineGraph = LineGraph()
    unbreakable_cnt_bonds_line_graph.single_line_graph(
        data_arrays=cnt_stress_strain_data_array,
        figure_size=(18, 10),
        line_labels=['CNT - Unbreakable Bonds', 'CNT - Breakable Bonds'],
        line_colours=['cyan', 'orange'],
        x_label=r'$Strain$ (%)',
        y_label=r'$Stress$ (GPa)',
        y_lim=(-50, 650),
        x_lim=(-1, 35),
        graph_title=r'$\bf{Unbreakable\ &\ Breakable Bonds\ Carbon\ Nanotube\ (CNT)\ Stress\ vs\ Strain}$',
        figure_text=r'$\bf{Fig\ 1}$ Stress-strain curves for unbreakable & breakable CNTs',
        figure_text_font_size=17.5,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
