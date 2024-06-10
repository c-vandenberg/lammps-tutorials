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
    base_dir: str = os.getcwd()
    cnt_unbreakable_bonds_path: str = os.path.join(base_dir, '../../cnt-unbreakable-bonds/data/raw/stress-strain/')
    cnt_stress_strain_data_array: List[ndarray] = []

    # File contains columns of data, the .T transpose operation switches them to rows, making it easier to
    # unpack into separate variables. Ignore time with `_` placeholder
    # Extract CNT unbreakable bonds strain
    _, cnt_unbreakable_bonds_strain = numpy.loadtxt(
        os.path.join(cnt_unbreakable_bonds_path, 'output_cnt_strain.dat')
    ).T

    # Extract CNT unbreakable bonds stress
    _, cnt_unbreakable_bonds_stress = numpy.loadtxt(
        os.path.join(cnt_unbreakable_bonds_path, 'output_cnt_stress.dat')
    ).T

    # Convert stress from Pa to GPa
    cnt_unbreakable_bonds_stress /= 1e9

    # Combine time, and CNT length into ndarray
    cnt_stress_strain_data_array.append(
        numpy.vstack((cnt_unbreakable_bonds_strain, cnt_unbreakable_bonds_stress))
    )

    cnt_breakable_bonds_path: str = os.path.join(base_dir, '../../cnt-breakable-bonds/data/raw/stress-strain/')

    # Extract CNT breakable bonds strain
    _, cnt_breakable_bonds_strain = numpy.loadtxt(
        os.path.join(cnt_breakable_bonds_path, 'output_cnt_strain.dat')
    ).T

    # Extract CNT breakable bonds stress
    _, cnt_breakable_bonds_stress = numpy.loadtxt(
        os.path.join(cnt_breakable_bonds_path, 'output_cnt_stress.dat')
    ).T

    # Convert strain from Pa to GPa
    cnt_breakable_bonds_stress /= 1e9

    # Combine CNT breakable bonds strain and stress into ndarray
    cnt_stress_strain_data_array.append(
        numpy.vstack((cnt_breakable_bonds_strain, cnt_breakable_bonds_stress))
    )

    # Create 'CNT stress v strain' line graph
    LineGraph.single_line_graph(
        data_arrays=cnt_stress_strain_data_array,
        figure_size=(18, 10),
        line_labels=['CNT - Unbreakable Bonds', 'CNT - Breakable Bonds'],
        line_colours=['cyan', 'orange'],
        x_label=r'$Strain$ (%)',
        y_label=r'$Stress$ (GPa)',
        x_lim=(-1, 35),
        y_lim=(-50, 650),
        graph_title=r'$\bf{Unbreakable\ &\ Breakable Bonds\ Carbon\ Nanotube\ (CNT)\ Stress\ vs\ Strain}$',
        figure_text=r'$\bf{Fig\ 1}$ Stress-strain curves for unbreakable & breakable CNTs',
        figure_text_font_size=17.5,
        figure_text_x_coord=0.5,
        figure_text_y_coord=0.005,
        font_size=20,
        tick_label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
