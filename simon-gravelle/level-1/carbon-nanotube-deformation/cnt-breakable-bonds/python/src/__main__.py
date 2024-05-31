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
    # Extract first-input-log.lammps log file data & instantiate lammps_logfile.File object
    log_file: lammps_logfile.File = lammps_logfile.File('../../logs/cnt-breakable-bonds-log.lammps')
    total_energy_vs_time_array: List[ndarray] = []

    # Extract first run time from lammps_logfile.File object and convert to ps
    time_first_run: ndarray = log_file.get('Step', run_num=0) / 2000
    # Extract first run total energy from lammps_logfile.File object
    total_energy_first_run: ndarray = log_file.get('TotEng', run_num=0)

    total_energy_vs_time_array.append(numpy.vstack((time_first_run, total_energy_first_run)))

    # Extract first run time from lammps_logfile.File object and convert to ps
    time_second_run: ndarray = log_file.get('Step', run_num=1) / 2000
    # Extract first run total energy from lammps_logfile.File object
    total_energy_second_run: ndarray = log_file.get('TotEng', run_num=1)

    total_energy_vs_time_array.append(numpy.vstack((time_second_run, total_energy_second_run)))

    # Instantiate line graph object and create 'CNT system total energy vs time' line graph
    unbreakable_cnt_bonds_line_graph: LineGraph = LineGraph()
    unbreakable_cnt_bonds_line_graph.single_line_graph(
        data_arrays=total_energy_vs_time_array,
        figure_size=(18, 10),
        line_colours=['orange', 'cyan'],
        x_label=r'$t$ (ps)',
        y_label=r'$E_{tot}$ (Kcal/mol)',
        y_lim=(-5000, -3800),
        x_lim=(0, 150),
        graph_title=r'$\bf{Breakable\ Bonds\ Carbon\ Nanotube\ (CNT)\ System\ Total\ Energy\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of total energy of CNT system as a function of time.',
        figure_text_font_size=17.5,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
