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
    # Instantiate line graph object
    pure_h2o_line_graph: LineGraph = LineGraph()

    # File contains columns of data, the .T transpose operation switches them to rows, making it easier to
    # unpack into separate variables
    # Extract H2O temperature data
    time, h2o_temperature = numpy.loadtxt(
        '../../data/raw/temperature-vs-time/ave_H2O_temperature.dat'
    ).T

    # Subtracts the first element from every element in the time array, ensuring the first element is 0/the time starts
    # at zero
    time -= time[0]

    # Convert time from femtoseconds to picoseconds
    time /= 1000

    # Extract H2O volume data
    _, h2o_volume = numpy.loadtxt(
        '../../data/raw/volume-vs-time/ave_H2O_volume.dat'
    ).T

    # Convert H2O volume from cubic Å to cubic nm
    h2o_volume /= 1000

    # Extract time and H2O density data
    _, h2o_density = numpy.loadtxt(
        '../../data/raw/density-vs-time/ave_H2O_density.dat'
    ).T

    # Convert density from molecule/Å-3 to molecule/nm-3
    h2o_density *= 1000

    # Combine time, and H2O temperature into ndarray
    h2o_temperature_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, h2o_temperature)),
    ]

    # Combine time, and H2O volume into ndarray
    h2o_volume_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, h2o_volume)),
    ]

    # Combine time, and H2O density into ndarray
    h2o_density_vs_time_data_array: List[ndarray] = [
        numpy.vstack((time, h2o_density)),
    ]

    # Create 'H2O temperature vs time' line graph
    pure_h2o_line_graph.single_line_graph(
        data_arrays=h2o_temperature_vs_time_data_array,
        figure_size=(18, 10),
        line_colours=['cyan'],
        x_label=r'$t$ (ps)',
        y_label=r'$Temperature$ (K)',
        y_lim=(25, 325),
        x_lim=(0, 20),
        graph_title=r'$\bf{Pure\ H_{2}O\ Temperature\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 1}$ Evolution of pure H$_{2}$O temperature as a function of time',
        figure_text_font_size=17.5,
        font_size=20,
        label_size=20,
        line_width=3.5
    )

    # Create 'H2O volume vs time' line graph
    pure_h2o_line_graph.single_line_graph(
        data_arrays=h2o_volume_vs_time_data_array,
        figure_size=(18, 10),
        line_colours=['cyan'],
        x_label=r'$t$ (ps)',
        y_label=r'$Volume$ (nm$^{3}$)',
        y_lim=(30, 90),
        x_lim=(0, 20),
        graph_title=r'$\bf{Pure\ H_{2}O\ Volume\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 2}$ Evolution of pure H$_{2}$O volume as a function of time',
        figure_text_font_size=17.5,
        font_size=20,
        label_size=20,
        line_width=3.5
    )

    # Create 'H2O density vs time' line graph
    pure_h2o_line_graph.single_line_graph(
        data_arrays=h2o_density_vs_time_data_array,
        figure_size=(18, 10),
        line_colours=['cyan'],
        x_label=r'$t$ (ps)',
        y_label=r'$p$ (nm$^{-3}$)',
        y_lim=(12, 36),
        x_lim=(0, 20),
        graph_title=r'$\bf{Pure\ H_{2}O\ Density\ vs\ Time}$',
        figure_text=r'$\bf{Fig\ 3}$ Evolution of pure H$_{2}$O volume as a function of time',
        figure_text_font_size=17.5,
        font_size=20,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
