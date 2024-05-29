#!/usr/bin/env python3
import lammps_logfile
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
        figure_text_font_size=20,
        font_size=20,
        label_size=20,
        line_width=4

    )

    log_file: lammps_logfile.File = lammps_logfile.File('../log.lammps')
    timestep: float = 0.005
    pe_vs_time: List[ndarray] = []
    ke_vs_time: List[ndarray] = []

    energy_min_time: ndarray = log_file.get("Step", run_num=0) * timestep
    energy_min_potential_energy: ndarray = log_file.get("PotEng", run_num=0)
    energy_min_kinetic_energy: ndarray = log_file.get("KinEng", run_num=0)

    pe_vs_time.append(numpy.vstack((energy_min_time, energy_min_potential_energy)))
    ke_vs_time.append(numpy.vstack((energy_min_time, energy_min_kinetic_energy)))

    molecular_dynamics_time: ndarray = log_file.get("Step", run_num=1) * timestep
    molecular_dynamics_potential_energy: ndarray = log_file.get("PotEng", run_num=1)
    molecular_dynamics_kinetic_energy: ndarray = log_file.get("KinEng", run_num=1)

    pe_vs_time.append(numpy.vstack((molecular_dynamics_time, molecular_dynamics_potential_energy)))
    ke_vs_time.append(numpy.vstack((molecular_dynamics_time, molecular_dynamics_kinetic_energy)))

    lennard_jones_fluid_interactions_plot.single_line_graph(
        data_arrays=pe_vs_time,
        figure_size=(18, 10),
        line_labels=['Energy Minimization', 'Molecular Dynamics'],
        line_colours=['orange', 'cyan'],
        x_label=r'$t$',
        y_label=r'$p_{e}$',
        y_lim=(-1.0, 0.0),
        x_lim=(0, 20),
        graph_title='Lennard-Jones Fluid Potential Energy vs Time (During Energy Minimization & Molecular Dynamics '
                    'Simulation)',
        figure_text=r'$\bf{Fig\ 2}$ Potential energy ($p_{e}$) rapidly decreases during energy minimization (orange),'
                    r'increases once the molecular dynamics simulation starts (blue) before plateauing at ~ -0.25',
        figure_text_font_size=12,
        font_size=15,
        label_size=20,
        line_width=3.5
    )

    lennard_jones_fluid_interactions_plot.single_line_graph(
        data_arrays=ke_vs_time,
        figure_size=(18, 10),
        line_labels=['Energy Minimization', 'Molecular Dynamics'],
        line_colours=['orange', 'cyan'],
        x_label=r'$t$',
        y_label=r'$p_{e}$',
        y_lim=(-0.05, 2.0),
        x_lim=(0, 20),
        graph_title='Lennard-Jones Fluid Kinetic Energy vs Time (During Energy Minimization & Molecular Dynamics '
                    'Simulation)',
        figure_text=r'$\bf{Fig\ 3}$ Kinetic energy ($k_{e}$) is zero during energy minimization (orange) due to no '
                    r'atomic or molecular motion. It increases once the molecular dynamics simulation starts (blue) '
                    r'before plateauing at ~ -1.5',
        figure_text_font_size=12,
        font_size=15,
        label_size=20,
        line_width=3.5
    )


if __name__ == '__main__':
    main()
