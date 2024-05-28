#!/usr/bin/env python3

from visualisation.solvated_peg_trajectory_plot import SolvatedPEGTrajectoryPlot
from MDAnalysis import AtomGroup, Universe
from constants.solvated_peg_constants import SolvatedPEGConstants


def main():
    # Instantiate MD Universe object with `../data/raw/solvated_PEG.data` molecular topology data
    # & `../data/raw/solvated_PEG_dump.lammpstrj` simulation trajectory coordinates
    md_universe: Universe = Universe(
        '../data/raw/solvated_PEG.data',
        '../data/raw/solvated_PEG_dump.lammpstrj',
        topology_format='data',
        format='lammpsdump'
    )

    # Instantiate dedicated class for plotting polymer in water data
    solvated_peg_trajectory_plot: SolvatedPEGTrajectoryPlot = SolvatedPEGTrajectoryPlot()

    # Group PEG atoms together into AtomGroup object
    peg_molecule: AtomGroup = md_universe.select_atoms(SolvatedPEGConstants.PEG_MOLECULE_ATOMS)

    # Extract temporal evolution data of first type 4 hydrogen atom in PEG molecule
    peg_hydrogen_atom_4_traj: list = solvated_peg_trajectory_plot.get_first_atom_temporal_evolution_from_type(
        md_universe,
        peg_molecule,
        'type 4'
    )

    # Extract temporal evolution data into separate list variables for timestep frames and x, y, z coordinates
    timestep_frames: list = [position['timestep'] for position in peg_hydrogen_atom_4_traj]
    x_coordinates: list = [position['x'] for position in peg_hydrogen_atom_4_traj]
    y_coordinates: list = [position['y'] for position in peg_hydrogen_atom_4_traj]
    z_coordinates: list = [position['z'] for position in peg_hydrogen_atom_4_traj]

    # Plot x and y coordinates occupied by first type 4 hydrogen atom during simulation on 2-D scatter plot
    solvated_peg_trajectory_plot.two_dimensional_scatter_plot(
        x_coordinates,
        y_coordinates,
        'x (Å)',
        'y (Å)',
        'X and Y Coordinates of Hydrogen Atom 4 During Equilibration'
    )

    # Plot x, y, z coordinates on 3-D scatter plot, with colour bar representing timestep frames
    solvated_peg_trajectory_plot.three_dimensional_scatter_plot(
        x_coordinates,
        y_coordinates,
        z_coordinates,
        timestep_frames,
        'x (Å)',
        'y (Å)',
        'z (Å)',
        'Timestep Frame',
        '3D Coordinates of Hydrogen Atom 4 During Equilibration'
    )


if __name__ == '__main__':
    main()
