from visualisation.solvated_peg_trajectory_plot import SolvatedPEGTrajectoryPlot
from MDAnalysis import AtomGroup, Universe
from constants.solvated_peg_constants import SolvatedPEGConstants


def main():
    md_universe: Universe = Universe(
        '../data/raw/solvated_PEG.data',
        '../data/raw/solvated_PEG_dump.lammpstrj',
        topology_format='data',
        format='lammpsdump'
    )

    solvated_peg_trajectory_plot: SolvatedPEGTrajectoryPlot = SolvatedPEGTrajectoryPlot()

    peg_molecule: AtomGroup = md_universe.select_atoms(SolvatedPEGConstants.PEG_MOLECULE_ATOMS)

    peg_hydrogen_atom_4_traj: list = solvated_peg_trajectory_plot.get_first_atom_temporal_evolution_from_type(
        md_universe,
        peg_molecule,
        'type 4'
    )

    timestep_frames: list = [position['timestep'] for position in peg_hydrogen_atom_4_traj]
    x_coordinates: list = [position['x'] for position in peg_hydrogen_atom_4_traj]
    y_coordinates: list = [position['y'] for position in peg_hydrogen_atom_4_traj]
    z_coordinates: list = [position['z'] for position in peg_hydrogen_atom_4_traj]

    solvated_peg_trajectory_plot.two_dimensional_scatter_plot(
        x_coordinates,
        y_coordinates,
        'x (Å)',
        'y (Å)',
        'X and Y Coordinates of Hydrogen Atom 4 During Equilibration'
    )

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
