from modules.solvated_peg_trajectory_plot import SolvatedPEGTrajectoryPlot
from MDAnalysis import AtomGroup
from constants.solvated_peg_constants import SolvatedPEGConstants


def main():
    solvated_peg_traj_plot: SolvatedPEGTrajectoryPlot = SolvatedPEGTrajectoryPlot()
    solvated_peg_traj_plot.md_universe(
        '../solvated_PEG.data',
        '../solvated_PEG_dump.lammpstrj',
        'data',
        'lammpsdump'
    )

    peg_molecule: AtomGroup = solvated_peg_traj_plot.get_atom_group(SolvatedPEGConstants.PEG_MOLECULE_ATOMS)

    peg_hydrogen_atom_4_traj: list = solvated_peg_traj_plot.get_first_atom_temporal_evolution_from_type(
        peg_molecule,
        'type 4'
    )

    timestep_frames: list = [position['timestep'] for position in peg_hydrogen_atom_4_traj]
    x_coordinates: list = [position['x'] for position in peg_hydrogen_atom_4_traj]
    y_coordinates: list = [position['y'] for position in peg_hydrogen_atom_4_traj]
    z_coordinates: list = [position['z'] for position in peg_hydrogen_atom_4_traj]

    solvated_peg_traj_plot.two_dimensional_scatter_plot(
        x_coordinates,
        y_coordinates,
        'x (Å)',
        'y (Å)',
        'X and Y Coordinates of Hydrogen Atom 4 During Equilibration'
    )

    solvated_peg_traj_plot.three_dimensional_scatter_plot(
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
