import os
import sys
from typing import List, Dict, Union
from MDAnalysis import AtomGroup, Universe

sys.path.append(
    os.getenv('LAMMPS_MD_ANALYSIS_BASE_DIRECTORY')
)

from modules.scatter_plot import ScatterPlot


class SolvatedPEGTrajectoryPlot(ScatterPlot):
    @staticmethod
    def get_first_atom_temporal_evolution_from_type(md_universe: Universe, molecule: AtomGroup, atom_type: str):
        """
        Get the temporal evolution of the first atom of a specified type in a molecule.

        Parameters
        ----------
        md_universe : Universe
            The MDAnalysis Universe object containing the simulation data.
        molecule : AtomGroup
            The AtomGroup representing the molecule from which to select the atom.
        atom_type : str
            The type of atom to select (e.g., 'type 1').

        Returns
        -------
        List[Dict[str, Union[int, float]]]
            A list of dictionaries containing the timestep and the x, y, z coordinates of the atom at each timestep.

        Raises
        ------
        ValueError
            If the atom type is not found in the molecule or another MDAnalysis error occurs.
        """
        try:
            first_atom: AtomGroup = molecule.select_atoms(atom_type)[0]
        except Exception as e:
            raise ValueError(e)

        position_vs_time: List[Dict[str, Union[int, float]]] = []

        for timestep in md_universe.trajectory:
            x, y, z = first_atom.position
            position_vs_time.append({'timestep': timestep.frame, 'x': x, 'y': y, 'z': z})

        return position_vs_time
