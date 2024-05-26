#!/usr/bin/env python3

import sys
from typing import List, Dict, Union
from MDAnalysis import AtomGroup, Universe
from constants.solvated_peg_constants import SolvatedPEGConstants

sys.path.append(
    SolvatedPEGConstants.BASE_DIRECTORY
)

from src.modules.scatter_plot import ScatterPlot


class SolvatedPEGTrajectoryPlot(ScatterPlot):
    def get_first_atom_temporal_evolution_from_type(self, md_universe: Universe, molecule: AtomGroup, atom_type: str):
        try:
            first_atom: AtomGroup = molecule.select_atoms(atom_type)[0]
        except Exception as e:
            raise ValueError(e)

        position_vs_time: List[Dict[str, Union[int, float]]] = []

        for timestep in md_universe.trajectory:
            x, y, z = first_atom.position
            position_vs_time.append({'timestep': timestep.frame, 'x': x, 'y': y, 'z': z})

        return position_vs_time
