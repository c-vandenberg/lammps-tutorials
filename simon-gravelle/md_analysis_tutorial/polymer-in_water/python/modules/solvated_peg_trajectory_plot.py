#!/usr/bin/env python3

from typing import List, Dict, Union
import sys
from MDAnalysis import AtomGroup

sys.path.append(
    '/home/chris-vdb/Computational-Chemistry/lammps-tutorials/simon-gravelle/md_analysis_tutorial/python'
)

from src.modules.scatter_plot import ScatterPlot
from src.modules.md_universe import MDUniverse


class SolvatedPEGTrajectoryPlot(MDUniverse, ScatterPlot):
    def get_first_atom_temporal_evolution_from_type(self, molecule: AtomGroup, atom_type: str):
        try:
            first_atom: AtomGroup = molecule.select_atoms(atom_type)[0]
        except Exception as e:
            raise ValueError(e)

        position_vs_time: List[Dict[str, Union[int, float]]] = []

        for timestep in self._md_universe.trajectory:
            x, y, z = first_atom.position
            position_vs_time.append({'timestep': timestep.frame, 'x': x, 'y': y, 'z': z})

        return position_vs_time
