#!/usr/bin/env python3

from typing import List, Dict, Union
import sys

sys.path.append(
    '/home/chris-vdb/Computational-Chemistry/lammps-tutorials/simon-gravelle/md_analysis_tutorial/python/src'
)

import MDAnalysis as MDAnalysis
from MDAnalysis import AtomGroup
import numpy
import matplotlib.pyplot as pyplot
from modules.scatter_plot import ScatterPlot


class SolvatedPEGTrajectoryPlot(ScatterPlot):
    def __init__(self):
        self._md_universe = None

    def md_universe(
            self,
            topology_file_path: str,
            trajectory_file_path: str,
            topology_file_format: str,
            trajectory_file_format: str
    ):
        try:
            self._md_universe: MDAnalysis.Universe = MDAnalysis.Universe(
                topology_file_path,
                trajectory_file_path,
                topology_format=topology_file_format,
                format=trajectory_file_format
            )
        except Exception as e:
            raise ValueError(e)

    def get_atom_group(self, molecule_atoms: str) -> AtomGroup:
        try:
            return self._md_universe.select_atoms(molecule_atoms)
        except Exception as e:
            raise ValueError(e)

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
