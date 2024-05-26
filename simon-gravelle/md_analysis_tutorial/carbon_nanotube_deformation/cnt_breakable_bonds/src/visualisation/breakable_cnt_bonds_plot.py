#!/usr/bin/env python3

import sys
import numpy
from MDAnalysis import AtomGroup, Universe
from constants.breakable_cnt_bonds_constants import BreakableCNTBondsConstants

sys.path.append(
    BreakableCNTBondsConstants.BASE_DIRECTORY
)

from src.modules.line_graph import LineGraph


class BreakableCNTBondsPlot(LineGraph):
    @staticmethod
    def extract_bond_lengths_bond_numbers(
            md_universe: Universe,
            cnt_atom_group: AtomGroup,
            bond_length_vs_timestep_path: str,
            bond_number_vs_timestep_path: str
    ):
        bond_length_vs_timestep_frame: list = []
        bond_number_vs_timestep_frame: list = []

        for timestep in md_universe.trajectory:
            frame = timestep.frame
            current_timestep_bond_lengths: list = []
            for index_1, index_2 in cnt_atom_group.atoms.bonds.indices:
                position_1 = md_universe.atoms.positions[md_universe.atoms.indices == index_1]
                position_2 = md_universe.atoms.positions[md_universe.atoms.indices == index_2]

                bond_length: float = numpy.sqrt(numpy.sum((position_1 - position_2) ** 2))

                if bond_length < 1.8:
                    current_timestep_bond_lengths.append(bond_length)

            mean_bond_length = numpy.mean(current_timestep_bond_lengths)
            number_of_bonds: float = round(
                len(current_timestep_bond_lengths) / 2)  # Divide by two to avoid counting twice

            bond_length_vs_timestep_frame.append([frame, mean_bond_length])
            bond_number_vs_timestep_frame.append([frame, number_of_bonds])

        # Save 'bond length vs timestep frame' and 'bond number vs timestep frame' data to file
        numpy.savetxt(bond_length_vs_timestep_path, bond_length_vs_timestep_frame)
        numpy.savetxt(bond_number_vs_timestep_path, bond_number_vs_timestep_frame)
