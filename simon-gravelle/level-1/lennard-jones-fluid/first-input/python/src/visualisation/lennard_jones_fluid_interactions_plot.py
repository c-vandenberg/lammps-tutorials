import os
import sys
import numpy
from typing import List, Tuple
from numpy import ndarray

sys.path.append(
    os.getenv('LAMMPS_MD_ANALYSIS_BASE_DIRECTORY')
)

from src.modules.line_graph import LineGraph


class LennardJonesFluidInteractionsPlot(LineGraph):
    @staticmethod
    def calculate_lennard_jones_potentials(r: ndarray, lj_parameters: List[Tuple[float, float]]):
        """
        Compute Lennard-Jones potentials for a given range and list of (sigma, epsilon) pairs.

        Parameters
        ----------
        r : ndarray
            Array of inter-particle distances.
        lj_parameters : List[Tuple[float, float]]
            List of (sigma, epsilon) pairs.

        Returns
        -------
        List[ndarray]
            List of computed Lennard-Jones potentials for each (sigma, epsilon) pair.
        """
        lj_potentials: List[ndarray] = []
        for sigma, epsilon in lj_parameters:
            lj_potential = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
            lj_potentials.append(numpy.vstack((r, lj_potential)))
        return lj_potentials
