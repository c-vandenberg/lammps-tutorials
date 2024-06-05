import os
import sys
import numpy
from typing import Tuple
from numpy import ndarray
from MDAnalysis import AtomGroup, Universe

sys.path.append(
    os.getenv('LAMMPS_MD_ANALYSIS_BASE_DIRECTORY')
)

from modules.line_graph import LineGraph


class BreakableCNTBondsPlot(LineGraph):
    @staticmethod
    def extract_mean_bond_lengths_bond_numbers(md_universe: Universe, cnt_atom_group: AtomGroup,
                                               bond_length_vs_timestep_path: str, bond_number_vs_timestep_path: str):
        """
        Extract mean bond lengths and bond numbers from a given molecular dynamics universe and atom group,
        and save the results to specified file paths.

        Parameters
        ----------
        md_universe : Universe
            The MDAnalysis Universe object containing the simulation data.
        cnt_atom_group : AtomGroup
            The AtomGroup representing the carbon nanotube (CNT) from which to extract bond lengths and numbers.
        bond_length_vs_timestep_path : str
            File path to save the bond length versus timestep frame data.
        bond_number_vs_timestep_path : str
            File path to save the bond number versus timestep frame data.

        Returns
        -------
        None
        """
        bond_length_vs_timestep_frame: list = []
        bond_number_vs_timestep_frame: list = []

        for timestep in md_universe.trajectory:
            frame: int = timestep.frame
            current_timestep_bond_lengths: list = []

            # `cnt_atom_group.atoms.bonds.indices` is a list of tuples and the expression `index_1, index_2`
            # unpacks each tuple (index_1 and index_2 represent the indices of the two atoms in each bond)
            for index_1, index_2 in cnt_atom_group.atoms.bonds.indices:
                position_1: ndarray = md_universe.atoms.positions[md_universe.atoms.indices == index_1]
                position_2: ndarray = md_universe.atoms.positions[md_universe.atoms.indices == index_2]

                bond_length: float = numpy.sqrt(numpy.sum((position_1 - position_2) ** 2))

                if bond_length < 1.8:
                    current_timestep_bond_lengths.append(bond_length)

            mean_bond_length: float = numpy.mean(current_timestep_bond_lengths)
            number_of_bonds: int = round(
                len(current_timestep_bond_lengths) / 2)  # Divide by two to avoid counting twice

            bond_length_vs_timestep_frame.append([frame, mean_bond_length])
            bond_number_vs_timestep_frame.append([frame, number_of_bonds])

        # Save 'bond length vs timestep frame' and 'bond number vs timestep frame' data to file
        numpy.savetxt(bond_length_vs_timestep_path, bond_length_vs_timestep_frame)
        numpy.savetxt(bond_number_vs_timestep_path, bond_number_vs_timestep_frame)

    @staticmethod
    def extract_bond_length_distributions(md_universe: Universe, cnt_atom_group: AtomGroup,
                                          first_bond_length_distribution_path: str,
                                          second_bond_length_distribution_path: str,
                                          bond_length_cutoff: float, number_of_bins: int,
                                          bond_length_range: Tuple[float, float], first_frame_range: Tuple[int, int],
                                          second_frame_range: Tuple[int, int]):
        """
        Extract bond length distributions from a given molecular dynamics universe and atom group,
        and save the results to specified file paths.

        Parameters
        ----------
        md_universe : Universe
            The MDAnalysis Universe object containing the simulation data.
        cnt_atom_group : AtomGroup
            The AtomGroup representing the carbon nanotube (CNT) from which to extract bond lengths.
        first_bond_length_distribution_path : str
            File path to save the first bond length distribution.
        second_bond_length_distribution_path : str
            File path to save the second bond length distribution.
        bond_length_cutoff : float
            Maximum bond length to consider.
        number_of_bins : int
            Number of bins to use in the histogram.
        bond_length_range : Tuple[float, float]
            Range of bond lengths to consider in the histogram.
        first_frame_range : Tuple[int, int]
            Start and end frame indices for the first bond length distribution.
        second_frame_range : Tuple[int, int]
            Start and end frame indices for the second bond length distribution.

        Returns
        -------
        None
        """
        bond_length_distributions = []
        for timestep in md_universe.trajectory:
            frame: int = timestep.frame
            current_timestep_bond_lengths: list = []
            for index_1, index_2 in cnt_atom_group.atoms.bonds.indices:
                position_1: ndarray = md_universe.atoms.positions[md_universe.atoms.indices == index_1]
                position_2: ndarray = md_universe.atoms.positions[md_universe.atoms.indices == index_2]

                bond_length: float = numpy.sqrt(numpy.sum((position_1 - position_2) ** 2))

                if bond_length < bond_length_cutoff:
                    current_timestep_bond_lengths.append(bond_length)

            if frame > 0:  # Ignore first frame
                # Histogram calculation of all bond lengths and 50 bins ranging from 1.3 to 1.65
                # Variable `bond_length_histo` contains counts of bond lengths in each bin
                bond_length_histo, bin_edges = numpy.histogram(current_timestep_bond_lengths, bins=number_of_bins,
                                                               range=bond_length_range)

                # Convert bin edges to bin centers by averaging each pair of adjacent bin edges
                bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

                # Normalize histogram counts by calculating bin widths & dividing histogram counts by total number of
                # bond lengths * the bin width
                bin_width = bin_edges[1] = bin_edges[0]
                bond_length_histo = bond_length_histo / (numpy.sum(bond_length_histo) * bin_width)

                # Store bin center and normalized histogram counts by stacking vertically using `numpy.vstack`, and
                # appending to `bond_length_distributions` list
                bond_length_distributions.append(numpy.vstack([bin_centers, bond_length_histo]))

        # Slice `bond_length_distributions` list to extract bond length distributions for first frame range
        # Calculate mean of these distributions along the vertical axis
        first_bond_length_distribution = numpy.mean(
            bond_length_distributions[first_frame_range[0]:first_frame_range[1]], axis=0
        )

        # Slice `bond_length_distributions` list to extract bond length distributions between second frame range
        second_deformation_bond_length_distribution = numpy.mean(
            bond_length_distributions[second_frame_range[0]:second_frame_range[1]], axis=0
        )

        numpy.savetxt(first_bond_length_distribution_path, first_bond_length_distribution.T)
        numpy.savetxt(second_bond_length_distribution_path, second_deformation_bond_length_distribution.T)
