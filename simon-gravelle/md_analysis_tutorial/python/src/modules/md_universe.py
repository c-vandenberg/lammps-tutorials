import MDAnalysis as MDAnalysis
from MDAnalysis import AtomGroup

class MDUniverse:
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