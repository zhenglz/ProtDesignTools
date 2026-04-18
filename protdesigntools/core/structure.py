
import os
from typing import Dict, List, Optional, Any
from protdesigntools.core.sequence import Sequence
import numpy as np

class Structure(Sequence):
    """
    Class to handle protein structures. 
    Inherits from Sequence to manage the sequence aspect of the primary chain 
    or a combined sequence.
    """
    
    def __init__(self, pdb_path: str, name: Optional[str] = None):
        self.pdb_path = pdb_path
        self.name = name or os.path.basename(pdb_path).split('.')[0]
        self.chains: Dict[str, str] = {} # chain_id -> sequence
        self.coords: Dict[str, np.ndarray] = {} # chain_id -> atomic coordinates (CA atoms)
        self._parse_pdb()
        
        # Initialize base Sequence with the first chain's sequence
        first_chain_id = sorted(self.chains.keys())[0] if self.chains else ""
        super().__init__(self.chains.get(first_chain_id, ""), name=self.name)

    def _parse_pdb(self):
        """Simple PDB parser to extract sequences and CA coordinates per chain."""
        if not os.path.exists(self.pdb_path):
            return

        current_chain = None
        chain_seqs = {}
        chain_coords = {}
        
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        with open(self.pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    chain_id = line[21].strip() or "A"
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    if chain_id not in chain_seqs:
                        chain_seqs[chain_id] = []
                        chain_coords[chain_id] = []
                    
                    chain_seqs[chain_id].append(aa_map.get(res_name, 'X'))
                    chain_coords[chain_id].append([x, y, z])
        
        for cid in chain_seqs:
            self.chains[cid] = "".join(chain_seqs[cid])
            self.coords[cid] = np.array(chain_coords[cid])

    def get_chain_sequence(self, chain_id: str) -> str:
        return self.chains.get(chain_id, "")

    def get_all_chains(self) -> List[str]:
        return list(self.chains.keys())

    def save_pdb(self, output_path: str):
        """Save the structure to a PDB file (placeholder for now)."""
        # In a real tool, we might use Biopython or similar to save/modify PDBs
        import shutil
        shutil.copy(self.pdb_path, output_path)

    def get_rmsd(self, other: 'Structure', chain_id: str = "A") -> float:
        """Calculate RMSD between two structures for a given chain."""
        if chain_id not in self.coords or chain_id not in other.coords:
            return -1.0
        
        c1 = self.coords[chain_id]
        c2 = other.coords[chain_id]
        
        if len(c1) != len(c2):
            # Simple truncation for demonstration; real tools use alignment
            min_len = min(len(c1), len(c2))
            c1 = c1[:min_len]
            c2 = c2[:min_len]
            
        diff = c1 - c2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
