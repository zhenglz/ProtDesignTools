
from typing import List, Dict, Optional, Union
import re

class Sequence:
    """Class to handle protein sequences and mutations."""
    
    def __init__(self, sequence: str, name: str = "WT"):
        self.wt_sequence = sequence
        self.current_sequence = sequence
        self.name = name
        self.mutations: List[str] = []

    def apply_mutations(self, mutations_str: str, separator: str = ","):
        """
        Apply mutations to the sequence.
        Format: 'A12G, S30C' (1-based indexing)
        """
        if not mutations_str or mutations_str.upper() in ["WT", "NATIVE"]:
            self.current_sequence = self.wt_sequence
            self.mutations = []
            return self.current_sequence

        mut_list = [m.strip() for m in mutations_str.split(separator)]
        seq_chars = list(self.wt_sequence)
        
        applied_muts = []
        for mut in mut_list:
            # Parse mutation string like 'A12G'
            match = re.match(r"([A-Z])(\d+)([A-Z])", mut)
            if not match:
                raise ValueError(f"Invalid mutation format: {mut}")
            
            wt_aa, pos, mut_aa = match.groups()
            pos = int(pos)
            
            if pos < 1 or pos > len(seq_chars):
                raise ValueError(f"Position {pos} out of range for sequence length {len(seq_chars)}")
            
            if seq_chars[pos-1] != wt_aa:
                print(f"Warning: Original AA at position {pos} is {seq_chars[pos-1]}, but mutation specifies {wt_aa}")
            
            seq_chars[pos-1] = mut_aa
            applied_muts.append(mut)
            
        self.current_sequence = "".join(seq_chars)
        self.mutations = applied_muts
        return self.current_sequence

    def get_sequence(self) -> str:
        return self.current_sequence

    def __len__(self):
        return len(self.current_sequence)

    def __str__(self):
        return self.current_sequence

    @classmethod
    def from_fasta(cls, fasta_path: str):
        # Basic fasta parser
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            name = lines[0].strip().lstrip('>')
            seq = "".join([line.strip() for line in lines[1:]])
        return cls(seq, name=name)
