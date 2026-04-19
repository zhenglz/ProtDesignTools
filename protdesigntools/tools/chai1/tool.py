
import os
import argparse
import logging
import tempfile
import glob
import numpy as np
import shutil
from typing import Dict, Any, List, Optional
from protdesigntools.core.base_tool import BaseTool
from protdesigntools.core.structure import Structure
from time import sleep

# We might need BioPython for CIF to PDB conversion
try:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.PDBIO import PDBIO
    from Bio.PDB.PDBParser import PDBParser
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

logger = logging.getLogger(__name__)

class Chai1(BaseTool):
    """
    Chai-1 for structure prediction.
    Inputs: Sequence string(s), FASTA file.
    Outputs: Predicted PDB structure, pLDDT, iPTM (if complex).
    """

    @classmethod
    def get_cli_parser(cls) -> argparse.ArgumentParser:
        parser = super().get_cli_parser()
        parser.add_argument("--pdb_path", type=str, help="Input PDB file to extract sequences from")
        parser.add_argument("--fasta_path", type=str, help="Input FASTA file containing sequences")
        parser.add_argument("--sequence", type=str, help="Single protein sequence string")
        parser.add_argument("--ligand", type=str, help="Optional ligand SMILES for complex prediction")
        parser.add_argument("--target_dna_seq", type=str, help="Target DNA sequence to replace DNA chains in PDB")
        return parser

    def _prepare_fasta(self, input_params: Dict[str, Any], temp_dir: str) -> str:
        """Prepare the input FASTA file for Chai-1."""
        fasta_out = os.path.join(temp_dir, "input.fasta")
        
        if input_params.get("fasta_path") and os.path.exists(input_params["fasta_path"]):
            shutil.copy(input_params["fasta_path"], fasta_out)
            return fasta_out
            
        pdb_path = input_params.get("pdb_path")
        target_dna_seq = input_params.get("target_dna_seq")
        
        if pdb_path and os.path.exists(pdb_path):
            # Parse chains from PDB and generate multi-chain fasta
            struct = Structure(pdb_path)
            
            with open(fasta_out, "w") as f:
                # Naive heuristic: if sequence contains only ATCG, consider it DNA
                for chain_id, seq in struct.chains.items():
                    is_dna = set(seq).issubset({'A', 'T', 'C', 'G', 'N', 'U'})
                    
                    if is_dna and target_dna_seq:
                        # Simple replacement logic
                        f.write(f">DNA|{chain_id}\n{target_dna_seq}\n")
                        # For a full implementation, we'd handle complementary strands properly
                    elif is_dna:
                        f.write(f">DNA|{chain_id}\n{seq}\n")
                    else:
                        f.write(f">protein|{chain_id}\n{seq}\n")
                        
            return fasta_out

        seq_str = input_params.get("sequence")
        if not seq_str:
            raise ValueError("Chai-1 requires either fasta_path, sequence, or pdb_path.")
            
        with open(fasta_out, "w") as f:
            f.write(f">protein|protein\n{seq_str}\n")
            
            # Chai-1 supports ligands in fasta via special headers
            ligand = input_params.get("ligand")
            if ligand:
                f.write(f">ligand|ligand\n{ligand}\n")
                
        return fasta_out

    def _convert_cif_to_pdb(self, cif_path: str, pdb_path: str):
        """Convert CIF to PDB using BioPython if available."""
        if not HAS_BIOPYTHON:
            logger.warning("BioPython not installed. Cannot convert CIF to PDB. Copying CIF instead.")
            shutil.copy(cif_path, pdb_path)
            return

        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("chai_model", cif_path)
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_path)
        except Exception as e:
            logger.error(f"Failed to convert CIF to PDB: {e}")
            shutil.copy(cif_path, pdb_path)

    def _extract_plddt(self, struct_path: str) -> float:
        """Extract average pLDDT from the B-factor column of CA atoms."""
        plddt_sum = 0.0
        count = 0
        
        if not os.path.exists(struct_path):
            return 0.0
            
        with open(struct_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        # B-factor is usually columns 60-66
                        b_factor = float(line[60:66].strip())
                        plddt_sum += b_factor
                        count += 1
                    except ValueError:
                        pass
                        
        return plddt_sum / count if count > 0 else 0.0

    def _extract_iptm(self, npz_path: str) -> float:
        """Extract iPTM from Chai-1 scores .npz file."""
        if not os.path.exists(npz_path):
            return 0.0
            
        try:
            data = np.load(npz_path, allow_pickle=True)
            if 'iptm' in data:
                # iptm might be a scalar or a 0-d array
                return float(data['iptm'])
            elif 'ptm' in data: # Fallback to ptm if iptm not present
                return float(data['ptm'])
        except Exception as e:
            logger.error(f"Error reading {npz_path}: {e}")
            
        return 0.0

    def _parse_chai_output(self, temp_dir: str, output_dir: str) -> Dict[str, Any]:
        """Find and parse the best model from Chai-1 outputs."""
        # Chai-1 usually outputs to a subfolder or directly in the given dir
        # Look for .cif files
        cif_files = glob.glob(os.path.join(temp_dir, "*.cif"))
        if not cif_files:
            # Check subdirectories
            cif_files = glob.glob(os.path.join(temp_dir, "*", "*.cif"))
            
        if not cif_files:
            raise RuntimeError(f"No .cif output found in {temp_dir}")
            
        # Prioritize model_idx_0 if it exists
        best_cif = cif_files[0]
        for f in cif_files:
            if "model_idx_0" in f:
                best_cif = f
                break
                
        base_name = os.path.splitext(os.path.basename(best_cif))[0]
        
        # Convert to PDB
        final_pdb = os.path.join(output_dir, f"{base_name}.pdb")
        self._convert_cif_to_pdb(best_cif, final_pdb)
        
        # Calculate pLDDT from the converted PDB (or CIF)
        plddt = self._extract_plddt(final_pdb if os.path.exists(final_pdb) else best_cif)
        
        # Look for corresponding npz file for iPTM
        npz_files = glob.glob(os.path.join(os.path.dirname(best_cif), "*.npz"))
        best_npz = None
        for f in npz_files:
            if "scores" in f and "model_idx_0" in f:
                best_npz = f
                break
        if not best_npz and npz_files:
            best_npz = npz_files[0]
            
        iptm = self._extract_iptm(best_npz) if best_npz else 0.0
        
        return {
            "predicted_pdb": final_pdb,
            "plddt": plddt,
            "iptm": iptm,
            "combined_score": (plddt / 100.0 + iptm) / 2.0 if iptm > 0 else plddt / 100.0
        }

    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        logger.info("Running Chai-1 Structure Prediction")
        
        if "output_dir" in input_params:
            self.output_dir = input_params["output_dir"]
        os.makedirs(self.output_dir, exist_ok=True)
        
        script_path = self.config.get("script_path", "/sugon_store/pub_data/tools/chai-lab/run.sh")
        
        results = {
            "tool": "Chai-1",
            "status": "failed"
        }
        
        # Keep tempdir active outside the block if running asynchronously or just don't clean it up immediately.
        # Alternatively, create a standard directory.
        temp_dir = os.path.join(self.output_dir, "chai_input_fasta")
        os.makedirs(temp_dir, exist_ok=True)
        
        fasta_path = self._prepare_fasta(input_params, temp_dir)
        # sleep 2s
        sleep(2)
        
        # The Chai-1 CLI usually takes the fasta and output directory
        args = [
            fasta_path,
            self.output_dir
        ]
        
        cmd = self.build_command(script_path, args)
        
        job_id = self.execute(cmd, job_name="chai1_pred")
        
        # Parse outputs directly from self.output_dir
        parsed_data = self._parse_chai_output(self.output_dir, self.output_dir)
        
        results.update(parsed_data)
        results["job_id"] = job_id
        results["output_dir"] = self.output_dir
        results["status"] = "success"

        return results

if __name__ == "__main__":
    Chai1.cli()
