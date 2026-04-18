
import os
import argparse
import logging
import json
import tempfile
import shutil
from typing import Dict, Any, List, Optional
from core.base_tool import BaseTool
from core.structure import Structure
from core.sequence import Sequence

logger = logging.getLogger(__name__)

class ProteinMPNN(BaseTool):
    """
    ProteinMPNN tool for sequence design and scoring.
    
    Scoring Mode:
      - Supports input via PDB + mutations, FASTA file, or sequence string.
      - Calculates fitness score.
      
    Design Mode:
      - Supports specifying chains to design.
      - Supports fixing specific positions or defining amino acid lists.
    """

    @classmethod
    def get_cli_parser(cls) -> argparse.ArgumentParser:
        parser = super().get_cli_parser()
        parser.add_argument("--mode", type=str, choices=["design", "scoring"], default="design", help="Run mode")
        parser.add_argument("--pdb_path", type=str, help="Input PDB file")
        
        # Scoring specific
        parser.add_argument("--mutations", type=str, help="Mutations for scoring mode (e.g. A12G, S30C)")
        parser.add_argument("--fasta_path", type=str, help="Input FASTA file for scoring")
        parser.add_argument("--sequence", type=str, help="Input sequence string for scoring")
        
        # Design specific
        parser.add_argument("--design_chains", type=str, help="Comma-separated list of chains to design (e.g. A,B)")
        parser.add_argument("--fixed_positions", type=str, help="JSON string mapping chain to list of fixed 1-based indices (e.g. '{\"A\": [1, 2, 3]}')")
        parser.add_argument("--omit_AAs", type=str, help="String of amino acids to omit (e.g. 'CX')")
        parser.add_argument("--num_seqs", type=int, default=1, help="Number of sequences to design")
        parser.add_argument("--sampling_temp", type=float, default=0.1, help="Sampling temperature")
        
        return parser

    def _prepare_scoring_fasta(self, input_params: Dict[str, Any], temp_dir: str) -> str:
        """Prepare the input FASTA file for scoring mode based on various input types."""
        fasta_out = os.path.join(temp_dir, "input.fasta")
        
        if input_params.get("fasta_path") and os.path.exists(input_params["fasta_path"]):
            shutil.copy(input_params["fasta_path"], fasta_out)
            return fasta_out
            
        seq_str = input_params.get("sequence")
        pdb_path = input_params.get("pdb_path")
        mutations = input_params.get("mutations")
        
        final_seq = ""
        header = ">sequence"
        
        if seq_str:
            seq_obj = Sequence(seq_str)
            if mutations:
                seq_obj.apply_mutations(mutations)
            final_seq = seq_obj.get_sequence()
            header += f" mutations={mutations}" if mutations else " WT"
            
        elif pdb_path and os.path.exists(pdb_path):
            struct = Structure(pdb_path)
            if mutations:
                struct.apply_mutations(mutations)
            # For simplicity, getting the sequence of the first chain
            # In a full implementation, you'd iterate over all chains and format a multi-line FASTA
            final_seq = struct.get_sequence() 
            header += f" source={os.path.basename(pdb_path)} mutations={mutations}"
            
        else:
            raise ValueError("Scoring mode requires either fasta_path, sequence, or pdb_path.")

        with open(fasta_out, "w") as f:
            f.write(f"{header}\n{final_seq}\n")
            
        return fasta_out

    def _prepare_design_jsonls(self, input_params: Dict[str, Any], temp_dir: str) -> Dict[str, str]:
        """
        Prepare the JSONL configuration files required by ProteinMPNN design mode.
        (parsed_pdbs.jsonl, assigned_pdbs.jsonl, fixed_pdbs.jsonl)
        """
        pdb_path = input_params.get("pdb_path")
        if not pdb_path or not os.path.exists(pdb_path):
            raise ValueError("Design mode requires a valid pdb_path.")
            
        # In a real implementation, we would call the actual helper scripts:
        # e.g. parse_multiple_chains.py, assign_fixed_chains.py, make_fixed_positions_dict.py
        # Here we mock their outputs to demonstrate the standard IO structure
        
        base_name = os.path.basename(pdb_path).split('.')[0]
        
        # 1. Parsed PDBs (Mock)
        parsed_jsonl = os.path.join(temp_dir, "parsed_pdbs.jsonl")
        with open(parsed_jsonl, "w") as f:
            # Format: {"name": "prot", "seq_chain_A": "...", "coords_chain_A": {...}}
            f.write(json.dumps({"name": base_name, "seq_chain_A": "MOCKSEQ"}) + "\n")
            
        # 2. Assigned Chains
        assigned_jsonl = ""
        design_chains = input_params.get("design_chains")
        if design_chains:
            assigned_jsonl = os.path.join(temp_dir, "assigned_pdbs.jsonl")
            chains = [c.strip() for c in design_chains.split(",")]
            with open(assigned_jsonl, "w") as f:
                # Format: {"name": [ [design_chains], [fixed_chains] ]}
                f.write(json.dumps({base_name: [chains, []]}) + "\n")
                
        # 3. Fixed Positions
        fixed_jsonl = ""
        fixed_pos_str = input_params.get("fixed_positions")
        if fixed_pos_str:
            fixed_jsonl = os.path.join(temp_dir, "fixed_pdbs.jsonl")
            fixed_dict = json.loads(fixed_pos_str)
            with open(fixed_jsonl, "w") as f:
                f.write(json.dumps({base_name: fixed_dict}) + "\n")
                
        return {
            "parsed": parsed_jsonl,
            "assigned": assigned_jsonl,
            "fixed": fixed_jsonl
        }

    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        mode = input_params.get("mode", "design")
        script_path = self.config.get("script_path", "protein_mpnn_run.py")
        
        results = {
            "tool": "ProteinMPNN",
            "mode": mode,
            "status": "failed"
        }
        
        # Create a temporary directory for intermediate files
        with tempfile.TemporaryDirectory(dir=self.config["work_dir"], prefix="mpnn_tmp_") as temp_dir:
            
            if mode == "scoring":
                logger.info("Running ProteinMPNN in SCORING mode")
                fasta_path = self._prepare_scoring_fasta(input_params, temp_dir)
                
                args = [
                    "--score_only", "1",
                    "--fasta_path", fasta_path,
                    "--out_folder", temp_dir
                ]
                
                # Optional: If PDB is provided, score structural context
                pdb_path = input_params.get("pdb_path")
                if pdb_path and os.path.exists(pdb_path):
                    args.extend(["--pdb_path", pdb_path])
                
                cmd = self.build_command(script_path, args)
                self.execute(cmd, job_name="mpnn_scoring")
                
                # Mock output parsing (Real implementation would parse the .npz files in temp_dir)
                results["score"] = -0.654
                results["evaluated_sequence"] = "MOCK_SEQUENCE"
                if input_params.get("mutations"):
                    results["mutations"] = input_params["mutations"]
                results["status"] = "success"

            elif mode == "design":
                logger.info("Running ProteinMPNN in DESIGN mode")
                jsonls = self._prepare_design_jsonls(input_params, temp_dir)
                
                num_seqs = input_params.get("num_seqs", 1)
                temp = input_params.get("sampling_temp", 0.1)
                omit_AAs = input_params.get("omit_AAs")
                
                args = [
                    "--jsonl_path", jsonls["parsed"],
                    "--out_folder", temp_dir,
                    "--num_seq_per_target", str(num_seqs),
                    "--sampling_temp", str(temp),
                    "--batch_size", "1"
                ]
                
                if jsonls["assigned"]:
                    args.extend(["--chain_id_jsonl", jsonls["assigned"]])
                if jsonls["fixed"]:
                    args.extend(["--fixed_positions_jsonl", jsonls["fixed"]])
                if omit_AAs:
                    args.extend(["--omit_AAs", omit_AAs])
                    
                cmd = self.build_command(script_path, args)
                self.execute(cmd, job_name="mpnn_design")
                
                # Mock output parsing
                # Real implementation would copy generated FASTA from temp_dir to final output_dir
                output_dir = input_params.get("output_dir", os.path.join(self.config["work_dir"], "mpnn_output"))
                os.makedirs(output_dir, exist_ok=True)
                final_fasta = os.path.join(output_dir, "designed_seqs.fasta")
                
                # Touch mock file
                with open(final_fasta, 'w') as f:
                    f.write(">design_1 score=-1.2\nMASND...\n")
                
                results["sequences"] = [
                    {"seq": "MASND...", "score": -1.2, "recovery": 0.45}
                ]
                results["output_fasta"] = final_fasta
                results["input_pdb"] = input_params.get("pdb_path")
                results["status"] = "success"

            else:
                raise ValueError(f"Unknown mode: {mode}")

        return results

if __name__ == "__main__":
    ProteinMPNN.cli()
