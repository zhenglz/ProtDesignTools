
import os
import argparse
import logging
import json
import tempfile
import shutil
import glob
import numpy as np
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
      - Calculates fitness score by reading .npz output.
      
    Design Mode:
      - Supports specifying chains to design.
      - Supports fixing specific positions or defining amino acid lists.
      - Parses .fa output to extract designed sequences and metrics.
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
        parser.add_argument("--omit_AAs", type=str, help="String of amino acids to omit globally (e.g. 'CX')")
        parser.add_argument("--omit_AA_per_pos", type=str, help="JSON string mapping chain to dict of pos: omit_AAs (e.g. '{\"A\": {\"1\": \"C\", \"2\": \"FWY\"}}')")
        parser.add_argument("--bias_AA_per_pos", type=str, help="JSON string mapping chain to dict of pos: dict of AA: bias (e.g. '{\"A\": {\"1\": {\"A\": 10.0}}}')")
        parser.add_argument("--tied_positions", type=str, help="JSON string for symmetric design. e.g. '{\"A\": [1,2], \"B\": [1,2]}'")
        parser.add_argument("--num_seqs", type=int, default=1, help="Number of sequences to design")
        parser.add_argument("--sampling_temp", type=float, default=0.1, help="Sampling temperature")
        
        return parser

    def _prepare_scoring_fasta(self, input_params: Dict[str, Any], temp_dir: str) -> str:
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
            final_seq = struct.get_sequence() 
            header += f" source={os.path.basename(pdb_path)} mutations={mutations}"
            
        else:
            raise ValueError("Scoring mode requires either fasta_path, sequence, or pdb_path.")

        with open(fasta_out, "w") as f:
            f.write(f"{header}\n{final_seq}\n")
            
        return fasta_out

    def _prepare_design_jsonls(self, input_params: Dict[str, Any], temp_dir: str) -> Dict[str, str]:
        pdb_path = input_params.get("pdb_path")
        if not pdb_path or not os.path.exists(pdb_path):
            raise ValueError("Design mode requires a valid pdb_path.")
            
        base_name = os.path.basename(pdb_path).split('.')[0]
        
        # 1. Parsed PDBs (Mocking the structure parsing for the wrapper)
        parsed_jsonl = os.path.join(temp_dir, "parsed_pdbs.jsonl")
        struct = Structure(pdb_path)
        parsed_data = {"name": base_name}
        for chain_id, seq in struct.chains.items():
            parsed_data[f"seq_chain_{chain_id}"] = seq
            # coords should ideally be included here if running full MPNN natively
        with open(parsed_jsonl, "w") as f:
            f.write(json.dumps(parsed_data) + "\n")
            
        # 2. Assigned Chains
        assigned_jsonl = ""
        design_chains = input_params.get("design_chains")
        if design_chains:
            assigned_jsonl = os.path.join(temp_dir, "assigned_pdbs.jsonl")
            chains = [c.strip() for c in design_chains.split(",")]
            all_chains = list(struct.chains.keys())
            fixed_chains = [c for c in all_chains if c not in chains]
            with open(assigned_jsonl, "w") as f:
                f.write(json.dumps({base_name: [chains, fixed_chains]}) + "\n")
                
        # 3. Fixed Positions
        fixed_jsonl = ""
        fixed_pos_str = input_params.get("fixed_positions")
        if fixed_pos_str:
            fixed_jsonl = os.path.join(temp_dir, "fixed_pdbs.jsonl")
            fixed_dict = json.loads(fixed_pos_str)
            with open(fixed_jsonl, "w") as f:
                f.write(json.dumps({base_name: fixed_dict}) + "\n")
                
        # 4. Omit AA per position
        omit_aa_jsonl = ""
        omit_aa_pos_str = input_params.get("omit_AA_per_pos")
        if omit_aa_pos_str:
            omit_aa_jsonl = os.path.join(temp_dir, "omit_AA.jsonl")
            omit_aa_dict = json.loads(omit_aa_pos_str)
            with open(omit_aa_jsonl, "w") as f:
                f.write(json.dumps({base_name: omit_aa_dict}) + "\n")
                
        # 5. Bias AA per position
        bias_aa_jsonl = ""
        bias_aa_pos_str = input_params.get("bias_AA_per_pos")
        if bias_aa_pos_str:
            bias_aa_jsonl = os.path.join(temp_dir, "bias_AA.jsonl")
            bias_aa_dict = json.loads(bias_aa_pos_str)
            with open(bias_aa_jsonl, "w") as f:
                f.write(json.dumps({base_name: bias_aa_dict}) + "\n")
                
        # 6. Tied positions (Symmetric design)
        tied_jsonl = ""
        tied_pos_str = input_params.get("tied_positions")
        if tied_pos_str:
            tied_jsonl = os.path.join(temp_dir, "tied_pdbs.jsonl")
            tied_dict = json.loads(tied_pos_str)
            with open(tied_jsonl, "w") as f:
                f.write(json.dumps({base_name: tied_dict}) + "\n")
                
        return {
            "parsed": parsed_jsonl,
            "assigned": assigned_jsonl,
            "fixed": fixed_jsonl,
            "omit_aa_pos": omit_aa_jsonl,
            "bias_aa_pos": bias_aa_jsonl,
            "tied": tied_jsonl
        }

    def _parse_scoring_output(self, temp_dir: str) -> float:
        """Parse the .npz file generated by ProteinMPNN in scoring mode."""
        npz_files = glob.glob(os.path.join(temp_dir, "score_only", "*.npz"))
        if not npz_files:
            raise RuntimeError(f"No .npz output found in {temp_dir}/score_only")
        
        # Read the first (and usually only) npz file
        npz_path = npz_files[0]
        try:
            dat = np.load(npz_path, allow_pickle=True)
            if 'score' in dat:
                return float(np.mean(dat['score']))
            else:
                logger.warning(f"'score' array not found in {npz_path}")
                return 0.0
        except Exception as e:
            logger.error(f"Error reading {npz_path}: {e}")
            raise

    def _parse_design_output(self, temp_dir: str) -> List[Dict[str, Any]]:
        """Parse the .fa file generated by ProteinMPNN in design mode."""
        seqs_dir = os.path.join(temp_dir, "seqs")
        fasta_files = glob.glob(os.path.join(seqs_dir, "*.fa"))
        
        if not fasta_files:
            raise RuntimeError(f"No .fa output found in {seqs_dir}")
            
        fasta_path = fasta_files[0]
        results = []
        
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            
        current_header = ""
        current_seq = ""
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    # Parse the previous entry
                    entry = self._parse_design_header(current_header, current_seq)
                    if entry:
                        results.append(entry)
                current_header = line
                current_seq = ""
            else:
                current_seq += line
                
        # Parse the last entry
        if current_header:
            entry = self._parse_design_header(current_header, current_seq)
            if entry:
                results.append(entry)
                
        return results

    def _parse_design_header(self, header: str, seq: str) -> Optional[Dict[str, Any]]:
        """
        Parse ProteinMPNN design fasta headers.
        Example 1 (WT): >protein, score=2.1936, fixed_chains=[], designed_chains=['A'], model_name=v_48_020
        Example 2 (Design): >T=0.1, sample=0, score=1.9629, seq_recovery=0.3333
        """
        header = header[1:] # Remove '>'
        entry = {"seq": seq}
        
        # We generally want to extract designed sequences, which have T=...
        if "T=" in header and "sample=" in header:
            parts = [p.strip() for p in header.split(",")]
            for part in parts:
                if "=" in part:
                    k, v = part.split("=", 1)
                    try:
                        entry[k] = float(v) if "." in v else int(v)
                    except ValueError:
                        entry[k] = v
            entry["type"] = "design"
            return entry
        elif "score=" in header and "model_name=" in header:
            # WT sequence info
            parts = [p.strip() for p in header.split(",")]
            for part in parts:
                if "=" in part:
                    k, v = part.split("=", 1)
                    if k == "score":
                        entry["score"] = float(v)
            entry["type"] = "native"
            return entry
            
        return None

    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        mode = input_params.get("mode", "design")
        script_path = self.config.get("script_path", "protein_mpnn_run.py")
        
        results = {
            "tool": "ProteinMPNN",
            "mode": mode,
            "status": "failed"
        }
        
        with tempfile.TemporaryDirectory(dir=self.work_dir, prefix="mpnn_tmp_") as temp_dir:
            
            if mode == "scoring":
                logger.info("Running ProteinMPNN in SCORING mode")
                fasta_path = self._prepare_scoring_fasta(input_params, temp_dir)
                
                args = [
                    "--score_only", "1",
                    "--fasta_path", fasta_path,
                    "--out_folder", temp_dir
                ]
                
                pdb_path = input_params.get("pdb_path")
                if pdb_path and os.path.exists(pdb_path):
                    args.extend(["--pdb_path", pdb_path])
                
                cmd = self.build_command(script_path, args)
                self.execute(cmd, job_name="mpnn_scoring")
                
                # Parse the .npz output
                score = self._parse_scoring_output(temp_dir)
                
                results["score"] = score
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
                if jsonls["omit_aa_pos"]:
                    args.extend(["--omit_AA_jsonl", jsonls["omit_aa_pos"]])
                if jsonls["bias_aa_pos"]:
                    args.extend(["--bias_AA_jsonl", jsonls["bias_aa_pos"]])
                if jsonls["tied"]:
                    args.extend(["--tied_positions_jsonl", jsonls["tied"]])
                if omit_AAs:
                    args.extend(["--omit_AAs", omit_AAs])
                    
                cmd = self.build_command(script_path, args)
                self.execute(cmd, job_name="mpnn_design")
                
                # Parse the .fa output
                parsed_seqs = self._parse_design_output(temp_dir)
                
                # Copy the generated FASTA to the final output_dir
                output_dir = input_params.get("output_dir", os.path.join(self.work_dir, "output"))
                os.makedirs(output_dir, exist_ok=True)
                
                seqs_dir = os.path.join(temp_dir, "seqs")
                fasta_files = glob.glob(os.path.join(seqs_dir, "*.fa"))
                final_fasta = ""
                if fasta_files:
                    final_fasta = os.path.join(output_dir, os.path.basename(fasta_files[0]))
                    shutil.copy(fasta_files[0], final_fasta)
                
                results["sequences"] = parsed_seqs
                if final_fasta:
                    results["output_fasta"] = final_fasta
                results["input_pdb"] = input_params.get("pdb_path")
                results["status"] = "success"

            else:
                raise ValueError(f"Unknown mode: {mode}")

        return results

if __name__ == "__main__":
    ProteinMPNN.cli()
