
import os
import argparse
import logging
import json
import tempfile
import glob
import shutil
from typing import Dict, Any, List, Optional
from protdesigntools.core.base_tool import BaseTool
from protdesigntools.core.structure import Structure

logger = logging.getLogger(__name__)

class RFDiffusion(BaseTool):
    """
    RFDiffusion (RFD3) for protein structure design.
    Inputs: PDB path, contig string (e.g. '10-100'), or use auto-detection.
    Outputs: Designed structure (.cif.gz or .pdb).
    """

    @classmethod
    def get_cli_parser(cls) -> argparse.ArgumentParser:
        parser = super().get_cli_parser()
        parser.add_argument("--pdb_path", type=str, required=True, help="Input PDB file")
        parser.add_argument("--contig", type=str, help="Contig string defining design constraints (e.g., '10-100')")
        parser.add_argument("--length", type=str, help="Length constraint (e.g., '50-100')")
        parser.add_argument("--fixed_atoms", type=str, help="JSON dict of fixed atoms")
        parser.add_argument("--num_designs", type=int, default=1, help="Number of designs to generate")
        parser.add_argument("--job_key", type=str, default="dsDNA_protein_design", help="JSON key for the job (e.g. dsDNA_protein_design)")
        parser.add_argument("--output_dir", type=str, help="Output directory")
        return parser

    def _prepare_inference_json(self, input_params: Dict[str, Any], output_dir: str) -> str:
        pdb_path = input_params.get("pdb_path")
        contig = input_params.get("contig")
        length = input_params.get("length", "")
        fixed_atoms_str = input_params.get("fixed_atoms")
        
        fixed_atoms = {}
        if fixed_atoms_str:
            if isinstance(fixed_atoms_str, str):
                fixed_atoms = json.loads(fixed_atoms_str)
            elif isinstance(fixed_atoms_str, dict):
                fixed_atoms = fixed_atoms_str
            
        # Fallback to auto-detect if no contig is provided
        if not contig:
            logger.warning("No contig provided. Generating unconditional scaffold of length 100.")
            contig = "100"
            
        job_key = input_params.get("job_key", "dsDNA_protein_design")
            
        json_data = {
            job_key: {
                "input": os.path.abspath(pdb_path),
                "contig": contig,
                "length": length,
                "select_fixed_atoms": fixed_atoms,
                "is_non_loopy": True,
                "dialect": 2
            }
        }
        
        json_path = os.path.join(output_dir, "inference.json")
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
            
        return json_path

    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        logger.info("Running RFDiffusion (RFD3) Design")
        
        if "output_dir" in input_params:
            self.output_dir = input_params["output_dir"]
        os.makedirs(self.output_dir, exist_ok=True)
        
        # We assume script_path is actually the source_rfd3.sh script
        source_script = self.config.get("script_path", "/sugon_store/zhengliangzhen/apps/RFDiffusion3/source_rfd3.sh")
        
        results = {
            "tool": "RFDiffusion",
            "status": "failed"
        }
        
        json_path = self._prepare_inference_json(input_params, self.output_dir)
        num_designs = input_params.get("num_designs", 1)
        
        cmd = f"source {source_script} && rfd3 design out_dir={os.path.abspath(self.output_dir)} inputs={os.path.abspath(json_path)}"
        if num_designs > 1:
            cmd += f" inference.num_designs={num_designs}"
            
        job_id = self.execute(cmd, job_name="rfd_design")
        
        # Parse output: look for generated .cif.gz files
        generated_files = glob.glob(os.path.join(self.output_dir, "*.cif.gz"))
        if not generated_files:
            logger.error(f"No .cif.gz output found in {self.output_dir}")
            return results
            
        results["generated_files"] = generated_files
        results["job_id"] = job_id
        results["output_dir"] = self.output_dir
        results["status"] = "success"

        return results

if __name__ == "__main__":
    RFDiffusion.cli()
