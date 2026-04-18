
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
        parser.add_argument("--output_dir", type=str, help="Output directory")
        return parser

    def _prepare_inference_json(self, input_params: Dict[str, Any], temp_dir: str) -> str:
        pdb_path = input_params.get("pdb_path")
        contig = input_params.get("contig")
        length = input_params.get("length", "")
        fixed_atoms_str = input_params.get("fixed_atoms")
        
        fixed_atoms = {}
        if fixed_atoms_str:
            fixed_atoms = json.loads(fixed_atoms_str)
            
        # Fallback to auto-detect if no contig is provided (stub for Structure method)
        if not contig:
            # Here you would typically call a helper to detect gaps in the PDB
            # For now we'll just set a generic contig for unconditional generation or 
            # assume the user provides it.
            logger.warning("No contig provided. Generating unconditional scaffold of length 100.")
            contig = "100"
            
        json_data = {
            "rfd3_design_task": {
                "input": os.path.abspath(pdb_path),
                "contig": contig,
                "length": length,
                "select_fixed_atoms": fixed_atoms,
                "is_non_loopy": True,
                "dialect": 2
            }
        }
        
        json_path = os.path.join(temp_dir, "inference.json")
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
            
        return json_path

    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        logger.info("Running RFDiffusion (RFD3) Design")
        
        output_dir = input_params.get("output_dir", os.path.join(self.work_dir, "output"))
        os.makedirs(output_dir, exist_ok=True)
        
        # We assume script_path is actually the source_rfd3.sh script
        source_script = self.config.get("script_path", "source_rfd3.sh")
        
        results = {
            "tool": "RFDiffusion",
            "status": "failed"
        }
        
        with tempfile.TemporaryDirectory(dir=self.work_dir, prefix="rfd_tmp_") as temp_dir:
            json_path = self._prepare_inference_json(input_params, temp_dir)
            num_designs = input_params.get("num_designs", 1)
            
            # The command from zDBProd runs `source script.sh && rfd3 design out_dir=... inputs=...`
            # Since this is a chained command, we build it directly
            cmd = f"source {source_script} && rfd3 design out_dir={os.path.abspath(temp_dir)} inputs={os.path.abspath(json_path)}"
            if num_designs > 1:
                cmd += f" inference.num_designs={num_designs}"
                
            job_id = self.execute(cmd, job_name="rfd_design")
            
            # Parse output: look for generated .cif.gz files
            generated_files = glob.glob(os.path.join(temp_dir, "*.cif.gz"))
            if not generated_files:
                logger.error(f"No .cif.gz output found in {temp_dir}")
                return results
                
            output_files = []
            for gf in generated_files:
                dest = os.path.join(output_dir, os.path.basename(gf))
                shutil.copy(gf, dest)
                output_files.append(dest)
                
            results["generated_files"] = output_files
            results["job_id"] = job_id
            results["output_dir"] = output_dir
            results["status"] = "success"

        return results

if __name__ == "__main__":
    RFDiffusion.cli()
