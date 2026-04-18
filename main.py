
import argparse
import json
import sys
import os

# Add current directory to path so core and tools can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from core.config import config as global_config
from tools.proteinmpnn.tool import ProteinMPNN
from tools.chai1.tool import Chai1
from tools.esm2.tool import ESM2
from tools.dlkcat.tool import DLKcat
from tools.tmalign.tool import TMalign
from tools.pythia.tool import Pythia
from tools.autodock_vina.tool import AutoDockVina
from tools.openmm.tool import OpenMMSimulation

def main():
    parser = argparse.ArgumentParser(description="ProtDesignTools: A Modular Protein Design Toolkit")
    
    # Global Configurations
    parser.add_argument("--global_config", type=str, help="Path to global JSON config file")
    parser.add_argument("--global_work_dir", type=str, help="Global working directory")
    parser.add_argument("--global_exec_mode", type=str, choices=["local", "slurm"], help="Global execution mode")
    
    # Tool Execution
    parser.add_argument("tool", type=str, help="Tool name to run (e.g., proteinmpnn, chai1, esm2, dlkcat, tmalign, pythia, vina, openmm)")
    parser.add_argument("--config", type=str, required=True, help="Path to tool-specific JSON config file")
    
    args = parser.parse_args()
    
    # Setup Global Config
    if args.global_config:
        global_config.load_from_file(args.global_config)
    if args.global_work_dir:
        global_config.update({"global_work_dir": args.global_work_dir})
    if args.global_exec_mode:
        global_config.update({"global_exec_mode": args.global_exec_mode})
    
    # Registry of available tools
    tool_map = {
        "proteinmpnn": ProteinMPNN,
        "chai1": Chai1,
        "esm2": ESM2,
        "dlkcat": DLKcat,
        "tmalign": TMalign,
        "pythia": Pythia,
        "vina": AutoDockVina,
        "openmm": OpenMMSimulation
    }
    
    if args.tool.lower() not in tool_map:
        print(f"Error: Tool '{args.tool}' not found. Available tools: {list(tool_map.keys())}")
        sys.exit(1)
        
    tool_class = tool_map[args.tool.lower()]
    tool_instance = tool_class()
    
    print(f"Running {args.tool} with config {args.config} in work_dir {tool_instance.work_dir}...")
    output_path = tool_instance.run_with_json(args.config)
    print(f"Done. Results saved to {output_path}")

if __name__ == "__main__":
    main()
