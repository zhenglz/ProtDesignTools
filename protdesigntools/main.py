
import argparse
import json
import sys
import os

# Add current directory to path so core and tools can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from protdesigntools.core.config import config as global_config
from protdesigntools.tools.proteinmpnn.tool import ProteinMPNN
from protdesigntools.tools.chai1.tool import Chai1
from protdesigntools.tools.esm2.tool import ESM2
from protdesigntools.tools.dlkcat.tool import DLKcat
from protdesigntools.tools.tmalign.tool import TMalign
from protdesigntools.tools.pythia.tool import Pythia
from protdesigntools.tools.autodock_vina.tool import AutoDockVina
from protdesigntools.tools.openmm.tool import OpenMMSimulation

def main():
    parser = argparse.ArgumentParser(description="ProtDesignTools: A Modular Protein Design Toolkit")
    
    # Global Configurations
    parser.add_argument("--global_config", type=str, help="Path to global JSON config file")
    parser.add_argument("--global_exec_mode", type=str, choices=["local", "slurm"], help="Global execution mode")
    
    # Tool Execution
    parser.add_argument("tool", type=str, help="Tool name to run (e.g., proteinmpnn, chai1, esm2, dlkcat, tmalign, pythia, vina, openmm)")
    parser.add_argument("--config", type=str, help="Optional: Path to tool-specific JSON config file (overrides global defaults)")
    
    args, unknown = parser.parse_known_args()
    
    # Setup Global Config
    if args.global_config:
        global_config.load_from_file(args.global_config)
    if args.global_exec_mode:
        global_config.update({"global_exec_mode": args.global_exec_mode})
    
    # Parse remaining dynamic parameters
    input_params = {}
    i = 0
    while i < len(unknown):
        if unknown[i].startswith("--"):
            key = unknown[i][2:]
            if i + 1 < len(unknown) and not unknown[i+1].startswith("--"):
                input_params[key] = unknown[i+1]
                i += 2
            else:
                input_params[key] = True
                i += 1
        else:
            i += 1
    
    # Registry of available tools
    tool_map = {
        "proteinmpnn": ProteinMPNN,
        "chai1": Chai1,
        "esm2": ESM2,
        "dlkcat": DLKcat,
        "tmalign": TMalign,
        "pythia": Pythia,
        "vina": AutoDockVina,
        "openmm": OpenMMSimulation,
        "rfdiffusion": RFDiffusion
    }
    
    if args.tool.lower() not in tool_map:
        print(f"Error: Tool '{args.tool}' not found. Available tools: {list(tool_map.keys())}")
        sys.exit(1)
        
    tool_class = tool_map[args.tool.lower()]
    tool_instance = tool_class(config_path=args.config)
    
    print(f"Running {args.tool} in work_dir {tool_instance.work_dir}...")
    output_path = tool_instance(output_json=None, **input_params)
    print(f"Done. Results:\n{json.dumps(output_path, indent=4)}")

if __name__ == "__main__":
    main()
