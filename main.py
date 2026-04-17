
import argparse
import json
import sys
import os

# Add current directory to path so core and tools can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from tools.proteinmpnn.tool import ProteinMPNN
# Import other tools as they are implemented

def main():
    parser = argparse.ArgumentParser(description="ProtDesignTools: A Modular Protein Design Toolkit")
    parser.add_argument("tool", type=str, help="Tool name to run (e.g., proteinmpnn)")
    parser.add_argument("--config", type=str, required=True, help="Path to JSON config file")
    
    args = parser.parse_args()
    
    tool_map = {
        "proteinmpnn": ProteinMPNN,
        # "af3": AlphaFold3,
        # "chai1": Chai1,
    }
    
    if args.tool.lower() not in tool_map:
        print(f"Error: Tool '{args.tool}' not found. Available tools: {list(tool_map.keys())}")
        sys.exit(1)
        
    tool_class = tool_map[args.tool.lower()]
    tool_instance = tool_class()
    
    print(f"Running {args.tool} with config {args.config}...")
    output_path = tool_instance.run_with_json(args.config)
    print(f"Done. Results saved to {output_path}")

if __name__ == "__main__":
    main()
