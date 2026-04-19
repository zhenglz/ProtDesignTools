
import os
import json
import logging
import argparse
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Union, List
from protdesigntools.core.task_manager import get_manager, TaskManager
from protdesigntools.core.sequence import Sequence
from protdesigntools.core.structure import Structure
from protdesigntools.core.config import config as global_config

logger = logging.getLogger(__name__)

class BaseTool(ABC):
    """
    Base class for all protein design tools.
    Supports both functional call (import) and script call (CLI).
    """
    
    def __init__(self, config_path: Optional[str] = None, global_config_path: Optional[str] = None, **kwargs):
        """
        Initialize the tool. 
        Configurations can be loaded from a JSON file or passed directly via kwargs.
        """
        self.tool_name = self.__class__.__name__
        
        # 0. Update global config if explicitly provided
        if global_config_path and os.path.exists(global_config_path):
            global_config.load_from_file(global_config_path)
            
        # Default tool-specific config (can be overridden by global tools section)
        self.config: Dict[str, Any] = {
            "python_env": "",
            "script_path": "",
            "binary_path": "",
            "tool_work_dir": self.tool_name.lower(), # Default subfolder name
            "exec_mode": "",  # Will fallback to global if empty
            "slurm_params": {}
        }
        
        # 1. Load tool specific defaults from the global config
        global_tools_config = global_config.get("tools", {})
        if self.tool_name in global_tools_config:
            self.config.update(global_tools_config[self.tool_name])
        
        # 2. Load from specific config file if provided (overrides global)
        if config_path and os.path.exists(config_path):
            with open(config_path, 'r') as f:
                user_config = json.load(f)
                self.config.update(user_config)
                
        # 3. Override with explicit kwargs (highest priority)
        self.config.update(kwargs)
        
        # Determine the output directory for the tool results
        self.output_dir = kwargs.get("output_dir", self.config.get("output_dir", f"./{self.tool_name.lower()}_output"))
        os.makedirs(self.output_dir, exist_ok=True)

    def get_exec_mode(self) -> str:
        """Get the execution mode (local or slurm), falling back to global config if not set locally."""
        mode = self.config.get("exec_mode")
        if not mode:
            mode = global_config.global_exec_mode
        return mode

    @abstractmethod
    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Main logic for the tool.
        Must be implemented by subclasses.
        """
        pass

    def __call__(self, **kwargs) -> Dict[str, Any]:
        """
        Allow the tool instance to be called directly like a function.
        """
        runtime_params = self.config.copy()
        runtime_params.update(kwargs)
        return self.run(runtime_params)

    def run_with_json(self, json_path: str) -> str:
        """Run the tool using a JSON configuration file and save the result to a JSON."""
        with open(json_path, 'r') as f:
            params = json.load(f)
        
        runtime_config = self.config.copy()
        runtime_config.update(params)
        
        results = self.run(runtime_config)
        
        output_json = kwargs.pop("output_json", os.path.join(self.output_dir, f"{self.tool_name}_output.json"))
        with open(output_json, 'w') as f:
            json.dump(results, f, indent=4)
        
        return output_json

    def build_command(self, script_or_bin: str, args: List[str]) -> str:
        """Helper to build execution command with appropriate python environment"""
        cmd_parts = []
        if self.config.get("python_env"):
            # Assuming conda
            cmd_parts.append(f"source activate {self.config['python_env']} &&")
            
        if script_or_bin.endswith(".py"):
            cmd_parts.append("python")
            
        cmd_parts.append(script_or_bin)
        cmd_parts.extend(args)
        return " ".join(cmd_parts)

    def execute(self, command: str, job_name: Optional[str] = None, **kwargs) -> Any:
        """
        Execute a command using the TaskManager.
        """
        job_name = job_name or f"{self.tool_name}_job"
        mode = self.get_exec_mode()
        
        # Setup Slurm params by merging global defaults with tool-specific params and runtime kwargs
        slurm_kwargs = global_config.get("default_slurm_params", {}).copy()
        slurm_kwargs.update(self.config.get("slurm_params", {}))
        slurm_kwargs.update(kwargs)
        print("Running cmd: ", cmd)
        manager = get_manager(mode, work_dir=self.output_dir, **slurm_kwargs)
        job_id = manager.submit(command, job_name, **slurm_kwargs)
        
        if kwargs.get("wait", True):
            manager.wait_for_jobs([job_id])
        return job_id

    @classmethod
    def get_cli_parser(cls) -> argparse.ArgumentParser:
        """
        Returns a base ArgumentParser for the tool.
        """
        parser = argparse.ArgumentParser(description=f"Run {cls.__name__} Tool")
        parser.add_argument("--global_config", type=str, help="Path to global JSON config file")
        parser.add_argument("--config", type=str, help="Path to input JSON config (optional)")
        parser.add_argument("--tool_work_dir", type=str, help="Subdirectory name for this tool's outputs")
        parser.add_argument("--output_dir", type=str, help="Specific output directory (overrides default nested path)")
        parser.add_argument("--exec_mode", type=str, choices=["local", "slurm"], help="Execution mode (overrides global)")
        parser.add_argument("--output_json", type=str, help="Path to save output JSON")
        return parser

    @classmethod
    def cli(cls):
        """Standard CLI entry point for the tool."""
        parser = cls.get_cli_parser()
        args, unknown = parser.parse_known_args()
        
        kwargs = {k: v for k, v in vars(args).items() if v is not None}
        
        config_path = kwargs.pop("config", None)
        global_config_path = kwargs.pop("global_config", None)
        output_json = kwargs.pop("output_json", None)
        
        tool = cls(config_path=config_path, global_config_path=global_config_path, **kwargs)
        
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
                
        results = tool(output_json=output_json, **input_params)
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=4)
            print(f"Results saved to {output_json}")
        else:
            print(json.dumps(results, indent=4))
