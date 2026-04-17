
import os
import json
import logging
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Union
from core.task_manager import get_manager, TaskManager

logger = logging.getLogger(__name__)

class BaseTool(ABC):
    """
    Base class for all protein design tools.
    Supports both functional call (import) and script call (CLI).
    """
    
    def __init__(self, config_path: Optional[str] = None):
        self.config: Dict[str, Any] = {
            "python_env": "",
            "script_path": "",
            "binary_path": "",
            "work_dir": "./work_dir",
            "exec_mode": "local",  # local or slurm
            "slurm_params": {}
        }
        
        # Override default config with user provided config
        if config_path and os.path.exists(config_path):
            with open(config_path, 'r') as f:
                user_config = json.load(f)
                self.config.update(user_config)
        
        self.tool_name = self.__class__.__name__
        os.makedirs(self.config["work_dir"], exist_ok=True)

    @abstractmethod
    def run(self, input_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Main entry point for the tool.
        input_params: Dictionary containing all inputs (paths, parameters, etc.)
        returns: Dictionary containing all outputs (scores, paths to files, etc.)
        """
        pass

    def run_with_json(self, json_path: str) -> str:
        """Run the tool using a JSON configuration file and save the result to a JSON."""
        with open(json_path, 'r') as f:
            params = json.load(f)
        
        # Merge class-level config with runtime params
        # This allows params to override tool defaults (like exec_mode)
        runtime_config = self.config.copy()
        runtime_config.update(params)
        
        results = self.run(runtime_config)
        
        output_json = runtime_config.get("output_json", os.path.join(self.config["work_dir"], f"{self.tool_name}_output.json"))
        with open(output_json, 'w') as f:
            json.dump(results, f, indent=4)
        
        return output_json

    def build_command(self, script_or_bin: str, args: List[str]) -> str:
        """Helper to build execution command with appropriate python environment"""
        cmd_parts = []
        if self.config.get("python_env"):
            # Assuming conda or virtualenv is activated via this path
            cmd_parts.append(f"source activate {self.config['python_env']} && python")
        elif script_or_bin.endswith(".py"):
            cmd_parts.append("python")
            
        cmd_parts.append(script_or_bin)
        cmd_parts.extend(args)
        return " ".join(cmd_parts)

    def execute(self, command: str, job_name: Optional[str] = None, **kwargs) -> Any:
        """
        Execute a command using the TaskManager.
        Useful for tools that call external binaries or scripts.
        """
        job_name = job_name or f"{self.tool_name}_job"
        mode = self.config.get("exec_mode", "local")
        
        # Slurm parameters
        slurm_kwargs = self.config.get("slurm_params", {})
        slurm_kwargs.update(kwargs)
        
        manager = get_manager(mode, work_dir=self.config["work_dir"], **slurm_kwargs)
        job_id = manager.submit(command, job_name, **slurm_kwargs)
        
        # If running synchronously
        if kwargs.get("wait", True):
            manager.wait_for_jobs([job_id])
            # In a real tool, we would check the status and handle logs/errors
            return job_id
        return job_id

    @classmethod
    def cli(cls):
        """CLI entry point for the tool."""
        import argparse
        parser = argparse.ArgumentParser(description=f"Run {cls.__name__}")
        parser.add_argument("--config", type=str, required=True, help="Path to input JSON config")
        args = parser.parse_args()
        
        tool = cls()
        tool.run_with_json(args.config)
