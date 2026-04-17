
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
        self.config: Dict[str, Any] = {}
        if config_path and os.path.exists(config_path):
            with open(config_path, 'r') as f:
                self.config = json.load(f)
        
        self.tool_name = self.__class__.__name__

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
        
        results = self.run(params)
        
        output_json = params.get("output_json", "output.json")
        with open(output_json, 'w') as f:
            json.dump(results, f, indent=4)
        
        return output_json

    def execute(self, command: str, mode: str = "local", job_name: Optional[str] = None, **kwargs) -> Any:
        """
        Execute a command using the TaskManager.
        Useful for tools that call external binaries or scripts.
        """
        job_name = job_name or f"{self.tool_name}_job"
        manager = get_manager(mode, **kwargs)
        job_id = manager.submit(command, job_name, **kwargs)
        
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
