
import os
import json
from typing import Dict, Any, Optional

class GlobalConfig:
    """
    Singleton class to manage global configurations across all tools.
    """
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(GlobalConfig, cls).__new__(cls)
            cls._instance.config = {
                "global_work_dir": "./work_dir",
                "global_exec_mode": "local", # can be overridden by specific tools
                "default_slurm_partition": "AMD"
            }
        return cls._instance

    def update(self, config_dict: Dict[str, Any]):
        """Update global settings from a dictionary."""
        self.config.update(config_dict)
        
    def load_from_file(self, json_path: str):
        """Update global settings from a JSON file."""
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                self.update(json.load(f))

    def get(self, key: str, default: Any = None) -> Any:
        return self.config.get(key, default)

    @property
    def global_work_dir(self) -> str:
        return self.config["global_work_dir"]

    @property
    def global_exec_mode(self) -> str:
        return self.config["global_exec_mode"]

# Global instance
config = GlobalConfig()
