
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
                "global_exec_mode": "local",
                "default_slurm_params": {}
            }
            
            # Auto-load the default global config from data/global_config.json if it exists
            default_config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "global_config.json")
            cls._instance.load_from_file(default_config_path)
            
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
    def global_exec_mode(self) -> str:
        return self.config["global_exec_mode"]

# Global instance
config = GlobalConfig()
