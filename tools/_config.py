#!/usr/bin/env python3
"""
Shared configuration loader for ProtDesignTools.

Loads tool dependency paths from a JSON config file, with fallback
to the repo's default data/config.json.
"""

import json
import os


def find_repo_root():
    """Find the repository root by walking up from this file looking for .git."""
    current = os.path.dirname(os.path.abspath(__file__))
    for _ in range(10):
        if os.path.exists(os.path.join(current, ".git")):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    # Fallback: assume repo root is the parent of tools/
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_config(config_path=None):
    """Load tool dependency configuration from a JSON file.

    Resolution order:
      1. If `config_path` is given and the file exists, use it.
      2. Otherwise try ``<repo_root>/data/config.json``.
      3. If neither exists, raise FileNotFoundError.

    Parameters
    ----------
    config_path : str or None
        Explicit path to a config JSON file.

    Returns
    -------
    dict
        The parsed configuration dictionary.
    """
    # 1. Explicit path
    if config_path and os.path.exists(config_path):
        with open(config_path) as f:
            return json.load(f)
    if config_path:
        raise FileNotFoundError(
            f"Config file not found: {config_path}"
        )

    # 2. Default repo config
    repo_root = find_repo_root()
    default_path = os.path.join(repo_root, "data", "config.json")
    if os.path.exists(default_path):
        with open(default_path) as f:
            return json.load(f)

    # 3. Nothing works
    raise FileNotFoundError(
        "No tool config file found. Create data/config.json in the repo root "
        "or pass an explicit --config path."
    )


if __name__ == "__main__":
    # Quick self-test when run directly
    try:
        cfg = load_config()
        print("Config loaded successfully:")
        print(json.dumps(cfg, indent=2))
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
