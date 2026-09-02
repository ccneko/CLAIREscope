"""CLAIREscope configuration loader and validator."""
import os
import yaml
from typing import Dict, Any, List

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(PACKAGE_DIR, "config")

def get_platform_path(win_path: str, wsl_path: str) -> str:
    """Resolve OS-appropriate directory path between Windows and WSL."""
    if os.path.exists(win_path):
        return win_path
    if os.path.exists(wsl_path):
        return wsl_path
    return win_path if os.name == 'nt' else wsl_path

def load_yaml_config(file_path: str, default: Any = None) -> Any:
    """Safely load a YAML configuration file with fallback default."""
    if not os.path.exists(file_path):
        return default if default is not None else {}
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
            return data if data is not None else default
    except Exception as e:
        print(f"Warning: Failed to load YAML at {file_path}: {e}")
        return default if default is not None else {}

def load_projects_config() -> Dict[str, Any]:
    """Load the multi-project registry configuration."""
    cfg_file = os.path.join(CONFIG_DIR, "projects.yaml")
    return load_yaml_config(cfg_file, default={})

def load_settings_config() -> Dict[str, Any]:
    """Load application global UI, typography, and plotting settings."""
    cfg_file = os.path.join(CONFIG_DIR, "settings.yaml")
    return load_yaml_config(cfg_file, default={})

def load_signatures_config() -> Dict[str, Any]:
    """Load curated gene signatures and pathway panels."""
    cfg_file = os.path.join(CONFIG_DIR, "signatures.yaml")
    return load_yaml_config(cfg_file, default={})
