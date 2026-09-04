"""CLAIREscope configuration loader and validator with user-override hierarchy."""
import os
import yaml
from typing import Dict, Any, List

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(PACKAGE_DIR, "config")
DEFAULTS_DIR = os.path.join(CONFIG_DIR, "defaults")
USER_DIR = os.path.join(CONFIG_DIR, "user")

def get_platform_path(win_path: str, wsl_path: str) -> str:
    """Resolve OS-appropriate directory path between Windows and WSL."""
    if os.path.exists(win_path):
        return win_path
    if os.path.exists(wsl_path):
        return wsl_path
    return win_path if os.name == 'nt' else wsl_path

def get_config_file_path(filename: str) -> str:
    """Resolve configuration file path: checks config/user/ first, then falls back to config/defaults/."""
    user_path = os.path.join(USER_DIR, filename)
    if os.path.exists(user_path) and os.path.getsize(user_path) > 0:
        return user_path
    default_path = os.path.join(DEFAULTS_DIR, filename)
    if os.path.exists(default_path):
        return default_path
    legacy_path = os.path.join(CONFIG_DIR, filename)
    if os.path.exists(legacy_path):
        return legacy_path
    return default_path

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
    """Load project registry (user override or defaults)."""
    cfg_file = get_config_file_path("projects.yaml")
    return load_yaml_config(cfg_file, default={})

def load_settings_config() -> Dict[str, Any]:
    """Load application global UI and plotting settings."""
    cfg_file = get_config_file_path("settings.yaml")
    return load_yaml_config(cfg_file, default={})

def load_signatures_config() -> Dict[str, Any]:
    """Load curated gene signatures and pathway panels."""
    cfg_file = get_config_file_path("signatures.yaml")
    return load_yaml_config(cfg_file, default={})

def load_pathways_config() -> Dict[str, List[str]]:
    """Load curated biological pathways database for ORA."""
    cfg_file = get_config_file_path("pathways.yaml")
    return load_yaml_config(cfg_file, default={})

def load_markers_config() -> Dict[str, Dict[str, List[str]]]:
    """Load canonical cell-type marker dictionary."""
    cfg_file = get_config_file_path("markers.yaml")
    return load_yaml_config(cfg_file, default={})

def load_css_styles() -> str:
    """Load external CSS stylesheet (user override or defaults)."""
    css_file = get_config_file_path("style.css")
    if os.path.exists(css_file):
        with open(css_file, "r", encoding="utf-8") as f:
            return f.read()
    return ""
