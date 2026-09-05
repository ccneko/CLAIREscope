"""CLAIREscope configuration loader and validator with user-override hierarchy."""
import os
import re
import yaml
from typing import Dict, Any, List, Tuple

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
    fallback_root_path = os.path.join(CONFIG_DIR, filename)
    if os.path.exists(fallback_root_path):
        return fallback_root_path
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

def save_user_project_config(project_key: str, project_dict: Dict[str, Any]) -> str:
    """Save or update a project definition in config/user/projects.yaml."""
    os.makedirs(USER_DIR, exist_ok=True)
    user_projects_file = os.path.join(USER_DIR, "projects.yaml")
    existing_data = load_yaml_config(user_projects_file, default={})
    if not isinstance(existing_data, dict):
        existing_data = {}
    existing_data[project_key] = project_dict
    with open(user_projects_file, "w", encoding="utf-8") as f:
        yaml.dump(existing_data, f, sort_keys=False, allow_unicode=True)
    return user_projects_file

def get_next_new_project_name(existing_projects: Dict[str, Any]) -> Tuple[str, str]:
    """
    Calculate the next available 'New Project N' name and identifier key.
    If 'New Project 1' exists, increments to 'New Project 2', etc.
    """
    existing_names = [p.get("name", "") for p in existing_projects.values() if isinstance(p, dict)]
    existing_keys = list(existing_projects.keys())
    
    i = 1
    while True:
        candidate_name = f"New Project {i}"
        candidate_key = f"PROJ_NEW_{i:03d}"
        name_conflict = any(candidate_name.lower() == str(name).lower() for name in existing_names)
        key_conflict = candidate_key in existing_keys
        if not name_conflict and not key_conflict:
            return candidate_name, candidate_key
        i += 1

def scan_project_datasets(proj_base: str, scan_subdirs: List[str] = None, max_depth: int = 4) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Find all .h5ad dataset files from the project root folder path and subdirectories.
    Prioritizes explicit scan_subdirs, then performs bounded tree search from proj_base.
    """
    active_datasets: Dict[str, str] = {}
    all_datasets: Dict[str, str] = {}
    found_paths = set()
    
    # Load dataset settings for hidden datasets
    dataset_cfg_file = get_config_file_path("dataset_config.yaml")
    cfg = load_yaml_config(dataset_cfg_file, default={})
    hidden_list = cfg.get("hidden_datasets", []) if isinstance(cfg, dict) else []
    
    if not proj_base or not os.path.exists(proj_base):
        return active_datasets, all_datasets
        
    candidate_dirs = []
    if scan_subdirs:
        for sub in scan_subdirs:
            p = os.path.abspath(os.path.join(proj_base, sub))
            if os.path.exists(p) and p not in candidate_dirs:
                candidate_dirs.append(p)
                
    if proj_base not in candidate_dirs:
        candidate_dirs.append(proj_base)
        
    # 1. First pass: scan explicit candidate directories
    for d in candidate_dirs:
        if not os.path.exists(d):
            continue
        try:
            for entry in sorted(os.listdir(d)):
                if entry.endswith(".h5ad") and not entry.startswith("."):
                    filepath = os.path.join(d, entry)
                    if filepath not in found_paths and os.path.isfile(filepath):
                        found_paths.add(filepath)
                        rel_dir = os.path.relpath(d, proj_base)
                        if rel_dir == ".":
                            tag = "root"
                        else:
                            tag = os.path.basename(d) if os.path.dirname(rel_dir) == "" else rel_dir.replace("\\", "/")
                        ds_name = f"{entry[:-5]} ({tag})"
                        all_datasets[ds_name] = filepath
                        if ds_name not in hidden_list:
                            active_datasets[ds_name] = filepath
        except Exception as e:
            print(f"Notice: scanning candidate dir {d}: {e}")

    # 2. Second pass: search tree from proj_base up to max_depth
    try:
        base_depth = len(os.path.abspath(proj_base).rstrip(os.sep).split(os.sep))
        for root, dirs, files in os.walk(proj_base):
            # Exclude non-data & system hidden directories
            dirs[:] = [d for d in dirs if not d.startswith(".") and d not in ["__pycache__", "node_modules", ".git", ".antigravity"]]
            cur_depth = len(os.path.abspath(root).rstrip(os.sep).split(os.sep))
            if cur_depth - base_depth > max_depth:
                dirs.clear()
                continue
            for f in sorted(files):
                if f.endswith(".h5ad") and not f.startswith("."):
                    filepath = os.path.join(root, f)
                    if filepath not in found_paths and os.path.isfile(filepath):
                        found_paths.add(filepath)
                        rel_dir = os.path.relpath(root, proj_base)
                        tag = "root" if rel_dir == "." else rel_dir.replace("\\", "/")
                        ds_name = f"{f[:-5]} ({tag})"
                        all_datasets[ds_name] = filepath
                        if ds_name not in hidden_list:
                            active_datasets[ds_name] = filepath
    except Exception as e:
        print(f"Notice: scanning tree under {proj_base}: {e}")

    return active_datasets, all_datasets

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
