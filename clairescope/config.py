"""CLAIREscope configuration loader and validator with user-override hierarchy."""
import os
import re
import yaml
from typing import Dict, Any, List, Tuple, Optional

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

def save_settings_config(settings_dict: Dict[str, Any]) -> str:
    """Save user modified application settings into config/user/settings.yaml."""
    os.makedirs(USER_DIR, exist_ok=True)
    target_file = os.path.join(USER_DIR, "settings.yaml")
    with open(target_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(settings_dict, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    return target_file

def load_dataset_config() -> Dict[str, Any]:
    """Load dataset preferences (default dataset, hidden datasets, custom datasets)."""
    cfg_file = get_config_file_path("dataset_config.yaml")
    return load_yaml_config(cfg_file, default={})

def save_dataset_config(dataset_cfg: Dict[str, Any]) -> str:
    """Save dataset preferences into config/user/dataset_config.yaml."""
    os.makedirs(USER_DIR, exist_ok=True)
    target_file = os.path.join(USER_DIR, "dataset_config.yaml")
    with open(target_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(dataset_cfg, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    return target_file

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

def load_annotation_colors(col_name: Optional[str] = None) -> Dict[str, Any]:
    """Load user or dataset custom annotation category colors."""
    ds_cfg = load_dataset_config()
    colors = ds_cfg.get("annotation_colors", {})
    if not colors:
        app_settings = load_settings_config()
        colors = app_settings.get("annotation_colors", {})
    if col_name:
        return colors.get(col_name, {})
    return colors

def save_annotation_colors(col_name: str, color_map: Dict[str, str]) -> str:
    """Persist custom annotation category colors to user dataset_config.yaml."""
    ds_cfg = load_dataset_config()
    all_colors = ds_cfg.setdefault("annotation_colors", {})
    all_colors[col_name] = color_map
    return save_dataset_config(ds_cfg)

def import_dataset_from_yaml(yaml_source: Any, persist: bool = True) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Import a dataset definition and optionally its parent project from a YAML file or string.
    Supports:
      1. Single dataset object under 'dataset' key.
      2. List of dataset objects under 'datasets' key.
      3. Associated project metadata under 'project' key.
      4. Default preferences under 'settings' key.
      5. Custom annotation colors under 'annotation_colors'.
    """
    data = None
    if isinstance(yaml_source, dict):
        data = yaml_source
    elif isinstance(yaml_source, str):
        clean_src = yaml_source.strip()
        if ("\n" not in clean_src) and os.path.isfile(clean_src):
            try:
                with open(clean_src, "r", encoding="utf-8") as f:
                    data = yaml.safe_load(f)
            except Exception as e:
                return False, f"Failed to read YAML file at {clean_src}: {e}", {}
        else:
            try:
                data = yaml.safe_load(yaml_source)
            except Exception as e:
                return False, f"Failed to parse YAML text: {e}", {}
    
    if not isinstance(data, dict):
        return False, "Invalid YAML format: root configuration must be a dictionary.", {}

    projects = load_projects_config()
    ds_cfg = load_dataset_config()
    app_settings = load_settings_config()

    imported_datasets = []
    target_project_key = None

    # Process settings block
    if "settings" in data and isinstance(data["settings"], dict):
        s_block = data["settings"]
        if "default_project" in s_block and s_block["default_project"]:
            target_project_key = s_block["default_project"]
            app_settings["default_project"] = s_block["default_project"]
            ds_cfg["default_project"] = s_block["default_project"]
        if "default_dataset" in s_block and s_block["default_dataset"]:
            ds_cfg["default_dataset"] = s_block["default_dataset"]
        if "default_annotation_col" in s_block and s_block["default_annotation_col"]:
            app_settings["default_annotation_col"] = s_block["default_annotation_col"]

    # Process project block
    if "project" in data and isinstance(data["project"], dict):
        p_block = data["project"]
        p_key = p_block.get("key") or p_block.get("id") or target_project_key
        if p_key:
            target_project_key = p_key
            if p_key not in projects:
                projects[p_key] = {}
            for k, v in p_block.items():
                if k != "key":
                    projects[p_key][k] = v

    # Extract dataset entries
    ds_entries = []
    if "dataset" in data and isinstance(data["dataset"], dict):
        ds_entries.append(data["dataset"])
    if "datasets" in data and isinstance(data["datasets"], list):
        for d in data["datasets"]:
            if isinstance(d, dict):
                ds_entries.append(d)

    # Process dataset entries
    for ds in ds_entries:
        ds_name = ds.get("name")
        ds_path = ds.get("path") or ds.get("wsl_path") or ds.get("win_path")
        p_key = ds.get("project_key") or ds.get("project_id") or target_project_key
        
        target_proj = None
        if p_key and p_key in projects:
            target_proj = projects[p_key]
            target_project_key = p_key
        else:
            for pk, pval in projects.items():
                if isinstance(pval, dict) and pval.get("id") == p_key:
                    target_proj = pval
                    target_project_key = pk
                    break
        
        if not target_proj and projects:
            target_project_key = list(projects.keys())[0]
            target_proj = projects[target_project_key]

        win_p = ds.get("win_path", ds_path)
        wsl_p = ds.get("wsl_path", ds_path)
        resolved_ds_path = get_platform_path(win_p, wsl_p)

        proj_base = ""
        if target_proj:
            w_base = target_proj.get("win_base", target_proj.get("paths", {}).get("windows", ""))
            l_base = target_proj.get("wsl_base", target_proj.get("paths", {}).get("wsl", ""))
            proj_base = get_platform_path(w_base, l_base)

        if resolved_ds_path and not os.path.isabs(resolved_ds_path) and proj_base:
            resolved_ds_path = os.path.abspath(os.path.join(proj_base, resolved_ds_path))

        if resolved_ds_path and os.path.isfile(resolved_ds_path):
            file_basename = os.path.basename(resolved_ds_path)
            sub_tag = "custom"
            if proj_base and os.path.exists(proj_base):
                rel_dir = os.path.relpath(os.path.dirname(resolved_ds_path), proj_base)
                sub_tag = "root" if rel_dir == "." else rel_dir.replace("\\", "/")
                if rel_dir != "." and target_proj:
                    s_subs = target_proj.setdefault("scan_subdirs", [])
                    norm_sub = rel_dir.replace("\\", "/")
                    if norm_sub not in s_subs:
                        s_subs.insert(0, norm_sub)

            if not ds_name:
                ds_name = f"{file_basename[:-5]} ({sub_tag})"

            if ds.get("is_default", False):
                ds_cfg["default_dataset"] = ds_name
                if target_proj:
                    target_proj["default_preload"] = file_basename
                if target_project_key:
                    app_settings["default_project"] = target_project_key
                    ds_cfg["default_project"] = target_project_key

            if "hidden_datasets" in ds_cfg and ds_name in ds_cfg["hidden_datasets"]:
                ds_cfg["hidden_datasets"].remove(ds_name)

            custom_list = ds_cfg.setdefault("custom_datasets", [])
            existing = [c for c in custom_list if isinstance(c, dict) and (c.get("name") == ds_name or c.get("path") == resolved_ds_path)]
            if not existing:
                custom_list.append({
                    "name": ds_name,
                    "path": resolved_ds_path,
                    "win_path": win_p if win_p else resolved_ds_path,
                    "wsl_path": wsl_p if wsl_p else resolved_ds_path,
                    "project_key": target_project_key,
                })

            imported_datasets.append({
                "name": ds_name,
                "path": resolved_ds_path,
                "project": target_project_key,
                "is_default": ds.get("is_default", False)
            })

    # Process annotation_colors if present
    anno_colors = data.get("annotation_colors")
    if not anno_colors and "dataset" in data and isinstance(data["dataset"], dict):
        anno_colors = data["dataset"].get("annotation_colors")
    if anno_colors and isinstance(anno_colors, dict):
        cfg_annos = ds_cfg.setdefault("annotation_colors", {})
        for col_k, col_v in anno_colors.items():
            if isinstance(col_v, dict):
                cfg_annos.setdefault(col_k, {}).update(col_v)
        app_settings.setdefault("annotation_colors", {}).update(cfg_annos)

    if persist:
        os.makedirs(USER_DIR, exist_ok=True)
        user_projects_file = os.path.join(USER_DIR, "projects.yaml")
        with open(user_projects_file, "w", encoding="utf-8") as f:
            yaml.safe_dump(projects, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
        save_dataset_config(ds_cfg)
        save_settings_config(app_settings)

    msg = f"Successfully imported {len(imported_datasets)} dataset(s) into project '{target_project_key}'."
    if ds_cfg.get("default_dataset"):
        msg += f" Default dataset set to '{ds_cfg['default_dataset']}'."

    return True, msg, {
        "imported_datasets": imported_datasets,
        "default_project": app_settings.get("default_project"),
        "default_dataset": ds_cfg.get("default_dataset"),
    }
