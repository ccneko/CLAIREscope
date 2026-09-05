"""Unit tests for YAML configuration loading, new project creation, and dynamic dataset scanning."""
import os
import shutil
import tempfile
import pytest
from clairescope.config import (
    load_projects_config,
    load_settings_config,
    load_signatures_config,
    load_pathways_config,
    load_markers_config,
    load_css_styles,
    get_config_file_path,
    save_user_project_config,
    get_next_new_project_name,
    scan_project_datasets,
)

def test_config_file_path_resolution():
    settings_path = get_config_file_path("settings.yaml")
    assert os.path.exists(settings_path)
    assert ("defaults" in settings_path or "user" in settings_path)

def test_load_projects_config():
    projects = load_projects_config()
    assert isinstance(projects, dict)
    assert len(projects) > 0

def test_get_next_new_project_name():
    dummy_projects = {
        "PROJ_1": {"name": "New Project 1"},
        "PROJ_2": {"name": "New Project 2"},
    }
    next_name, next_key = get_next_new_project_name(dummy_projects)
    assert next_name == "New Project 3"
    assert next_key == "PROJ_NEW_003"

def test_scan_project_datasets_with_tempdir():
    temp_dir = tempfile.mkdtemp()
    try:
        # Create mock h5ad files at root and in subfolders
        root_h5ad = os.path.join(temp_dir, "dataset_root.h5ad")
        with open(root_h5ad, "w") as f:
            f.write("mock")
            
        sub_dir = os.path.join(temp_dir, "out", "2026-09-05_test")
        os.makedirs(sub_dir, exist_ok=True)
        sub_h5ad = os.path.join(sub_dir, "dataset_sub.h5ad")
        with open(sub_h5ad, "w") as f:
            f.write("mock")

        active, all_ds = scan_project_datasets(temp_dir, scan_subdirs=["out/2026-09-05_test"])
        assert len(all_ds) == 2
        assert any("dataset_root" in k for k in all_ds.keys())
        assert any("dataset_sub" in k for k in all_ds.keys())
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def test_load_settings_config():
    settings = load_settings_config()
    assert isinstance(settings, dict)
    assert "app" in settings
    assert "ui" in settings
    assert settings["ui"]["top_padding"] == "3.0rem"
    assert settings["ui"]["column_gap"] == "0.5rem"

def test_load_signatures_config():
    signatures = load_signatures_config()
    assert isinstance(signatures, dict)
    assert "human" in signatures
    assert "mouse" in signatures
    assert "Adherens Junction Complex" in signatures["human"]

def test_load_pathways_config():
    pathways = load_pathways_config()
    assert isinstance(pathways, dict)
    assert "Hallmark: Epithelial Mesenchymal Transition" in pathways
    assert "Epidermal: Basal Stem & Hemidesmosome" in pathways

def test_load_markers_config():
    markers = load_markers_config()
    assert isinstance(markers, dict)
    assert "Human" in markers
    assert "Mouse" in markers
    assert "Basal 1 (Quiescent / Anchored)" in markers["Human"]

def test_load_css_styles():
    css = load_css_styles()
    assert len(css) > 0
    assert "block-container" in css
    assert "3.0rem" in css
