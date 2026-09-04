"""Unit tests for YAML configuration loading and user override hierarchy."""
import os
import pytest
from clairescope.config import (
    load_projects_config,
    load_settings_config,
    load_signatures_config,
    load_pathways_config,
    load_markers_config,
    load_css_styles,
    get_config_file_path,
)

def test_config_file_path_resolution():
    # settings.yaml exists in defaults/
    settings_path = get_config_file_path("settings.yaml")
    assert os.path.exists(settings_path)
    assert ("defaults" in settings_path or "user" in settings_path)

def test_load_projects_config():
    projects = load_projects_config()
    assert isinstance(projects, dict)
    assert len(projects) > 0

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
