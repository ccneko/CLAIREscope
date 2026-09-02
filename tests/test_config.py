"""Unit tests for YAML configuration loading and validation."""
import pytest
from clairescope.config import load_projects_config, load_settings_config, load_signatures_config

def test_load_projects_config():
    projects = load_projects_config()
    assert isinstance(projects, dict)
    assert "D001_Natsuga_JEB_snRNAseq" in projects
    d001 = projects["D001_Natsuga_JEB_snRNAseq"]
    assert d001["id"] == "D001"
    assert "canonical_samples" in d001
    assert "sample_colors" in d001
    assert "default_signatures" in d001

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
