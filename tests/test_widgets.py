"""Unit tests for UI widget utilities."""
import pytest
from clairescope.ui.widgets import draggable_multiselect

def test_draggable_multiselect_filtering():
    options = ["basal", "spinous", "granular"]
    
    # When default has invalid elements, only valid elements are preserved
    res = draggable_multiselect("Test", options=options, default=["basal", "invalid_state"])
    assert "basal" in res
    assert "invalid_state" not in res

def test_draggable_multiselect_empty_fallback():
    options = ["basal", "spinous", "granular"]
    
    # When default has zero valid elements, it falls back to all options
    res = draggable_multiselect("Test", options=options, default=["completely_different_state"])
    assert res == options
