"""Custom interactive Streamlit UI widgets."""
import os
import streamlit as st
import streamlit.components.v1 as components

APP_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COMPONENT_DIR = os.path.join(APP_DIR, "components", "draggable_multiselect")

if os.path.exists(COMPONENT_DIR):
    _draggable_multiselect_comp = components.declare_component("draggable_multiselect", path=COMPONENT_DIR)
else:
    _draggable_multiselect_comp = None

def draggable_multiselect(label: str, options: list, default: list = None, key: str = None) -> list:
    """Sortable draggable multiselect chip component."""
    if default is None:
        default = []
    if _draggable_multiselect_comp is not None:
        val = _draggable_multiselect_comp(label=label, options=options, default=default, key=key)
        if val is not None:
            return val
    return st.multiselect(label, options=options, default=default, key=key)
