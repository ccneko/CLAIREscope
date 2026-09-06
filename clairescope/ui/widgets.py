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
    """Sortable draggable multiselect chip component with strict option filtering and fallback."""
    if default is None:
        default = []
    options_list = list(options) if options is not None else []
    default_list = list(default) if default is not None else []
    
    # Filter default to only valid options
    valid_default = [x for x in default_list if x in options_list]
    if not valid_default and options_list and default_list:
        valid_default = list(options_list)

    if _draggable_multiselect_comp is not None:
        val = _draggable_multiselect_comp(label=label, options=options_list, default=valid_default, key=key)
        if val is not None:
            valid_val = [x for x in val if x in options_list]
            if valid_val:
                return valid_val
            return valid_default
    
    val = st.multiselect(label, options=options_list, default=valid_default, key=key)
    valid_val = [x for x in val if x in options_list]
    return valid_val if valid_val else valid_default
