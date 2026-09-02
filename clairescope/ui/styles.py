"""UI styling, Matplotlib/Plotly typography, and CSS injection."""
import streamlit as st
import matplotlib
import plotly.io as pio
from clairescope.config import load_css_styles

def apply_global_styles():
    """Inject global CSS rules and typography configurations."""
    # Matplotlib styling
    matplotlib.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 13
    })
    
    # Plotly styling
    if "plotly_white" in pio.templates:
        pio.templates["plotly_white"].layout.font.family = "Segoe UI, Arial, sans-serif"
        pio.templates["plotly_white"].layout.font.size = 15
        pio.templates["plotly_white"].layout.title.font.size = 18
        pio.templates["plotly_white"].layout.legend.font.size = 14
        pio.templates["plotly_white"].layout.legend.title.font.size = 15
        pio.templates["plotly_white"].layout.xaxis.tickfont.size = 14
        pio.templates["plotly_white"].layout.yaxis.tickfont.size = 14
        pio.templates["plotly_white"].layout.xaxis.title.font.size = 16
        pio.templates["plotly_white"].layout.yaxis.title.font.size = 16
    pio.templates.default = "plotly_white"
    
    # Inject external stylesheet
    css_content = load_css_styles()
    if css_content:
        st.markdown(f"<style>{css_content}</style>", unsafe_allow_html=True)
