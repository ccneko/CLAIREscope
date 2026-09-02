"""UI styling, Matplotlib/Plotly typography, and CSS injection."""
import streamlit as st
import matplotlib
import plotly.io as pio

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
    
    # Streamlit CSS
    st.markdown("""
    <style>
        .block-container, [data-testid="block-container"], .main .block-container {
            padding-top: 3.0rem !important;
            padding-bottom: 2.0rem !important;
            padding-left: 2.2rem !important;
            padding-right: 2.2rem !important;
        }
        [data-testid="stSidebarContent"], section[data-testid="stSidebar"] > div {
            padding-top: 1.5rem !important;
        }
        .clairescope-title {
            font-size: 26pt;
            font-weight: 800;
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 2px;
        }
        div[data-testid="stImage"] img {
            max-width: 80% !important;
            margin: 0 auto !important;
            display: block !important;
        }
        .st-emotion-cache-tn0cau, [data-testid="stHorizontalBlock"] {
            gap: 0.5rem !important;
        }
    </style>
    """, unsafe_allow_html=True)
