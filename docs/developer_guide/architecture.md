# Architecture & Developer Guide

## System Overview
<span class="brand">CLAIREscope</span> is built as a lightweight Python application using Streamlit and Scanpy:

- **Core Framework**: Python 3.10+
- **Computational Engine**: Scanpy, SciPy, NumPy, Pandas
- **Visualization Engines**: Matplotlib, Seaborn, Plotly (WebGL)
- **Configuration Layer**: Hierarchical YAML (`config/defaults/` sanitized defaults + `config/user/` git-ignored overrides)

## Directory Structure
```
CLAIREscope/
├── app.py                  # Main Streamlit application
├── clairescope/            # Core library modules
│   ├── config.py           # Configuration loader & dataset scanner
│   ├── signatures.py       # Cutaneous biology gene signatures
│   └── stats.py            # Statistical testing functions
├── config/                 # YAML configuration definitions
│   ├── defaults/           # Public sanitized defaults
│   └── user/               # Private local overrides (git-ignored)
├── docs/                   # Documentation source files (MkDocs)
└── tests/                  # Unit test suite
```
