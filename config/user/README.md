# User Custom Configurations

Files placed in this directory will override the default configurations in `config/defaults/`.
This folder is automatically excluded from Git via `.gitignore`, ensuring your private lab dataset paths, custom signatures, and local settings are never published.

### Supported Overrides:
- `projects.yaml` : Private datasets, local file paths, and project-specific metadata.
- `signatures.yaml` : Custom gene signatures and pathway panels.
- `markers.yaml` : Custom cell-type marker gene lists.
- `pathways.yaml` : Custom biological pathway databases for ORA.
- `settings.yaml` : Custom UI typography, padding, and layout parameters.
- `style.css` : Custom CSS rules.
