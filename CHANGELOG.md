# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.1] - 2026-04-18

### Changed
- GUI migrated from Streamlit to Dash — launch with `python examples/gui_app.py`
- `gui` optional dependency updated: `dash>=2.14.0` + `dash-bootstrap-components>=1.5.0` (replaces `streamlit`)
- `launch_app.bat` and `stop_app.bat` updated for Dash process management
- Bumped `black` minimum from `>=23.0.0` to `>=24.3.0` in `requirements.txt`

### Fixed
- Removed stale `streamlit` dependency from `requirements.txt` and `environment.yml`
- Corrected all remaining Streamlit references in docs to Dash (`about.md`, `installation.md`, `contributing.md`, `tutorials/index.md`)
- Updated version number from 0.1.0 to 0.2.0 in `docs/about.md` (including BibTeX citation)
- Fixed Python minimum version in `docs/installation.md` (3.8 → 3.9)
- Updated author contact email to `The-Schultz-Lab@users.noreply.github.com` in `docs/about.md`
- Replaced broken `docs/ROADMAP.md` and `docs/PERFORMANCE_TESTING_GUIDE.md` links in `README.md` and `docs/contributing.md`
- Updated stale Streamlit cache reference in `tests/test_performance.py`
- Fixed logo not rendering on PyPI by replacing relative `logo.svg` path with absolute raw GitHub URL

## [0.2.0] - 2026-04-04

### Added
- **Vibrational Mode Visualization** - Complete system for visualizing molecular vibrations from quantum chemistry calculations
  - Three file format parsers: Gaussian (.log), ORCA (.out), Molden (.molden) with auto-detection
  - Three visualization modes:
    - Static displacement arrows using Plotly Cone traces
    - Animated vibrations with interactive controls (play/pause, frame slider)
    - Heatmap coloring by displacement magnitude
  - New `vibrations.py` module with comprehensive dataclasses and functions
  - Streamlit "Vibration Settings" section with file upload and interactive controls
  - 21 new tests achieving ~95% coverage of vibration module
  - Test fixtures for all three file formats (water molecule examples)
- Exported vibration functions in `__init__.py` for public API access
- Symbol-to-atomic-number mapping (`symbol_to_number`) in `atomProperties.py`
- Conda environment (`environment.yml`) replacing venv-based setup
- Security CI workflow: `pip-audit` dependency scanning + CodeQL static analysis
- Dependabot configuration for automated dependency updates
- Branch protection rules (required status checks, no force-push)
- PyPI publishing checklist (`docs/PYPI_PUBLISHING.md`)

### Changed

- Aromatic bond rendering now uses ring-center geometry for correct dashed bond offset direction
- Displacement arrows auto-scaled relative to molecular size for consistent visibility
- Animation caching replaced `@st.cache_data` with session-state caching and live progress bar
- Longer dashes (75% vs 60%) for aromatic bond rendering
- Updated README installation instructions for conda workflow
- Expanded test suite from 26 to 47 tests
- Bumped `requires-python` from `>=3.8` to `>=3.9`
- Fixed repo URLs in `pyproject.toml` (now correctly point to The-Schultz-Lab org)

## [0.1.0] - 2026-01-31

### Added
- Core 3D molecular visualization with Plotly and RDKit integration.
- Input support for SMILES, XYZ, MOL/PDB, and cube files.
- Streamlit GUI for interactive visualization.
- Test suite and CI workflows.
