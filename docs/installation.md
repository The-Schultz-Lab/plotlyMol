# Installation

## Requirements

plotlyMol requires **Python 3.9 or higher**.

## Install from Source

Currently, plotlyMol is best installed from source:

```bash
# Clone the repository
git clone https://github.com/The-Schultz-Lab/plotlyMol.git
cd plotlyMol

# Create a virtual environment (recommended)
python -m venv .venv

# Activate the virtual environment
# On Windows:
.venv\Scripts\activate
# On macOS/Linux:
source .venv/bin/activate

# Install the package
pip install -e .
```

## Dependencies

plotlyMol automatically installs the following dependencies:

### Core Runtime Dependencies

- **plotly** (≥5.0.0) - Interactive 3D plotting
- **numpy** (≥1.20.0) - Numerical arrays and operations
- **rdkit** (≥2022.3.1) - Cheminformatics toolkit
- **kaleido** (≥0.2.1) - Static image export

### Optional Dependencies

For development and testing:

```bash
pip install -r requirements.txt
```

This includes:
- **pytest** - Testing framework
- **black** - Code formatting
- **ruff** - Fast linting
- **mypy** - Type checking
- **dash** + **dash-bootstrap-components** - GUI application

## Verify Installation

Test that plotlyMol is correctly installed:

```python
from plotlymol3d import draw_3D_rep

# Create a simple molecule
fig = draw_3D_rep(smiles="CCO")
print("plotlyMol is installed correctly!")
```

## Troubleshooting

### RDKit Installation Issues

RDKit can be challenging to install on some systems. If you encounter issues:

**Using conda (recommended for RDKit):**
```bash
conda create -n plotlymol python=3.11
conda activate plotlymol
conda install -c conda-forge rdkit
pip install plotly numpy kaleido
pip install -e .
```

**Using pip on Windows:**
- Ensure you have Visual C++ Build Tools installed
- Or use pre-built wheels from [conda-forge](https://anaconda.org/conda-forge/rdkit)

### Kaleido Issues

If static image export fails:
```bash
pip install --upgrade kaleido
```

### Import Errors

If you get `ModuleNotFoundError`:
1. Ensure the virtual environment is activated
2. Try reinstalling: `pip install -e .`
3. Check that you're importing from `plotlymol3d` (with "3d")

### Platform-Specific Notes

**Windows:**
- Use `py -m venv .venv` if `python` command not found
- May need Visual C++ Build Tools for some dependencies

**macOS (Apple Silicon):**
- Some packages may require Rosetta 2
- Consider using conda for better ARM64 support

**Linux:**
- May need to install system libraries: `sudo apt-get install build-essential`

## Next Steps

Once installed, check out the [Quick Start](quickstart.md) guide to start visualizing molecules!
