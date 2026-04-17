# plotlyMol

<p align="center">
  <img src="logo.svg" alt="plotlyMol Logo" width="300"/>
</p>

[![Tests](https://github.com/The-Schultz-Lab/plotlyMol/actions/workflows/test.yml/badge.svg)](https://github.com/The-Schultz-Lab/plotlyMol/actions/workflows/test.yml)
[![Lint](https://github.com/The-Schultz-Lab/plotlyMol/actions/workflows/lint.yml/badge.svg)](https://github.com/The-Schultz-Lab/plotlyMol/actions/workflows/lint.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Interactive molecular visualizations with Plotly. Supports SMILES, XYZ, MOL/PDB, and cube orbitals.

## Features

- 3D ball-and-stick, stick, and VDW representations
- SMILES-to-3D embedding via RDKit
- XYZ, MOL/SDF (single), and PDB input support
- Cube file orbital isosurfaces
- **Vibrational mode visualization** (Gaussian, ORCA, Molden formats)
  - Static displacement arrows
  - Animated vibrations with interactive controls
  - Heatmap coloring by displacement magnitude
- Streamlit GUI for interactive exploration

## Installation

### From source (recommended for now)

```bash
git clone https://github.com/The-Schultz-Lab/plotlyMol.git
cd plotlyMol

# Create and activate the conda environment (includes all dependencies)
conda env create -f environment.yml
conda activate plotlymol

# Install the package in editable mode
pip install -e .
```

> **Note:** [Conda](https://docs.conda.io/en/latest/) is required. If you don't have it, install [Miniforge](https://github.com/conda-forge/miniforge) (recommended) or Miniconda. All packages are installed from the `conda-forge` channel.

### Updating the environment

If `environment.yml` changes after pulling new commits:

```bash
conda env update -f environment.yml --prune
```

## Quick start

```python
from plotlymol3d import draw_3D_rep

# Draw a molecule from SMILES
fig = draw_3D_rep(smiles="CCNCOCSC", mode="ball+stick", ambient=0.1)
fig.show()

# Draw from XYZ file
fig = draw_3D_rep(xyzfile="path/to/file.xyz", mode="ball+stick", ambient=0.1)
fig.show()
```

### Orbitals from cube files

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(
	cubefile="path/to/file.cube",
	molfile="path/to/file.mol",
	mode="ball+stick",
	ambient=0.1,
	cubedraw="orbitals",
	orbital_opacity=0.25,
	orbital_colors=["darkorange", "darkblue"],
)
fig.show()
```

### Vibrational mode visualization

Visualize molecular vibrations from quantum chemistry calculations. Supports Gaussian `.log`, ORCA `.out`, and Molden `.molden` files.

**Static displacement arrows:**

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(
    smiles="O",  # Water molecule
    vibration_file="water_freq.log",
    vibration_mode=1,  # First vibrational mode
    vibration_display="arrows",
    vibration_amplitude=1.5,
)
fig.show()
```

**Animated vibration:**

```python
from plotlymol3d import parse_vibrations, create_vibration_animation
from rdkit.Chem import MolFromSmiles, AddHs
from rdkit.Chem.AllChem import EmbedMolecule

# Parse vibration data
vib_data = parse_vibrations("water_freq.log")

# Create molecule
mol = MolFromSmiles("O")
mol = AddHs(mol)
EmbedMolecule(mol)

# Generate animation
fig = create_vibration_animation(
    vib_data=vib_data,
    mode_number=1,
    mol=mol,
    amplitude=0.5,
    n_frames=20,  # Smoother with more frames
    mode="ball+stick"
)
fig.show()
```

**Heatmap coloring by displacement:**

```python
from plotlymol3d import draw_3D_rep, parse_vibrations, add_vibrations_to_figure

# Create molecular figure
fig = draw_3D_rep(smiles="O", mode="ball+stick")

# Parse vibrations and add heatmap
vib_data = parse_vibrations("water_freq.log")
fig = add_vibrations_to_figure(
    fig=fig,
    vib_data=vib_data,
    mode_number=1,
    display_type="heatmap",
    heatmap_colorscale="Reds"
)
fig.show()
```

**Available parsers:**

```python
from plotlymol3d import (
    parse_gaussian_vibrations,  # Gaussian .log files
    parse_orca_vibrations,       # ORCA .out files
    parse_molden_vibrations,     # Molden .molden files
    parse_vibrations,            # Auto-detect format
)

# Auto-detect format from file extension
vib_data = parse_vibrations("calculation.log")

# Access mode data
for mode in vib_data.modes:
    print(f"Mode {mode.mode_number}: {mode.frequency:.1f} cm⁻¹")
    if mode.ir_intensity:
        print(f"  IR Intensity: {mode.ir_intensity:.1f} km/mol")
```

## GUI

Launch the Streamlit app for interactive controls:

```bash
streamlit run examples/gui_app.py
```

## Examples

### Demo Scripts

- Demo script: `python examples/demo_visualizations.py`
- Package data includes sample XYZ/MOL/CUBE files under src/plotlymol3d/

### Jupyter Notebooks

**Vibration Visualization:**

- [Vibration Visualization Basics](examples/vibration_visualization_basics.ipynb) - Getting started with vibrational modes
- [Advanced Vibration Analysis](examples/vibration_analysis_advanced.ipynb) - Batch processing, IR spectra, publication figures

**Performance Testing:**

- [Performance Benchmarking](examples/performance_benchmarking.ipynb) - Quantitative performance analysis and optimization

Launch notebooks:

```bash
jupyter notebook examples/
```

## Performance Testing

Quantitatively measure and optimize rendering performance:

**Standalone script:**

```bash
python tests/test_performance.py
```

**Interactive notebook:**

```bash
jupyter notebook examples/performance_benchmarking.ipynb
```

**Full guide:** [Performance Testing Guide](docs/PERFORMANCE_TESTING_GUIDE.md)

**Key metrics tracked:**

- Rendering time vs molecule size
- Resolution impact (8-64)
- Memory usage profiling
- Vibration parsing speed
- Animation generation performance

Use these tools to identify bottlenecks and optimize GUI responsiveness.

## Repository layout

```
plotlyMol/
├─ src/
│  └─ plotlymol3d/        # Library package code + sample data files
├─ examples/              # Demo scripts and GUI app
├─ tests/                 # Pytest suite
├─ docs/                  # Roadmap and documentation assets
├─ pyproject.toml          # Packaging and tooling configuration
├─ requirements.txt        # Consolidated dependencies
└─ README.md
```

## Roadmap

See the current roadmap in [docs/ROADMAP.md](docs/ROADMAP.md).
