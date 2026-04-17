# Quick Start

Get up and running with plotlyMol in minutes!

## Basic Usage

The main function you'll use is `draw_3D_rep()`, which creates interactive 3D molecular visualizations.

### From SMILES

The simplest way to visualize a molecule:

```python
from plotlymol3d import draw_3D_rep

# Visualize ethanol
fig = draw_3D_rep(smiles="CCO", mode="ball+stick")
fig.show()
```

plotlyMol automatically:
1. Parses the SMILES string
2. Generates 3D coordinates
3. Creates the interactive visualization

### From XYZ File

If you have molecular coordinates:

```python
fig = draw_3D_rep(
    xyzfile="path/to/molecule.xyz",
    mode="ball+stick",
    ambient=0.1
)
fig.show()
```

### From MOL/SDF or PDB Files

```python
# MOL/SDF file
fig = draw_3D_rep(
    molfile="molecule.mol",
    mode="ball+stick"
)

# PDB file
fig = draw_3D_rep(
    pdbfile="protein.pdb",
    mode="stick"
)
```

## Visualization Modes

plotlyMol supports three visualization modes:

### Ball and Stick (Default)

Classic molecular representation:

```python
fig = draw_3D_rep(smiles="c1ccccc1", mode="ball+stick")
fig.show()
```

- Atoms shown as spheres at VDW radii
- Bonds shown as cylinders
- Best for general molecular visualization

### Stick Only

Simplified representation emphasizing connectivity:

```python
fig = draw_3D_rep(smiles="c1ccccc1", mode="stick")
fig.show()
```

- Small atoms
- Prominent bonds
- Better for complex molecules

### Van der Waals (VDW)

Space-filling representation:

```python
fig = draw_3D_rep(smiles="c1ccccc1", mode="vdw")
fig.show()
```

- Only atoms shown (no bonds)
- Spheres at VDW radii
- Shows molecular volume

## Customizing Appearance

### Lighting

Adjust lighting for better visualization:

```python
fig = draw_3D_rep(
    smiles="CCO",
    mode="ball+stick",
    ambient=0.2,      # Ambient light (0.0-1.0)
    diffuse=0.8,      # Diffuse reflection
    specular=0.5,     # Specular highlights
    roughness=0.5     # Surface roughness
)
fig.show()
```

**Lighting Tips:**
- Lower `ambient` (0.1-0.3) for more dramatic lighting
- Higher `diffuse` (0.7-1.0) for matte appearance
- Higher `specular` (0.5-1.0) for shiny surfaces
- Lower `roughness` (0.1-0.3) for smooth, glossy look

### Resolution

Control sphere tessellation quality:

```python
fig = draw_3D_rep(
    smiles="CCO",
    resolution=64  # Higher = smoother (default: 32)
)
fig.show()
```

**Resolution Guidelines:**
- 16: Fast preview
- 32: Default, good balance
- 64: High quality for publication

### Background Color

```python
fig = draw_3D_rep(smiles="CCO", bgcolor="black")
fig.show()

# Or after creation
from plotlymol3d import format_figure
format_figure(fig, bgcolor="white", title="Ethanol")
```

## Working with Figures

### Save as HTML

```python
fig = draw_3D_rep(smiles="CCO")
fig.write_html("molecule.html")
```

The HTML file is fully interactive and self-contained.

### Export Static Image

```python
fig = draw_3D_rep(smiles="CCO")
fig.write_image("molecule.png", width=800, height=600)
```

Requires `kaleido` package (installed by default).

### Customize Further

The returned figure is a standard Plotly `Figure` object:

```python
fig = draw_3D_rep(smiles="CCO")

# Update layout
fig.update_layout(
    title="My Molecule",
    scene=dict(
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.5)
        )
    )
)

fig.show()
```

## Orbital Visualization

Visualize molecular orbitals from Gaussian cube files:

```python
fig = draw_3D_rep(
    molfile="molecule.mol",      # Molecular structure
    cubefile="orbital.cube",      # Cube file with orbital data
    mode="ball+stick",
    cubedraw="orbitals",
    orbital_isovalue=0.02,        # Isosurface threshold
    orbital_colors=["red", "blue"],  # Positive/negative phases
    orbital_opacity=0.3           # Transparency
)
fig.show()
```

## Using the GUI

<video autoplay loop muted playsinline style="width:100%;border-radius:8px;margin-bottom:1rem">
  <source src="../assets/dash-example.webm" type="video/webm">
</video>

Launch the interactive Dash app:

```bash
python examples/gui_app.py
```

Or on Windows, double-click `launch_app.bat`.

The GUI provides:

- Molecule search by name via PubChem
- Input via SMILES string or built-in sample library
- Visualization mode selection (ball+stick, stick, VDW)
- Lighting presets
- Formula and atom/bond count display

## Complete Example

Here's a comprehensive example showing multiple features:

```python
from plotlymol3d import draw_3D_rep

# Create visualization
fig = draw_3D_rep(
    smiles="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    mode="ball+stick",
    resolution=48,
    ambient=0.15,
    diffuse=0.85,
    specular=0.5,
    roughness=0.4,
    bgcolor="white"
)

# Add title
fig.update_layout(title="Aspirin (Acetylsalicylic Acid)")

# Save
fig.write_html("aspirin.html")
fig.write_image("aspirin.png", width=1200, height=900)

# Display
fig.show()
```

## Next Steps

- [User Guide](user-guide/basic-usage.md) - Comprehensive usage guide
- [API Reference](api/index.md) - Complete function documentation
- [Tutorials](tutorials/index.md) - Step-by-step tutorials
- [Examples](examples/index.md) - Gallery of examples

## Common Recipes

### High-Quality Publication Figure

```python
fig = draw_3D_rep(
    smiles="your_molecule",
    mode="ball+stick",
    resolution=64,
    ambient=0.1,
    diffuse=0.9,
    specular=0.6,
    roughness=0.3,
    bgcolor="white"
)
fig.write_image("publication.png", width=1600, height=1200, scale=2)
```

### Quick Preview

```python
fig = draw_3D_rep(smiles="your_molecule", mode="stick", resolution=16)
fig.show()
```

### Batch Processing

```python
smiles_list = ["CCO", "CC(C)O", "CCCO"]

for i, smiles in enumerate(smiles_list):
    fig = draw_3D_rep(smiles=smiles, mode="ball+stick")
    fig.write_html(f"molecule_{i}.html")
```
