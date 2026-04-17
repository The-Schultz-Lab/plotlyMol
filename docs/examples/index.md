# Examples Gallery

Collection of molecular visualization examples using plotlyMol.

## Basic Examples

### Simple Molecules

#### Ethanol (C₂H₆O)

Ball-and-stick representation of ethanol:

```python
from plotlymol3d import draw_3D_rep

fig = draw_3D_rep(
    smiles="CCO",
    mode="ball+stick",
    ambient=0.15
)
fig.show()
```

#### Benzene (C₆H₆)

Aromatic ring with dashed bonds:

```python
fig = draw_3D_rep(
    smiles="c1ccccc1",
    mode="ball+stick",
    ambient=0.1,
    bgcolor="white"
)
fig.show()
```

### Visualization Modes

Compare different rendering modes:

=== "Ball and Stick"
    ```python
    fig = draw_3D_rep(smiles="CC(C)O", mode="ball+stick")
    fig.show()
    ```

=== "Stick Only"
    ```python
    fig = draw_3D_rep(smiles="CC(C)O", mode="stick")
    fig.show()
    ```

=== "Van der Waals"
    ```python
    fig = draw_3D_rep(smiles="CC(C)O", mode="vdw")
    fig.show()
    ```

## Drug Molecules

### Aspirin (Acetylsalicylic Acid)

```python
fig = draw_3D_rep(
    smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
    mode="ball+stick",
    resolution=48,
    ambient=0.12,
    title="Aspirin"
)
fig.write_html("aspirin.html")
```

### Caffeine (C₈H₁₀N₄O₂)

```python
fig = draw_3D_rep(
    smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    mode="ball+stick",
    ambient=0.1,
    title="Caffeine"
)
fig.show()
```

### Ibuprofen

```python
fig = draw_3D_rep(
    smiles="CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    mode="ball+stick",
    title="Ibuprofen"
)
fig.show()
```

## Amino Acids

### Glycine

```python
fig = draw_3D_rep(
    smiles="C(C(=O)O)N",
    mode="ball+stick",
    ambient=0.15,
    title="Glycine"
)
fig.show()
```

### Phenylalanine

```python
fig = draw_3D_rep(
    smiles="C1=CC=C(C=C1)CC(C(=O)O)N",
    mode="ball+stick",
    title="Phenylalanine"
)
fig.show()
```

## Complex Molecules

### Glucose (β-D-Glucose)

```python
fig = draw_3D_rep(
    smiles="C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O",
    mode="ball+stick",
    resolution=64,
    ambient=0.1,
    title="β-D-Glucose"
)
fig.show()
```

### Cholesterol

```python
fig = draw_3D_rep(
    smiles="CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
    mode="ball+stick",
    resolution=48,
    ambient=0.1,
    title="Cholesterol"
)
fig.show()
```

## From Files

### XYZ File

```python
fig = draw_3D_rep(
    xyzfile="data/benzene.xyz",
    mode="ball+stick",
    ambient=0.1
)
fig.show()
```

### MOL File

```python
fig = draw_3D_rep(
    molfile="data/molecule.mol",
    mode="stick",
    resolution=48
)
fig.show()
```

### PDB File

```python
fig = draw_3D_rep(
    pdbfile="data/protein_fragment.pdb",
    mode="stick",
    ambient=0.2
)
fig.show()
```

## Orbital Visualization

### HOMO Visualization

```python
fig = draw_3D_rep(
    molfile="benzene.mol",
    cubefile="benzene_HOMO.cube",
    mode="stick",
    cubedraw="orbitals",
    orbital_isovalue=0.02,
    orbital_colors=["red", "blue"],
    orbital_opacity=0.3,
    ambient=0.1,
    bgcolor="white",
    title="Benzene HOMO"
)
fig.show()
```

### Electron Density

```python
fig = draw_3D_rep(
    molfile="molecule.mol",
    cubefile="density.cube",
    mode="ball+stick",
    cubedraw="orbitals",
    orbital_isovalue=0.01,
    orbital_colors=["purple", "orange"],
    orbital_opacity=0.25,
    title="Electron Density"
)
fig.show()
```

## Custom Styling

### High-Quality Publication Figure

```python
fig = draw_3D_rep(
    smiles="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    mode="ball+stick",
    resolution=64,
    ambient=0.08,
    diffuse=0.92,
    specular=0.6,
    roughness=0.25,
    bgcolor="white"
)

fig.update_layout(
    title=dict(
        text="Aspirin (Acetylsalicylic Acid)",
        font=dict(size=20, family="Arial")
    ),
    scene=dict(
        camera=dict(eye=dict(x=1.3, y=1.3, z=1.3))
    )
)

fig.write_image("aspirin_publication.png", width=1600, height=1200, scale=2)
```

### Dark Theme

```python
fig = draw_3D_rep(
    smiles="c1ccccc1",
    mode="ball+stick",
    ambient=0.3,
    bgcolor="rgb(17, 17, 17)",  # Dark background
    resolution=48
)

fig.update_layout(
    title=dict(
        text="Benzene",
        font=dict(color="white")
    ),
    scene=dict(
        xaxis=dict(
            backgroundcolor="rgb(17, 17, 17)",
            gridcolor="rgb(50, 50, 50)"
        ),
        yaxis=dict(
            backgroundcolor="rgb(17, 17, 17)",
            gridcolor="rgb(50, 50, 50)"
        ),
        zaxis=dict(
            backgroundcolor="rgb(17, 17, 17)",
            gridcolor="rgb(50, 50, 50)"
        )
    )
)

fig.show()
```

## Batch Processing

### Multiple Molecules

```python
molecules = {
    "Methanol": "CO",
    "Ethanol": "CCO",
    "Propanol": "CCCO",
    "Butanol": "CCCCO"
}

for name, smiles in molecules.items():
    fig = draw_3D_rep(
        smiles=smiles,
        mode="ball+stick",
        title=name
    )
    fig.write_html(f"{name.lower()}.html")
    print(f"Created {name}.html")
```

### Conformer Gallery

```python
# Generate multiple conformers (requires RDKit)
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("CCCCCC")
mol = Chem.AddHs(mol)

# Generate conformers
AllChem.EmbedMultipleConfs(mol, numConfs=5)
AllChem.MMFFOptimizeMoleculeConfs(mol)

# Visualize each conformer
for i in range(mol.GetNumConformers()):
    fig = draw_3D_rep(
        mol=mol,
        mode="stick",
        title=f"Hexane Conformer {i+1}"
    )
    fig.write_html(f"hexane_conformer_{i+1}.html")
```

## Interactive Features

### Custom Camera Angle

```python
fig = draw_3D_rep(smiles="c1ccccc1", mode="ball+stick")

fig.update_layout(
    scene_camera=dict(
        eye=dict(x=2, y=2, z=0.5),
        center=dict(x=0, y=0, z=0),
        up=dict(x=0, y=0, z=1)
    )
)

fig.show()
```

### Add Annotations

```python
fig = draw_3D_rep(smiles="CCO", mode="ball+stick")

fig.add_annotation(
    text="Hydroxyl Group",
    x=0.5, y=0.5, z=0.5,
    showarrow=True,
    arrowhead=2
)

fig.show()
```

## Complete Example Script

See `examples/demo_visualizations.py` in the repository for a comprehensive demo script showing all features.

```bash
# Run demo script
python examples/demo_visualizations.py
```

## Download Examples

All example code is available in the [GitHub repository](https://github.com/The-Schultz-Lab/plotlyMol/tree/main/examples).

## Contributing Examples

Have a cool visualization to share? We'd love to see it!

1. Create your example
2. Add it to `examples/` directory
3. Submit a pull request

See the [Contributing Guide](../contributing.md) for details.

---

Need help recreating these examples? Check the [Quick Start](../quickstart.md) guide or [API Reference](../api/index.md).
