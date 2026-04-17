# About plotlyMol

## Project Overview

plotlyMol is an open-source Python package for creating interactive 3D molecular visualizations. Born from the need for modern, web-based molecular graphics in chemistry education and research, plotlyMol combines the power of Plotly's interactive plotting with RDKit's comprehensive cheminformatics capabilities.

## Authors

**Jonathan Schultz**
North Carolina Central University, Assistant Professor of Chemistry

**Benjamin Lear**
(The Pennsylvania State University, Professor of Chemistry)

## Motivation

Traditional molecular visualization tools often suffer from:
- Limited interactivity
- Platform dependencies
- Difficult installation
- Static output formats

plotlyMol addresses these issues by:
- Providing fully interactive 3D graphics in the browser
- Using pure Python with minimal dependencies
- Generating standalone HTML visualizations
- Supporting modern web standards
- Offering a simple, intuitive API

## Technology Stack

plotlyMol is built on industry-standard scientific Python libraries:

### Core Dependencies

- **[Plotly](https://plotly.com/python/)** - Interactive plotting library
  - WebGL-based 3D rendering
  - Native interactivity (zoom, rotate, pan)
  - Export to HTML and static images

- **[RDKit](https://www.rdkit.org/)** - Cheminformatics toolkit
  - SMILES parsing and validation
  - 3D coordinate generation
  - Bond perception and molecular properties
  - File format readers (MOL, SDF, PDB)

- **[NumPy](https://numpy.org/)** - Numerical computing
  - Efficient array operations
  - Vector and matrix math
  - Mesh generation algorithms

### Optional Components

- **[Streamlit](https://streamlit.io/)** - Web app framework for GUI
- **[Kaleido](https://github.com/plotly/Kaleido)** - Static image export
- **pytest** - Testing framework
- **black/ruff** - Code quality tools

## Design Philosophy

### Simplicity

Simple things should be simple:
```python
from plotlymol3d import draw_3D_rep
fig = draw_3D_rep(smiles="CCO")
fig.show()
```

### Flexibility

Complex things should be possible:
```python
fig = draw_3D_rep(
    cubefile="orbital.cube",
    molfile="molecule.mol",
    mode="stick",
    cubedraw="orbitals",
    orbital_isovalue=0.02,
    orbital_colors=["red", "blue"],
    ambient=0.1
)
```

### Quality

- Comprehensive test suite
- CI/CD pipeline
- Type hints throughout
- Detailed documentation

### Open Source

- MIT License
- Public development
- Community contributions welcome
- Transparent roadmap

## Use Cases

### Education

- **Chemistry Courses**: Interactive molecular models for lectures
- **Lab Visualizations**: 3D representations of experimental results
- **Student Projects**: Easy-to-use tool for presentations

### Research

- **Publication Figures**: High-quality molecular graphics
- **Orbital Analysis**: Visualize quantum chemistry results
- **Structure Analysis**: Interactive exploration of complex molecules
- **Data Sharing**: Standalone HTML visualizations

### Industry

- **Drug Design**: Visualize potential drug candidates
- **Materials Science**: Explore molecular structures
- **Chemical Informatics**: Integrate with analysis pipelines

## Project History

### 2024 - Genesis

Initial development by Prof. Lear focused on creating a modern molecular visualization tool for computational chemistry research and teaching.

### 2025 - Refinement

- Added comprehensive test suite
- Implemented CI/CD pipeline
- Enhanced documentation
- Improved orbital visualization

### 2026 - Documentation Phase

- Created comprehensive docs with MkDocs
- Added API reference
- Developed tutorials and examples
- Prepared for public release

## Project Status

**Current Version**: 0.1.0 (Development)

**Status**: Pre-release

### Completed Phases

- Phase 1: Project Foundation
- Phase 2: Code Quality
- Phase 3: Testing & CI/CD
- Phase 4: Documentation

### Upcoming

- Phase 5: Feature Development
- Phase 6: Advanced Features
- Phase 7: Community & Distribution

## Contributing

We welcome contributions! See [Contributing Guide](contributing.md) for details.

### Ways to Contribute

- Report bugs
- Suggest features
- Improve documentation
- Submit pull requests
- Star the repository

## Acknowledgments

### Inspiration

- **[py3Dmol](https://github.com/avirshup/py3dmol)** - Python interface to 3Dmol.js
- **[nglview](https://github.com/arose/nglview)** - Jupyter widget for molecular visualization
- **[Jmol](http://jmol.sourceforge.net/)** - Pioneer in molecular visualization

### Tools & Libraries

Thanks to the developers of:
- Plotly for exceptional interactive graphics
- RDKit for comprehensive cheminformatics
- The Scientific Python ecosystem

### Community

- North Carolina Central University, Department of Chemistry and Biochemistry
- The Pennsylvania State University, Department of Chemistry
- The open-source chemistry community
- All contributors and users

## License

plotlyMol is released under the [MIT License](https://github.com/The-Schultz-Lab/plotlyMol/blob/main/LICENSE).

```
MIT License

Copyright (c) 2026 Jonathan Schultz & Benjamin Lear

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
```

## Contact

- **GitHub Issues**: [Report issues or request features](https://github.com/The-Schultz-Lab/plotlyMol/issues)
- **GitHub Discussions**: [Ask questions or share ideas](https://github.com/The-Schultz-Lab/plotlyMol/discussions)
- **Email**: jonathanschultzNU@users.noreply.github.com

## Citation

If you use plotlyMol in your research, please cite:

```bibtex
@software{plotlymol2026,
  title = {plotlyMol: Interactive 3D Molecular Visualizations},
  author = {Schultz, Jonathan and Lear, Benjamin},
  year = {2026},
  url = {https://github.com/The-Schultz-Lab/plotlyMol},
  note = {Version 0.1.0}
}
```

---

**Made at the Schultz and Lear Labs**
