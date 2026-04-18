# Contributing to plotlyMol

Thank you for your interest in contributing to plotlyMol! This guide will help you get started.

## Ways to Contribute

### Report Bugs

Found a bug? Please [open an issue](https://github.com/The-Schultz-Lab/plotlyMol/issues/new) with:

- Clear description of the problem
- Steps to reproduce
- Expected vs actual behavior
- System information (OS, Python version)
- Minimal code example

### Suggest Features

Have an idea? We'd love to hear it! [Open an issue](https://github.com/The-Schultz-Lab/plotlyMol/issues/new) describing:

- The feature and its use case
- Why it would be valuable
- Possible implementation approach (optional)

### Improve Documentation

Documentation improvements are always welcome:

- Fix typos or clarify explanations
- Add examples or tutorials
- Improve API documentation
- Translate documentation (future)

### Submit Code

Ready to code? Great! Follow the guidelines below.

## Development Setup

### 1. Fork and Clone

```bash
# Fork on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/plotlyMol.git
cd plotlyMol

# Add upstream remote
git remote add upstream https://github.com/The-Schultz-Lab/plotlyMol.git
```

### 2. Create Virtual Environment

```bash
# Create environment
python -m venv .venv

# Activate
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
# Install in editable mode with all dependencies
pip install -r requirements.txt
pip install -e .

# Install pre-commit hooks
pre-commit install
```

### 4. Create Branch

```bash
# Create feature branch
git checkout -b feature/my-new-feature

# Or bugfix branch
git checkout -b fix/bug-description
```

## Development Workflow

### Code Style

plotlyMol follows strict code quality standards:

#### Formatting

- **Black** (line length: 88)
- **Ruff** for linting
- **isort** for import sorting

```bash
# Format code
black src/plotlymol3d tests

# Check linting
ruff check src/plotlymol3d tests

# Auto-fix linting issues
ruff check --fix src/plotlymol3d tests
```

#### Type Hints

All functions should have type hints:

```python
from typing import Optional, List
import plotly.graph_objects as go

def my_function(
    param1: str,
    param2: Optional[int] = None
) -> go.Figure:
    """Function docstring."""
    ...
```

#### Docstrings

Use Google-style docstrings:

```python
def draw_molecule(smiles: str, mode: str = "ball+stick") -> go.Figure:
    """Create 3D molecular visualization.

    Args:
        smiles: SMILES string for the molecule.
        mode: Visualization mode ("ball+stick", "stick", "vdw").

    Returns:
        Plotly Figure object with 3D visualization.

    Raises:
        ValueError: If SMILES string is invalid.

    Example:
        >>> fig = draw_molecule("CCO", mode="ball+stick")
        >>> fig.show()
    """
    ...
```

### Testing

#### Writing Tests

- Place tests in `tests/` directory
- Use descriptive test names: `test_<function>_<scenario>`
- Include docstrings explaining what's being tested
- Use fixtures from `conftest.py` for sample data

```python
def test_draw_3D_rep_from_smiles():
    """Test basic visualization from SMILES string."""
    fig = draw_3D_rep(smiles="CCO", mode="ball+stick")

    assert fig is not None
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0  # Has traces
```

#### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=plotlymol3d --cov-report=term-missing

# Run specific test file
pytest tests/test_visualization.py

# Run specific test
pytest tests/test_visualization.py::test_draw_3D_rep_from_smiles

# Run with verbose output
pytest -v
```

#### Coverage Requirements

- Aim for >60% overall coverage
- New features should have >80% coverage
- Bug fixes should include regression tests

### Pre-commit Hooks

Hooks run automatically before each commit:

- Trailing whitespace removal
- End-of-file fixer
- YAML syntax check
- Black formatting
- Ruff linting

```bash
# Install hooks
pre-commit install

# Run manually
pre-commit run --all-files
```

### Continuous Integration

All pull requests must pass CI checks:

#### Test Workflow

- Runs on Ubuntu, macOS, Windows
- Tests Python 3.9, 3.10, 3.11, 3.12
- Runs full test suite with coverage

#### Lint Workflow

- Black formatting check
- Ruff linting
- Mypy type checking

## Pull Request Process

### 1. Prepare Your Changes

```bash
# Sync with upstream
git fetch upstream
git rebase upstream/main

# Run tests
pytest

# Check formatting
black --check src/plotlymol3d tests
ruff check src/plotlymol3d tests
mypy src/plotlymol3d
```

### 2. Commit Your Changes

Follow conventional commit format:

```bash
# Feature
git commit -m "Add support for MOL2 file format"

# Bug fix
git commit -m "Fix bond order detection for aromatic rings"

# Documentation
git commit -m "Update API documentation for draw_3D_rep"

# Tests
git commit -m "Add tests for XYZ file parsing"
```

### 3. Push and Create PR

```bash
# Push to your fork
git push origin feature/my-new-feature

# Create PR on GitHub
# Use the PR template
# Link related issues
```

### 4. PR Requirements

Your PR should:

- [ ] Have a clear, descriptive title
- [ ] Include purpose and changes in description
- [ ] Reference related issues (#123)
- [ ] Pass all CI checks
- [ ] Include tests for new features
- [ ] Update documentation if needed
- [ ] Follow code style guidelines

### 5. Review Process

- Maintainers will review your PR
- Address feedback and requested changes
- Keep the PR focused (one feature/fix)
- Be patient and respectful

## Code Organization

### Adding New Features

#### 1. Input Format Support

To add support for a new file format:

1. Add parser function in `plotlyMol3D.py`:

   ```python
   def newformat_to_rdkitmol(filepath: str) -> Chem.Mol:
       """Parse new format file to RDKit Mol."""
       ...
   ```

2. Update `draw_3D_rep()` to accept new parameter

3. Add tests in `tests/test_input_processing.py`

4. Update documentation

#### 2. Visualization Features

To add new visualization features:

1. Add function in appropriate module
2. Add corresponding tests
3. Update API documentation
4. Add example to quickstart guide

#### 3. GUI Features

To enhance the Dash GUI:

1. Modify `src/plotlymol3d/app.py`
2. Test manually with `python examples/gui_app.py`
3. Update user guide

## Documentation

### Building Docs Locally

```bash
# Install doc dependencies (already in requirements.txt)
pip install mkdocs mkdocs-material mkdocstrings[python]

# Serve locally (auto-reload)
mkdocs serve

# Build static site
mkdocs build
```

Visit `http://127.0.0.1:8000` to view docs.

### Adding Documentation

- **New pages**: Add `.md` files in `docs/`
- **Update nav**: Edit `mkdocs.yml` navigation
- **API docs**: Use mkdocstrings syntax in API pages
- **Examples**: Add to `docs/examples/` or `examples/`

## Community Guidelines

### Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Provide constructive feedback
- Focus on the code, not the person

### Communication

- **GitHub Issues**: Bug reports, feature requests
- **GitHub Discussions**: Questions, ideas, general discussion
- **Pull Requests**: Code contributions

### Getting Help

Stuck? Need guidance?

- Check existing issues and PRs
- Read the documentation
- Ask in GitHub Discussions
- Tag maintainers if urgent

## Release Process

(For maintainers)

### Version Numbering

We use [Semantic Versioning](https://semver.org/):

- **Major** (1.0.0): Breaking changes
- **Minor** (0.1.0): New features, backwards compatible
- **Patch** (0.0.1): Bug fixes

### Release Checklist

1. Update version in `pyproject.toml`
2. Update `CHANGELOG.md`
3. Run full test suite
4. Build documentation
5. Create release on GitHub
6. Build and upload to PyPI (future)

## Development Roadmap

See [GitHub Issues](https://github.com/The-Schultz-Lab/plotlyMol/issues) for planned features and improvements.

## Recognition

Contributors are recognized in:

- `CHANGELOG.md` for each release
- GitHub contributors page
- Special thanks in documentation

## Questions?

Have questions about contributing?

- Open a [Discussion](https://github.com/The-Schultz-Lab/plotlyMol/discussions)
- Check existing documentation
- Reach out to maintainers

Thank you for contributing to plotlyMol!
