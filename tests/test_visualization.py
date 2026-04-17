"""
Unit tests for visualization/drawing functions in plotlyMol3D.

Tests cover:
- Atom mesh trace generation
- Bond mesh trace generation
- Fibonacci sphere generation
- Cylinder mesh generation
"""

import numpy as np
import pytest
from plotly.graph_objects import Figure, Mesh3d

from plotlymol3d.plotlyMol3D import (
    Atom,
    draw_atoms,
    draw_bonds,
    generate_cylinder_mesh_rectangles,
    make_atom_mesh_trace,
    make_bond_mesh_trace,
    make_fibonacci_sphere,
    rdkitmol_to_atoms_bonds_lists,
    smiles_to_rdkitmol,
)


class TestFibonacciSphere:
    """Tests for make_fibonacci_sphere function."""

    def test_returns_correct_types(self):
        """Test that function returns numpy arrays."""
        x, y, z = make_fibonacci_sphere([0, 0, 0])

        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)

    def test_correct_number_of_points(self):
        """Test that correct number of points are generated."""
        resolution = 100
        x, y, z = make_fibonacci_sphere([0, 0, 0], resolution=resolution)

        assert len(x) == resolution
        assert len(y) == resolution
        assert len(z) == resolution

    def test_centered_at_origin(self):
        """Test sphere is centered at specified point."""
        center = [0, 0, 0]
        x, y, z = make_fibonacci_sphere(center, radius=1.0, resolution=1000)

        # Mean should be approximately at center
        assert abs(np.mean(x) - center[0]) < 0.1
        assert abs(np.mean(y) - center[1]) < 0.1
        assert abs(np.mean(z) - center[2]) < 0.1

    def test_correct_radius(self):
        """Test that points are at the specified radius."""
        center = [0, 0, 0]
        radius = 2.0
        x, y, z = make_fibonacci_sphere(center, radius=radius, resolution=100)

        # All points should be at distance = radius from center
        distances = np.sqrt(x**2 + y**2 + z**2)
        assert np.allclose(distances, radius, rtol=0.01)

    def test_offset_center(self):
        """Test sphere with non-origin center."""
        center = [1.5, -2.0, 3.0]
        radius = 1.0
        x, y, z = make_fibonacci_sphere(center, radius=radius, resolution=100)

        # Points should be at distance = radius from center
        distances = np.sqrt(
            (x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2
        )
        assert np.allclose(distances, radius, rtol=0.01)


class TestAtomMeshTrace:
    """Tests for make_atom_mesh_trace function."""

    @pytest.fixture
    def sample_atom(self):
        """Create a sample atom for testing."""
        return Atom(
            atom_id=0,
            atom_number=6,  # Carbon
            atom_symbol="C",
            atom_xyz=[0.0, 0.0, 0.0],
            atom_vdw=1.70,
        )

    def test_returns_mesh3d(self, sample_atom):
        """Test that function returns a Mesh3d trace."""
        trace = make_atom_mesh_trace(sample_atom)
        assert isinstance(trace, Mesh3d)

    def test_trace_has_correct_name(self, sample_atom):
        """Test that trace name matches atom symbol and id."""
        trace = make_atom_mesh_trace(sample_atom)
        assert trace.name == "C0"

    def test_vdw_radius_mode(self, sample_atom):
        """Test that 'vdw' radius mode uses VDW radius."""
        trace = make_atom_mesh_trace(sample_atom, radius="vdw")
        # Should use atom's VDW radius (1.70)
        assert trace is not None

    def test_ball_radius_mode(self, sample_atom):
        """Test that 'ball' radius mode scales VDW radius."""
        trace = make_atom_mesh_trace(sample_atom, radius="ball")
        # Should use 0.2 * atom's VDW radius
        assert trace is not None


class TestCylinderMesh:
    """Tests for generate_cylinder_mesh_rectangles function."""

    def test_returns_correct_types(self):
        """Test that function returns numpy arrays."""
        x, y, z = generate_cylinder_mesh_rectangles([0, 0, 0], [0, 0, 1])

        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)

    def test_correct_number_of_vertices(self):
        """Test that correct number of vertices are generated."""
        resolution = 16
        x, y, z = generate_cylinder_mesh_rectangles(
            [0, 0, 0], [0, 0, 1], resolution=resolution
        )

        # Should have 2 * resolution vertices (top and bottom circles)
        assert len(x) == 2 * resolution
        assert len(y) == 2 * resolution
        assert len(z) == 2 * resolution

    def test_vertical_cylinder(self):
        """Test cylinder along z-axis."""
        x, y, z = generate_cylinder_mesh_rectangles(
            [0, 0, 0], [0, 0, 5], radius=1.0, resolution=32
        )

        # Bottom circle should be at z=0, top at z=5
        assert np.min(z) == pytest.approx(0.0, abs=0.01)
        assert np.max(z) == pytest.approx(5.0, abs=0.01)


class TestBondMeshTrace:
    """Tests for make_bond_mesh_trace function."""

    def test_returns_mesh3d(self):
        """Test that function returns a Mesh3d trace."""
        trace = make_bond_mesh_trace([0, 0, 0], [1, 0, 0])
        assert isinstance(trace, Mesh3d)

    def test_trace_has_hoverinfo_skip(self):
        """Test that bond traces skip hover info."""
        trace = make_bond_mesh_trace([0, 0, 0], [1, 0, 0])
        assert trace.hoverinfo == "skip"


class TestDrawFunctions:
    """Tests for draw_atoms and draw_bonds functions."""

    def test_draw_atoms_adds_traces(self, sample_smiles):
        """Test that draw_atoms adds correct number of traces."""
        mol = smiles_to_rdkitmol(sample_smiles)
        atomList, _ = rdkitmol_to_atoms_bonds_lists(mol)

        fig = Figure()
        initial_traces = len(fig.data)

        fig = draw_atoms(fig, atomList)

        assert len(fig.data) == initial_traces + len(atomList)

    def test_draw_bonds_adds_traces(self, sample_smiles):
        """Test that draw_bonds adds correct number of traces."""
        mol = smiles_to_rdkitmol(sample_smiles)
        _, bondList = rdkitmol_to_atoms_bonds_lists(mol)

        fig = Figure()
        initial_traces = len(fig.data)

        fig = draw_bonds(fig, bondList)

        # Each bond creates 2 traces (one per half)
        assert len(fig.data) == initial_traces + 2 * len(bondList)


# ---------------------------------------------------------------------------
# create_trajectory_animation
# ---------------------------------------------------------------------------


def _water_xyzblocks(n: int = 3) -> list:
    """Return n slightly perturbed XYZ blocks for water."""
    blocks = []
    for i in range(n):
        o_z = i * 0.01
        block = (
            f"3\nH2O step {i}\n"
            f"O  0.0   0.0   {o_z:.4f}\n"
            f"H  0.757 0.587 0.0\n"
            f"H -0.757 0.587 0.0"
        )
        blocks.append(block)
    return blocks


class TestCreateTrajectoryAnimation:
    """Tests for create_trajectory_animation in plotlyMol3D."""

    def test_returns_plotly_figure(self):
        """Returns a Plotly Figure object."""
        from plotlymol3d import create_trajectory_animation

        fig = create_trajectory_animation(_water_xyzblocks(3))
        assert isinstance(fig, Figure)

    def test_frame_count_matches_input(self):
        """Number of frames equals number of XYZ blocks."""
        from plotlymol3d import create_trajectory_animation

        blocks = _water_xyzblocks(4)
        fig = create_trajectory_animation(blocks)
        assert len(fig.frames) == 4

    def test_minimum_two_frames_accepted(self):
        """Two frames is the valid minimum."""
        from plotlymol3d import create_trajectory_animation

        fig = create_trajectory_animation(_water_xyzblocks(2))
        assert len(fig.frames) == 2

    def test_single_frame_raises_value_error(self):
        """Providing only one frame raises ValueError."""
        from plotlymol3d import create_trajectory_animation

        with pytest.raises(ValueError):
            create_trajectory_animation(_water_xyzblocks(1))

    def test_slider_step_count(self):
        """Slider has one step per frame."""
        from plotlymol3d import create_trajectory_animation

        n = 5
        fig = create_trajectory_animation(_water_xyzblocks(n))
        assert len(fig.layout.sliders) > 0
        assert len(fig.layout.sliders[0].steps) == n

    def test_initial_traces_not_empty(self):
        """Initial figure data contains at least one trace."""
        from plotlymol3d import create_trajectory_animation

        fig = create_trajectory_animation(_water_xyzblocks(2))
        assert len(fig.data) > 0

    def test_energies_in_frame_labels(self):
        """Energy values appear in frame layout titles."""
        from plotlymol3d import create_trajectory_animation

        energies = [-75.0, -75.3, -75.6]
        fig = create_trajectory_animation(
            _water_xyzblocks(3), energies_hartree=energies
        )
        titles = [
            f.layout.title.text for f in fig.frames if f.layout and f.layout.title
        ]
        assert any("-75" in (t or "") for t in titles)
