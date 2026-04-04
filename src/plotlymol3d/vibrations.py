"""Molecular vibration visualization module.

This module provides functionality for parsing and visualizing molecular vibrations
from quantum chemistry calculations. Supports Gaussian, ORCA, and Molden file formats.

Classes:
    VibrationalMode: Data structure for a single vibrational mode
    VibrationalData: Container for complete vibrational analysis data

Functions:
    parse_gaussian_vibrations: Parse Gaussian .log/.out files
    parse_orca_vibrations: Parse ORCA .out files
    parse_molden_vibrations: Parse Molden .molden files
    parse_vibrations: Auto-detect format and parse
    create_displacement_arrows: Generate arrow traces for displacements
    create_vibration_animation: Create animated vibration
    create_heatmap_colored_figure: Apply heatmap coloring by displacement
    add_vibrations_to_figure: Main integration function
"""

import re
from dataclasses import dataclass
from typing import Optional

import numpy as np
import plotly.graph_objects as go


@dataclass
class VibrationalMode:
    """Single vibrational mode data.

    Attributes:
        mode_number: Mode index (1-based, following QM convention)
        frequency: Frequency in cm⁻¹ (negative for imaginary)
        ir_intensity: IR intensity (km/mol), None if not available
        displacement_vectors: Array of shape (n_atoms, 3) with Cartesian displacements
        is_imaginary: True if frequency is negative (imaginary frequency)
    """

    mode_number: int
    frequency: float
    ir_intensity: Optional[float]
    displacement_vectors: np.ndarray  # shape: (n_atoms, 3)
    is_imaginary: bool


@dataclass
class VibrationalData:
    """Complete vibrational analysis data.

    Attributes:
        coordinates: Atomic coordinates (n_atoms, 3) in Angstroms
        atomic_numbers: List of atomic numbers
        modes: List of VibrationalMode objects
        source_file: Original filename for reference
        program: Source program ("gaussian", "orca", "molden")
    """

    coordinates: np.ndarray  # shape: (n_atoms, 3)
    atomic_numbers: list[int]
    modes: list[VibrationalMode]
    source_file: str
    program: str

    def get_mode(self, mode_number: int) -> Optional[VibrationalMode]:
        """Retrieve mode by number.

        Args:
            mode_number: Mode number to retrieve (1-based)

        Returns:
            VibrationalMode object if found, None otherwise
        """
        for mode in self.modes:
            if mode.mode_number == mode_number:
                return mode
        return None

    def get_displacement_magnitudes(self, mode_number: int) -> np.ndarray:
        """Calculate displacement magnitude for each atom in a mode.

        Args:
            mode_number: Mode number to analyze (1-based)

        Returns:
            Array of displacement magnitudes for each atom
        """
        mode = self.get_mode(mode_number)
        if mode is None:
            return np.array([])
        return np.linalg.norm(mode.displacement_vectors, axis=1)


def parse_gaussian_vibrations(filepath: str) -> VibrationalData:
    """Parse vibrational data from Gaussian .log/.out file.

    Strategy:
    1. Find final "Standard orientation" section for coordinates
    2. Locate "Harmonic frequencies" section
    3. Parse frequency blocks (printed in groups of 3-5 modes)
    4. Extract frequencies, IR intensities, displacement vectors

    Args:
        filepath: Path to Gaussian output file

    Returns:
        VibrationalData object

    Raises:
        ValueError: If no vibrations found or file malformed
    """
    with open(filepath) as f:
        content = f.read()

    # 1. Extract coordinates from last "Standard orientation"
    # Pattern: Standard orientation section with coordinates between two dash lines
    coord_pattern = r"Standard orientation:.*?-{50,}.*?-{50,}\s*(.*?)\s*-{50,}"
    coord_matches = list(re.finditer(coord_pattern, content, re.DOTALL))
    if not coord_matches:
        raise ValueError("No coordinates found in Gaussian file")

    last_coords_block = coord_matches[-1].group(1)
    coords = []
    atomic_numbers = []

    # Parse coordinate lines: "  1   6   0   0.000000   0.000000   0.000000"
    # Use \s* at start to handle first line which may have stripped whitespace
    coord_line_pattern = r"\s*\d+\s+(\d+)\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    for match in re.finditer(coord_line_pattern, last_coords_block):
        atomic_num = int(match.group(1))
        x, y, z = map(float, match.groups()[1:])
        atomic_numbers.append(atomic_num)
        coords.append([x, y, z])

    coords = np.array(coords)
    n_atoms = len(coords)

    # 2. Find vibration section
    vib_pattern = (
        r"Harmonic frequencies \(cm\*\*-1\), IR intensities.*?(?=\n\s*-{3,}|\Z)"
    )
    vib_match = re.search(vib_pattern, content, re.DOTALL)
    if not vib_match:
        raise ValueError("No vibrational data found in Gaussian file")

    vib_section = vib_match.group(0)

    # 3. Parse frequency blocks (Gaussian prints 3 modes per block typically)
    modes = []

    # Find all mode number header lines
    # Pattern: "     1         2         3" (at start of line, followed by newline or symmetry labels)
    # Must NOT match atom lines like "     1   6     0.00..."
    # Mode headers have larger spacing between numbers and are followed by symmetry labels or newline
    mode_header_pattern = (
        r"^\s+(\d+)\s{5,}(\d+)(?:\s{5,}(\d+))?(?:\s{5,}(\d+))?(?:\s{5,}(\d+))?\s*$"
    )

    # Find all blocks by splitting at mode headers
    mode_blocks = []
    for match in re.finditer(mode_header_pattern, vib_section, re.MULTILINE):
        match.start()
        mode_nums = [int(x) for x in match.groups() if x is not None]

        # Find the next mode header or end of section
        next_match = re.search(
            mode_header_pattern, vib_section[match.end() :], re.MULTILINE
        )
        end_pos = match.end() + next_match.start() if next_match else len(vib_section)

        block_content = vib_section[match.end() : end_pos]
        mode_blocks.append((mode_nums, block_content))

    # Process each block
    for mode_nums, block in mode_blocks:

        # Extract frequencies
        freq_match = re.search(
            r"Frequencies --\s+([-\d.]+)\s+([-\d.]+)\s*(?:([-\d.]+))?\s*(?:([-\d.]+))?\s*(?:([-\d.]+))?",
            block,
        )
        if not freq_match:
            continue

        frequencies = [float(x) for x in freq_match.groups() if x is not None]

        # Extract IR intensities
        ir_match = re.search(
            r"IR Inten\s+--\s+([-\d.]+)\s+([-\d.]+)\s*(?:([-\d.]+))?\s*(?:([-\d.]+))?\s*(?:([-\d.]+))?",
            block,
        )
        ir_intensities: list[Optional[float]] = (
            [float(x) for x in ir_match.groups() if x is not None]
            if ir_match
            else [None] * len(frequencies)
        )

        # Extract displacement vectors
        disp_start = block.find("Atom  AN")
        if disp_start == -1:
            continue

        disp_section = block[disp_start:]
        disp_lines = disp_section.split("\n")[1:]  # Skip header

        # Initialize displacement arrays for this block's modes
        n_modes_in_block = len(frequencies)
        displacements = [np.zeros((n_atoms, 3)) for _ in range(n_modes_in_block)]

        # Parse displacement lines
        # Each line has one atom with X,Y,Z for each mode in the block
        # Format: "  1   6    X1  Y1  Z1    X2  Y2  Z2    X3  Y3  Z3  ..."
        atom_idx = 0
        line_idx = 0

        while atom_idx < n_atoms and line_idx < len(disp_lines):
            line = disp_lines[line_idx]
            if not line.strip():
                line_idx += 1
                continue

            parts = line.split()
            if len(parts) < 5:  # Need at least: atom_num, AN, X, Y, Z
                line_idx += 1
                continue

            # Skip atom number and atomic number (first 2 fields)
            displacement_values = parts[2:]

            # Parse X, Y, Z for each mode in the block
            # Each mode gets 3 consecutive values (X, Y, Z)
            for mode_idx in range(n_modes_in_block):
                start_idx = mode_idx * 3
                if start_idx + 2 < len(displacement_values):
                    displacements[mode_idx][atom_idx, 0] = float(
                        displacement_values[start_idx]
                    )
                    displacements[mode_idx][atom_idx, 1] = float(
                        displacement_values[start_idx + 1]
                    )
                    displacements[mode_idx][atom_idx, 2] = float(
                        displacement_values[start_idx + 2]
                    )

            atom_idx += 1
            line_idx += 1

        # Create VibrationalMode objects
        for _mode_idx, (mode_num, freq, ir_int, disp) in enumerate(
            zip(
                mode_nums[: len(frequencies)],
                frequencies,
                ir_intensities,
                displacements,
            )
        ):
            modes.append(
                VibrationalMode(
                    mode_number=mode_num,
                    frequency=freq,
                    ir_intensity=ir_int,
                    displacement_vectors=disp,
                    is_imaginary=(freq < 0),
                )
            )

    return VibrationalData(
        coordinates=coords,
        atomic_numbers=atomic_numbers,
        modes=modes,
        source_file=filepath,
        program="gaussian",
    )


def parse_vibrations(filepath: str) -> VibrationalData:
    """Auto-detect format and parse vibrational data.

    Detection strategy:
    - .log/.out with "Gaussian" in content → Gaussian
    - .out with "O   R   C   A" in content → ORCA
    - .molden or "[Molden Format]" in content → Molden

    Args:
        filepath: Path to vibration file

    Returns:
        VibrationalData object

    Raises:
        ValueError: If format cannot be detected or parsing fails
    """
    import os

    ext = os.path.splitext(filepath)[1].lower()

    with open(filepath) as f:
        first_kb = f.read(1024)

    if ext == ".molden" or "[Molden Format]" in first_kb:
        return parse_molden_vibrations(filepath)
    elif "Gaussian" in first_kb[:500]:
        return parse_gaussian_vibrations(filepath)
    elif "O   R   C   A" in first_kb[:500]:
        return parse_orca_vibrations(filepath)
    else:
        raise ValueError(f"Cannot detect vibration file format: {filepath}")


# Placeholder functions for Phase 3
def parse_orca_vibrations(filepath: str) -> VibrationalData:
    """Parse vibrational data from ORCA .out file.

    Strategy:
    1. Find "CARTESIAN COORDINATES (ANGSTROEM)" for geometry
    2. Locate "VIBRATIONAL FREQUENCIES" section
    3. Parse frequency list (skip first 6 translations/rotations)
    4. Extract displacement vectors from "NORMAL MODES" section
    5. ORCA prints modes in columns (6 modes per block)

    Args:
        filepath: Path to ORCA output file

    Returns:
        VibrationalData object

    Raises:
        ValueError: If no vibrations found or file malformed
    """
    with open(filepath) as f:
        content = f.read()

    # 1. Extract coordinates from "CARTESIAN COORDINATES (ANGSTROEM)"
    coord_pattern = r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*-+\s*(.*?)(?:\n\s*\n|\Z)"
    coord_match = re.search(coord_pattern, content, re.DOTALL)
    if not coord_match:
        raise ValueError("No coordinates found in ORCA file")

    coord_section = coord_match.group(1)
    coords = []
    atomic_numbers = []

    # Parse coordinate lines: "  O      0.000000    0.000000    0.119262"
    from .atomProperties import symbol_to_number

    coord_line_pattern = r"\s*([A-Z][a-z]?)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    for match in re.finditer(coord_line_pattern, coord_section):
        symbol = match.group(1)
        x, y, z = map(float, match.groups()[1:])
        atomic_num = symbol_to_number.get(symbol)
        if atomic_num is None:
            raise ValueError(f"Unknown element symbol: {symbol}")
        atomic_numbers.append(atomic_num)
        coords.append([x, y, z])

    coords = np.array(coords)
    n_atoms = len(coords)

    # 2. Find vibrational frequencies section
    freq_pattern = (
        r"VIBRATIONAL FREQUENCIES\s*-+.*?Scaling factor.*?\n(.*?)(?:\n\s*-+|\Z)"
    )
    freq_match = re.search(freq_pattern, content, re.DOTALL)
    if not freq_match:
        raise ValueError("No vibrational frequencies found in ORCA file")

    freq_section = freq_match.group(1)

    # Parse frequency lines: "   6:      1595.03 cm**-1"
    frequencies = []
    freq_line_pattern = r"\s*(\d+):\s+([-\d.]+)\s+cm"
    for match in re.finditer(freq_line_pattern, freq_section):
        mode_idx = int(match.group(1))
        freq = float(match.group(2))
        frequencies.append((mode_idx, freq))

    # Filter out first 6 modes (translations/rotations with freq ≈ 0)
    vibrational_freqs = [(idx, freq) for idx, freq in frequencies if abs(freq) > 1.0]

    # 3. Find normal modes section
    modes_pattern = r"NORMAL MODES\s*-+(.*?)(?:\n-{3,}|\Z)"
    modes_match = re.search(modes_pattern, content, re.DOTALL)
    if not modes_match:
        raise ValueError("No normal modes found in ORCA file")

    modes_section = modes_match.group(1)

    # Parse mode blocks (ORCA prints 6 modes per block)
    # Format:
    #                   6          7          8
    #  0   0.000000  0.000000  0.000000
    #  1   0.070000 -0.000000 -0.060000
    # ...

    # Find all mode blocks by header lines (integers only, no floats)
    mode_header_pattern = r"^\s+(\d+)(?:\s+(\d+))?(?:\s+(\d+))?(?:\s+(\d+))?(?:\s+(\d+))?(?:\s+(\d+))?\s*$"

    mode_blocks = []
    lines = modes_section.split("\n")
    i = 0

    while i < len(lines):
        line = lines[i]
        header_match = re.match(mode_header_pattern, line)

        if header_match:
            # Extract mode numbers
            mode_nums = [int(x) for x in header_match.groups() if x is not None]

            # Read displacement data
            displacement_data = []
            i += 1

            # Read next n_atoms * 3 lines (x, y, z for each atom)
            for _ in range(n_atoms * 3):
                if i >= len(lines):
                    break
                data_line = lines[i]

                # Parse line: " 0   0.000000  0.000000  0.000000"
                parts = data_line.split()
                if len(parts) >= 2:  # Row index + at least one value
                    values = [float(x) for x in parts[1:]]
                    displacement_data.append(values)
                i += 1

            if len(displacement_data) == n_atoms * 3:
                mode_blocks.append((mode_nums, displacement_data))
        else:
            i += 1

    # 4. Build VibrationalMode objects
    modes: list[VibrationalMode] = []

    # Organize displacement data by mode
    all_mode_displacements: dict[int, list] = {}

    for mode_nums, disp_data in mode_blocks:
        # disp_data is (n_atoms * 3) rows × len(mode_nums) columns
        for mode_idx, mode_num in enumerate(mode_nums):
            if mode_num not in all_mode_displacements:
                all_mode_displacements[mode_num] = []

            # Extract column for this mode
            for row in disp_data:
                if mode_idx < len(row):
                    all_mode_displacements[mode_num].append(row[mode_idx])

    # Create VibrationalMode objects for vibrational modes only
    for mode_idx, freq in vibrational_freqs:
        if mode_idx in all_mode_displacements:
            disp_vector = all_mode_displacements[mode_idx]

            # Reshape from flat list to (n_atoms, 3)
            if len(disp_vector) == n_atoms * 3:
                disp_array = np.array(disp_vector).reshape((n_atoms, 3))

                # Renumber modes starting from 1
                mode_number = len(modes) + 1

                modes.append(
                    VibrationalMode(
                        mode_number=mode_number,
                        frequency=freq,
                        ir_intensity=None,  # ORCA doesn't always include IR in std output
                        displacement_vectors=disp_array,
                        is_imaginary=(freq < 0),
                    )
                )

    return VibrationalData(
        coordinates=coords,
        atomic_numbers=atomic_numbers,
        modes=modes,
        source_file=filepath,
        program="orca",
    )


def parse_molden_vibrations(filepath: str) -> VibrationalData:
    """Parse vibrational data from Molden format file.

    Strategy:
    1. Parse [Atoms] section for coordinates and atomic numbers
    2. Parse [FREQ] section for frequencies
    3. Parse [INT] section for IR intensities (optional)
    4. Parse [FR-NORM-COORD] section for displacement vectors

    Molden format is well-structured with clear section markers.
    Coordinates can be in Angstroms ([Atoms] Angs) or Bohr ([Atoms] AU).

    Args:
        filepath: Path to .molden file

    Returns:
        VibrationalData object

    Raises:
        ValueError: If required sections are missing or malformed
    """
    with open(filepath) as f:
        content = f.read()

    # 1. Parse [Atoms] section
    # Format: [Atoms] Angs (or AU)
    #         O  1   8   0.000000  0.000000  0.119262
    atoms_pattern = r"\[Atoms\]\s+(Angs|AU)\s*\n(.*?)(?:\n\[|\Z)"
    atoms_match = re.search(atoms_pattern, content, re.DOTALL | re.IGNORECASE)
    if not atoms_match:
        raise ValueError("No [Atoms] section found in Molden file")

    unit = atoms_match.group(1).upper()
    atoms_section = atoms_match.group(2)

    coords = []
    atomic_numbers = []

    # Parse atom lines: "O  1   8   0.000000  0.000000  0.119262"
    atom_line_pattern = (
        r"\s*([A-Z][a-z]?)\s+\d+\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    )
    for match in re.finditer(atom_line_pattern, atoms_section):
        match.group(1)
        atomic_num = int(match.group(2))
        x, y, z = map(float, match.groups()[2:])

        # Convert from Bohr to Angstroms if needed
        if unit == "AU":
            bohr_to_angstrom = 0.529177
            x *= bohr_to_angstrom
            y *= bohr_to_angstrom
            z *= bohr_to_angstrom

        atomic_numbers.append(atomic_num)
        coords.append([x, y, z])

    coords = np.array(coords)
    n_atoms = len(coords)

    # 2. Parse [FREQ] section
    freq_pattern = r"\[FREQ\]\s*\n(.*?)(?:\n\[|\Z)"
    freq_match = re.search(freq_pattern, content, re.DOTALL)
    if not freq_match:
        raise ValueError("No [FREQ] section found in Molden file")

    freq_section = freq_match.group(1)
    frequencies = []

    for line in freq_section.strip().split("\n"):
        line = line.strip()
        if line:
            frequencies.append(float(line))

    # 3. Parse [INT] section (optional)
    int_pattern = r"\[INT\]\s*\n(.*?)(?:\n\[|\Z)"
    int_match = re.search(int_pattern, content, re.DOTALL)

    ir_intensities: list[Optional[float]] = []
    if int_match:
        int_section = int_match.group(1)
        for line in int_section.strip().split("\n"):
            line = line.strip()
            if line:
                ir_intensities.append(float(line))

    # If no IR intensities, use None for all modes
    if not ir_intensities:
        ir_intensities = [None] * len(frequencies)
    elif len(ir_intensities) < len(frequencies):
        # Pad with None if not enough intensities
        ir_intensities.extend([None] * (len(frequencies) - len(ir_intensities)))

    # 4. Parse [FR-NORM-COORD] section
    norm_coord_pattern = r"\[FR-NORM-COORD\]\s*\n(.*?)(?:\n\[|\Z)"
    norm_coord_match = re.search(norm_coord_pattern, content, re.DOTALL)
    if not norm_coord_match:
        raise ValueError("No [FR-NORM-COORD] section found in Molden file")

    norm_coord_section = norm_coord_match.group(1)

    # Parse vibration blocks
    # Format:
    # vibration 1
    # 0.000000 0.070000 0.000000
    # 0.000000 -0.560000 -0.430000
    # 0.000000 -0.560000 0.430000
    modes = []

    vibration_blocks = re.split(r"vibration\s+(\d+)", norm_coord_section)

    # vibration_blocks = ['', '1', '\n0.000...', '2', '\n0.000...', ...]
    for i in range(1, len(vibration_blocks), 2):
        if i + 1 >= len(vibration_blocks):
            break

        mode_number = int(vibration_blocks[i])
        block_content = vibration_blocks[i + 1].strip()

        # Parse displacement vectors (one line per atom, 3 values per line)
        displacement_vectors = []
        for line in block_content.split("\n"):
            line = line.strip()
            if line:
                values = list(map(float, line.split()))
                if len(values) == 3:
                    displacement_vectors.append(values)

        if len(displacement_vectors) == n_atoms:
            disp_array = np.array(displacement_vectors)

            # Get corresponding frequency and IR intensity
            freq_idx = mode_number - 1
            if freq_idx < len(frequencies):
                freq = frequencies[freq_idx]
                ir_int = (
                    ir_intensities[freq_idx] if freq_idx < len(ir_intensities) else None
                )

                modes.append(
                    VibrationalMode(
                        mode_number=mode_number,
                        frequency=freq,
                        ir_intensity=ir_int,
                        displacement_vectors=disp_array,
                        is_imaginary=(freq < 0),
                    )
                )

    return VibrationalData(
        coordinates=coords,
        atomic_numbers=atomic_numbers,
        modes=modes,
        source_file=filepath,
        program="molden",
    )


# Phase 2: Visualization Functions
def create_displacement_arrows(
    vib_data: VibrationalData,
    mode_number: int,
    amplitude: float = 1.0,
    arrow_scale: float = 1.0,
    color: str = "red",
    show_small_displacements: bool = False,
    displacement_threshold: float = 0.01,
) -> list[go.Cone]:
    """Create 3D arrow traces for vibrational displacements.

    Uses Plotly Cone traces for vector field visualization.

    Args:
        vib_data: VibrationalData object
        mode_number: Which mode to visualize (1-based)
        amplitude: Displacement amplitude multiplier
        arrow_scale: Visual scale for arrow size
        color: Arrow color
        show_small_displacements: If False, hide arrows below threshold
        displacement_threshold: Minimum displacement magnitude to show

    Returns:
        List of Plotly Cone traces

    Raises:
        ValueError: If mode not found
    """
    mode = vib_data.get_mode(mode_number)
    if mode is None:
        raise ValueError(f"Mode {mode_number} not found")

    coords = vib_data.coordinates

    # Calculate molecular size for auto-scaling
    # Use the span of coordinates as a reference scale
    coord_ranges = np.ptp(coords, axis=0)  # peak-to-peak (max - min) for each axis
    molecular_size = np.mean(coord_ranges)

    # Scale displacements relative to molecular size
    # This ensures arrows are visible but not overwhelming
    displacements = mode.displacement_vectors * amplitude

    # Calculate magnitudes for filtering
    magnitudes = np.linalg.norm(displacements, axis=1)

    # Auto-scale: normalize to reasonable fraction of molecular size
    max_magnitude = np.max(magnitudes) if len(magnitudes) > 0 else 1.0
    if max_magnitude > 0:
        # Scale so max arrow is ~15% of molecular size by default
        auto_scale = (0.15 * molecular_size) / max_magnitude
        displacements = displacements * auto_scale

    # Filter by threshold if requested
    if not show_small_displacements:
        mask = magnitudes > displacement_threshold
        coords_filtered = coords[mask]
        displacements_filtered = displacements[mask]
    else:
        coords_filtered = coords
        displacements_filtered = displacements

    if len(coords_filtered) == 0:
        return []

    # Create cone trace (single trace for all arrows)
    # sizeref controls cone size - smaller = bigger cones, use small value
    trace = go.Cone(
        x=coords_filtered[:, 0],
        y=coords_filtered[:, 1],
        z=coords_filtered[:, 2],
        u=displacements_filtered[:, 0],
        v=displacements_filtered[:, 1],
        w=displacements_filtered[:, 2],
        colorscale=[[0, color], [1, color]],
        sizemode="scaled",
        sizeref=arrow_scale * 0.3,  # Scale down: smaller sizeref = larger cones
        showscale=False,
        name=f"Mode {mode_number} ({mode.frequency:.1f} cm⁻¹)",
        hovertemplate=(
            f"<b>Mode {mode_number}</b><br>"
            f"Frequency: {mode.frequency:.2f} cm⁻¹<br>"
            "Displacement: %{u:.3f}, %{v:.3f}, %{w:.3f}<br>"
            "<extra></extra>"
        ),
    )

    return [trace]


def create_vibration_animation(
    vib_data: VibrationalData,
    mode_number: int,
    mol,
    amplitude: float = 0.5,
    n_frames: int = 20,
    mode: str = "ball+stick",
    resolution: int = 16,
    progress_callback=None,
) -> go.Figure:
    """Create animated vibration using Plotly frames.

    Generates sinusoidal motion along the vibrational mode displacement vectors:
    coords(t) = coords_equilibrium + amplitude · sin(2πt) · displacement

    Strategy:
    1. Get displacement vectors for the mode
    2. Generate coordinates for each frame using sinusoidal motion
    3. For each frame, create a temporary RDKit Mol with updated coords
    4. Generate molecular visualization for each frame
    5. Combine into Plotly animation with frames

    Args:
        vib_data: VibrationalData object
        mode_number: Which mode to animate
        mol: RDKit Mol object (for bond connectivity)
        amplitude: Vibration amplitude in Angstroms
        n_frames: Number of animation frames (more = smoother)
        mode: Visualization mode ("ball+stick", "stick")
        resolution: Sphere resolution for rendering (lower = faster)
        progress_callback: Optional callable(current, total) called after each frame

    Returns:
        Plotly Figure with animation frames

    Raises:
        ValueError: If mode not found
    """
    from rdkit import Chem

    mode_obj = vib_data.get_mode(mode_number)
    if mode_obj is None:
        raise ValueError(f"Mode {mode_number} not found")

    equilibrium_coords = vib_data.coordinates
    displacement_vectors = mode_obj.displacement_vectors

    # Import draw_3D_mol from the main module
    from .plotlyMol3D import draw_3D_mol

    # Generate frames
    frames = []
    frame_data = []

    for i in range(n_frames):
        # Calculate phase (0 to 2π)
        t = i / n_frames
        phase = 2 * np.pi * t

        # Calculate displaced coordinates
        # coords(t) = coords_eq + A·sin(phase)·displacement
        displaced_coords = (
            equilibrium_coords + amplitude * np.sin(phase) * displacement_vectors
        )

        # Create a temporary mol with updated coordinates
        mol_copy = Chem.Mol(mol)
        conf = mol_copy.GetConformer()

        for atom_idx in range(mol_copy.GetNumAtoms()):
            x, y, z = displaced_coords[atom_idx]
            conf.SetAtomPosition(atom_idx, (float(x), float(y), float(z)))

        # Generate molecular visualization for this frame
        # draw_3D_mol requires an empty figure as first argument
        empty_fig = go.Figure()
        fig_frame = draw_3D_mol(empty_fig, mol_copy, mode=mode, resolution=resolution)

        # Extract trace data for this frame
        frame_traces = list(fig_frame.data)

        # Create frame
        frame = go.Frame(
            data=frame_traces,
            name=f"frame_{i}",
            layout=go.Layout(
                title_text=f"Mode {mode_number}: {mode_obj.frequency:.1f} cm⁻¹ (Frame {i+1}/{n_frames})"
            ),
        )
        frames.append(frame)

        # Store data for initial display (frame 0)
        if i == 0:
            frame_data = frame_traces

        # Call progress callback if provided
        if progress_callback is not None:
            progress_callback(i + 1, n_frames)

    # Create figure with initial frame
    fig = go.Figure(data=frame_data, frames=frames)

    # Add animation controls
    fig.update_layout(
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "buttons": [
                    {
                        "label": "▶ Play/Loop",
                        "method": "animate",
                        "args": [
                            None,
                            {
                                "frame": {"duration": 50, "redraw": True},
                                "fromcurrent": False,  # Always start from beginning
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                    {
                        "label": "Pause",
                        "method": "animate",
                        "args": [
                            [None],
                            {
                                "frame": {"duration": 0, "redraw": False},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                    },
                ],
                "x": 0.1,
                "y": 0.0,
                "xanchor": "left",
                "yanchor": "bottom",
            }
        ],
        sliders=[
            {
                "active": 0,
                "steps": [
                    {
                        "args": [
                            [f"frame_{k}"],
                            {
                                "frame": {"duration": 0, "redraw": True},
                                "mode": "immediate",
                                "transition": {"duration": 0},
                            },
                        ],
                        "label": str(k + 1),
                        "method": "animate",
                    }
                    for k in range(n_frames)
                ],
                "x": 0.1,
                "len": 0.85,
                "xanchor": "left",
                "y": 0.0,
                "yanchor": "top",
                "pad": {"b": 10, "t": 50},
                "currentvalue": {
                    "visible": True,
                    "prefix": "Frame: ",
                    "xanchor": "right",
                    "font": {"size": 14},
                },
                "transition": {"duration": 0},
            }
        ],
        title=f"Vibrational Mode {mode_number}: {mode_obj.frequency:.1f} cm⁻¹",
        scene={
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "aspectmode": "data",
        },
    )

    return fig


def create_heatmap_colored_figure(
    fig: go.Figure,
    vib_data: VibrationalData,
    mode_number: int,
    colorscale: str = "Reds",
    show_colorbar: bool = True,
) -> go.Figure:
    """Color atoms by displacement magnitude in vibrational mode.

    Modifies existing atom (sphere) traces in the figure to use colors based
    on displacement magnitude. Atoms with larger displacements get "hotter" colors.

    Strategy:
    1. Get displacement magnitudes for each atom in the mode
    2. Normalize magnitudes to 0-1 range
    3. Find all Mesh3d traces in the figure (atoms and bonds)
    4. Match traces to atoms by comparing coordinates
    5. Apply colorscale to atom traces based on displacement
    6. Optionally add colorbar

    Args:
        fig: Existing molecular figure (from draw_3D_mol)
        vib_data: VibrationalData object
        mode_number: Which mode to visualize
        colorscale: Plotly colorscale name (e.g., "Reds", "Blues", "Viridis")
        show_colorbar: Whether to show color scale bar

    Returns:
        Modified figure with heatmap coloring

    Raises:
        ValueError: If mode not found
    """
    mode = vib_data.get_mode(mode_number)
    if mode is None:
        raise ValueError(f"Mode {mode_number} not found")

    # Get displacement magnitudes for each atom
    magnitudes = vib_data.get_displacement_magnitudes(mode_number)

    # Normalize to 0-1 range
    max_magnitude = np.max(magnitudes)
    if max_magnitude > 0:
        normalized_magnitudes = magnitudes / max_magnitude
    else:
        normalized_magnitudes = magnitudes

    # Match atoms to traces by comparing coordinates
    # Assumption: atom traces are Mesh3d traces where the centroid matches atom position
    coords = vib_data.coordinates

    for trace_idx, trace in enumerate(fig.data):
        if trace.type == "mesh3d" and trace.x is not None and len(trace.x) > 0:
            # Calculate centroid of the mesh
            centroid = np.array([np.mean(trace.x), np.mean(trace.y), np.mean(trace.z)])

            # Find closest atom
            distances = np.linalg.norm(coords - centroid, axis=1)
            closest_atom_idx = np.argmin(distances)

            # If distance is small (< 1.0 Å), this is likely an atom trace
            # Note: Threshold is generous to handle coordinate differences between
            # molecule generation methods (SMILES vs QM coords)
            if distances[closest_atom_idx] < 1.0:
                # Apply color based on displacement magnitude
                magnitude = normalized_magnitudes[closest_atom_idx]

                # Update trace with colorscale
                # Note: Plotly Mesh3d expects intensity values for colorscale
                # We need to set all vertices to the same intensity value
                n_vertices = len(trace.x)
                intensities = np.full(n_vertices, magnitude)

                # Clear existing color attributes to avoid conflicts
                fig.data[trace_idx].vertexcolor = None
                fig.data[trace_idx].facecolor = None

                # Set intensity-based coloring
                fig.data[trace_idx].intensity = intensities
                fig.data[trace_idx].colorscale = colorscale
                fig.data[trace_idx].showscale = (
                    show_colorbar and trace_idx == 0
                )  # Only first trace shows colorbar

                if show_colorbar and trace_idx == 0:
                    # Add colorbar configuration
                    fig.data[trace_idx].colorbar = {
                        "title": {"text": "Displacement<br>Magnitude"},
                        "tickmode": "linear",
                        "tick0": 0,
                        "dtick": 0.2,
                        "thickness": 15,
                        "len": 0.7,
                    }

    return fig


def add_vibrations_to_figure(
    fig: go.Figure,
    vib_data: VibrationalData,
    mode_number: int,
    display_type: str = "arrows",
    amplitude: float = 1.0,
    arrow_scale: float = 1.0,
    arrow_color: str = "red",
    heatmap_colorscale: str = "Reds",
    show_colorbar: bool = True,
) -> go.Figure:
    """Add vibrational mode visualization to existing molecular figure.

    Main entry point following draw_cube_orbitals pattern.

    Args:
        fig: Existing molecular figure
        vib_data: VibrationalData object
        mode_number: Which mode to visualize
        display_type: "arrows", "heatmap", or "both"
        amplitude: Displacement amplitude for arrows
        arrow_scale: Visual scale for arrows
        arrow_color: Color for displacement arrows
        heatmap_colorscale: Colorscale for heatmap mode
        show_colorbar: Show colorbar in heatmap mode

    Returns:
        Modified figure with vibration overlay

    Raises:
        ValueError: If mode not found
    """
    if display_type in ("arrows", "both"):
        arrow_traces = create_displacement_arrows(
            vib_data=vib_data,
            mode_number=mode_number,
            amplitude=amplitude,
            arrow_scale=arrow_scale,
            color=arrow_color,
        )
        for trace in arrow_traces:
            fig.add_trace(trace)

    if display_type in ("heatmap", "both"):
        fig = create_heatmap_colored_figure(
            fig=fig,
            vib_data=vib_data,
            mode_number=mode_number,
            colorscale=heatmap_colorscale,
            show_colorbar=show_colorbar,
        )

    # Update title
    mode = vib_data.get_mode(mode_number)
    if mode:
        freq_label = f"{mode.frequency:.1f} cm⁻¹"
        if mode.is_imaginary:
            freq_label = f"{freq_label} (imaginary)"

        current_title = fig.layout.title.text if fig.layout.title else ""
        fig.update_layout(title=f"{current_title}<br>Mode {mode_number}: {freq_label}")

    # Fix aspect ratio to prevent arrows from squishing/stretching the molecule view
    # Use "cube" mode to maintain equal axis scaling
    # Update only aspectmode/aspectratio without replacing entire scene
    if display_type in ("arrows", "both"):
        fig.update_layout(
            scene_aspectmode="cube",
            scene_aspectratio={"x": 1, "y": 1, "z": 1},
        )

    return fig
