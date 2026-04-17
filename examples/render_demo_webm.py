"""
Generate a rotating-molecule demo animation and save as docs/assets/demo.webm.

Usage:
    python examples/render_demo_webm.py
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from plotly.subplots import make_subplots  # noqa: E402

from plotlymol3d import (  # noqa: E402
    draw_3D_mol,
    format_figure,
    format_lighting,
    smiles_to_rdkitmol,
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # caffeine
MODE = "ball+stick"
RESOLUTION = 64

N_FRAMES = 120  # one full 360° rotation (smoother at 30fps = 4s)
FPS = 30
CAMERA_DISTANCE = 1.8  # eye distance from origin
CAMERA_ELEVATION = 0.4  # z component of eye (slight upward tilt)
WIDTH = 1920
HEIGHT = 1080

OUT_PATH = Path(__file__).resolve().parents[1] / "docs" / "assets" / "demo.webm"
FFMPEG = r"C:\Users\schul\ffmpeg-7.1.1-full_build\bin\ffmpeg.exe"

# ---------------------------------------------------------------------------
# Build base figure
# ---------------------------------------------------------------------------
print("Building molecule figure...")
mol = smiles_to_rdkitmol(SMILES)
fig = make_subplots()
fig = format_figure(fig)
fig = draw_3D_mol(fig, mol, mode=MODE, resolution=RESOLUTION)
fig = format_lighting(fig)
fig.update_layout(
    width=WIDTH,
    height=HEIGHT,
    margin={"l": 0, "r": 0, "t": 0, "b": 0},
    paper_bgcolor="white",
    scene={
        "bgcolor": "white",
        "xaxis": {
            "showticklabels": False,
            "showgrid": False,
            "zeroline": False,
            "visible": False,
        },
        "yaxis": {
            "showticklabels": False,
            "showgrid": False,
            "zeroline": False,
            "visible": False,
        },
        "zaxis": {
            "showticklabels": False,
            "showgrid": False,
            "zeroline": False,
            "visible": False,
        },
    },
)

# ---------------------------------------------------------------------------
# Render frames
# ---------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as tmpdir:
    print(f"Rendering {N_FRAMES} frames...")
    for i in range(N_FRAMES):
        angle = 2 * math.pi * i / N_FRAMES
        eye_x = CAMERA_DISTANCE * math.cos(angle)
        eye_y = CAMERA_DISTANCE * math.sin(angle)
        eye_z = CAMERA_ELEVATION
        fig.update_layout(scene_camera={"eye": {"x": eye_x, "y": eye_y, "z": eye_z}})
        frame_path = os.path.join(tmpdir, f"frame_{i:04d}.png")
        fig.write_image(frame_path, format="png")
        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{N_FRAMES}")

    # -----------------------------------------------------------------------
    # Encode with ffmpeg
    # -----------------------------------------------------------------------
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    print(f"Encoding {OUT_PATH} ...")
    cmd = [
        FFMPEG,
        "-y",
        "-framerate",
        str(FPS),
        "-i",
        os.path.join(tmpdir, "frame_%04d.png"),
        "-c:v",
        "libvpx-vp9",
        "-b:v",
        "0",
        "-crf",
        "24",
        "-pix_fmt",
        "yuva420p",  # supports transparency if background changes
        str(OUT_PATH),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ffmpeg error:\n", result.stderr)
        sys.exit(1)

print(f"Done -> {OUT_PATH}")
