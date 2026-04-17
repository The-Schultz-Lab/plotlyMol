#!/usr/bin/env python3
"""
Dash GUI for plotlyMol3D - Interactive 3D Molecular Viewer.

Run with:
    python examples/gui_app.py

Requirements:
    pip install plotlymol3d[gui]
"""

import sys
from pathlib import Path

try:
    from plotlymol3d.app import main
except ModuleNotFoundError:
    repo_root = Path(__file__).resolve().parents[1]
    src_path = repo_root / "src"
    sys.path.insert(0, str(src_path))
    from plotlymol3d.app import main

main()
