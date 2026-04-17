#!/usr/bin/env bash
# launch_app.command — macOS launcher for plotlyMol3D
#
# Double-click this file in Finder to launch the app, or run:
#   bash launch_app.command
#
# To make it double-clickable, run once in Terminal:
#   chmod +x launch_app.command

# Change to the directory containing this script (works when double-clicked)
cd "$(dirname "$0")"

# Use the virtual environment's Python if it exists, otherwise fall back to
# the system/conda Python on $PATH
VENV_PY=".venv/bin/python"
if [ -f "$VENV_PY" ]; then
    "$VENV_PY" examples/gui_app.py
else
    python3 examples/gui_app.py
fi
