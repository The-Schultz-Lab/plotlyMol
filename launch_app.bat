@echo off
setlocal
set "ROOT=%~dp0"
cd /d "%ROOT%"

set "CONDA_PY=C:\Users\schul\miniconda3\envs\plotlymol\python.exe"
if exist "%CONDA_PY%" (
    "%CONDA_PY%" examples\gui_app.py
) else (
    python examples\gui_app.py
)
