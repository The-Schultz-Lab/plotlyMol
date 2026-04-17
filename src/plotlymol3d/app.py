"""
Dash GUI for plotlyMol3D — Interactive 3D Molecular Visualizer.

Run with:
    python examples/gui_app.py
    or
    python src/plotlymol3d/app.py
"""

from __future__ import annotations

import os
import signal
import webbrowser
from threading import Timer
from typing import Any

import dash_bootstrap_components as dbc
import requests
from dash import Dash, Input, Output, State, dcc, html
from plotly.subplots import make_subplots

from plotlymol3d import draw_3D_mol, format_figure, format_lighting, smiles_to_rdkitmol

_HOST = "127.0.0.1"
_PORT = 8050
_VIEWER_H = "calc(100vh - 90px)"

_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_TIMEOUT = 10

_SAMPLE_MOLECULES: list[dict[str, str]] = [
    {"label": "— choose a sample —", "value": ""},
    {"label": "Water (H2O)", "value": "O"},
    {"label": "Ethanol", "value": "CCO"},
    {"label": "Benzene", "value": "c1ccccc1"},
    {"label": "Caffeine", "value": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
    {"label": "Aspirin", "value": "CC(=O)OC1=CC=CC=C1C(=O)O"},
    {"label": "Naphthalene", "value": "c1ccc2ccccc2c1"},
    {"label": "Glucose", "value": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
    {"label": "Alanine", "value": "CC(N)C(=O)O"},
    {"label": "Cyclohexane", "value": "C1CCCCC1"},
]

_MODES = [
    {"label": "Ball + Stick", "value": "ball+stick"},
    {"label": "Ball", "value": "ball"},
    {"label": "Stick", "value": "stick"},
    {"label": "Van der Waals", "value": "vdw"},
]

_LIGHTING_PRESETS: dict[str, dict] = {
    "soft": {"ambient": 0.4, "diffuse": 0.8, "specular": 0.1, "roughness": 0.8},
    "default": {"ambient": 0.0, "diffuse": 1.0, "specular": 0.0, "roughness": 1.0},
    "bright": {"ambient": 0.5, "diffuse": 0.8, "specular": 0.3, "roughness": 0.5},
    "metallic": {"ambient": 0.2, "diffuse": 0.7, "specular": 1.0, "roughness": 0.1},
    "dramatic": {"ambient": 0.0, "diffuse": 1.0, "specular": 0.6, "roughness": 0.2},
}

_LIGHTING_OPTIONS = [
    {"label": "Soft", "value": "soft"},
    {"label": "Default", "value": "default"},
    {"label": "Bright", "value": "bright"},
    {"label": "Metallic", "value": "metallic"},
    {"label": "Dramatic", "value": "dramatic"},
]


def _pubchem_name_to_smiles(name: str) -> tuple[str, int]:
    """Return (isomeric_smiles, cid) for a compound name via PubChem REST API."""
    resp = requests.get(
        f"{_PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON",
        timeout=_TIMEOUT,
    )
    if resp.status_code == 404:
        raise ValueError(
            f'"{name}" was not found in PubChem. Check spelling or try a SMILES string.'
        )
    resp.raise_for_status()
    cid = resp.json()["IdentifierList"]["CID"][0]

    for prop in ("IsomericSMILES", "CanonicalSMILES"):
        prop_resp = requests.get(
            f"{_PUBCHEM_BASE}/compound/cid/{cid}/property/{prop}/JSON",
            timeout=_TIMEOUT,
        )
        if prop_resp.ok:
            record = (
                prop_resp.json().get("PropertyTable", {}).get("Properties", [{}])[0]
            )
            # PubChem may return the value under a different key than requested
            smiles = next(
                (v for k, v in record.items() if k != "CID" and isinstance(v, str)),
                None,
            )
            if smiles:
                return smiles, cid
    raise ValueError(f'Could not retrieve a SMILES string for "{name}" (CID {cid}).')


def _empty_figure():
    import plotly.graph_objects as go

    fig = go.Figure()
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        autosize=True,
        scene={
            "bgcolor": "#f8f9fa",
            "xaxis": {"showticklabels": False, "showgrid": False, "zeroline": False},
            "yaxis": {"showticklabels": False, "showgrid": False, "zeroline": False},
            "zaxis": {"showticklabels": False, "showgrid": False, "zeroline": False},
        },
        margin={"l": 0, "r": 0, "t": 30, "b": 0},
    )
    return fig


def _mol_info_children(mol) -> list:
    atom_counts: dict[str, int] = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        atom_counts[sym] = atom_counts.get(sym, 0) + 1
    formula = "".join(
        f"{sym}{cnt if cnt > 1 else ''}" for sym, cnt in sorted(atom_counts.items())
    )
    return [
        html.Div(f"Formula: {formula}"),
        html.Div(f"Atoms: {mol.GetNumAtoms()}"),
        html.Div(f"Bonds: {mol.GetNumBonds()}"),
    ]


def _render(smiles: str, mode: str, lighting: str = "soft") -> tuple:
    """Render a SMILES string and return callback outputs for the viewer."""
    _show = {"display": "block"}
    _hide = {"display": "none"}
    _placeholder_visible = {"height": _VIEWER_H, "fontSize": "1.1rem"}

    def _err(msg: str) -> tuple[Any, ...]:
        return _empty_figure(), "", msg, _hide, _placeholder_visible

    if not smiles or not smiles.strip():
        return _err("Enter a SMILES string, search by name, or choose a sample.")
    try:
        mol = smiles_to_rdkitmol(smiles.strip())
    except Exception as exc:
        return _err(f"Invalid SMILES: {exc}")
    try:
        fig = make_subplots()
        fig = format_figure(fig)
        fig = draw_3D_mol(fig, mol, mode=mode, resolution=32)
        fig = format_lighting(
            fig, **_LIGHTING_PRESETS.get(lighting, _LIGHTING_PRESETS["soft"])
        )
        fig.update_layout(autosize=True, margin={"l": 0, "r": 0, "t": 30, "b": 0})
    except Exception as exc:
        return _err(f"Rendering error: {exc}")

    return fig, _mol_info_children(mol), "", _show, _hide


def _build_layout() -> dbc.Container:
    return dbc.Container(
        [
            # -- Header ---------------------------------------------------------
            dbc.Row(
                dbc.Col(
                    [
                        html.H4("plotlyMol3D", className="mb-0 text-primary"),
                        html.P(
                            "Interactive 3D Molecular Visualization",
                            className="text-muted mb-0 small",
                        ),
                    ],
                    className="py-2 border-bottom mb-2",
                )
            ),
            # -- Body -----------------------------------------------------------
            dbc.Row(
                [
                    # Controls panel --------------------------------------------
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    # PubChem search -------------------------
                                    html.H6("Search PubChem", className="card-title"),
                                    dbc.Label(
                                        "Molecule name",
                                        html_for="pubchem-input",
                                        className="small",
                                    ),
                                    dbc.InputGroup(
                                        [
                                            dbc.Input(
                                                id="pubchem-input",
                                                placeholder="e.g. aspirin, caffeine…",
                                                type="text",
                                                debounce=False,
                                            ),
                                            dbc.Button(
                                                "Search",
                                                id="pubchem-btn",
                                                color="info",
                                                n_clicks=0,
                                            ),
                                        ],
                                        className="mb-1",
                                    ),
                                    html.Div(
                                        id="pubchem-status",
                                        className="text-muted small mb-3",
                                    ),
                                    html.Hr(),
                                    # SMILES / sample -------------------------
                                    html.H6(
                                        "Or enter directly", className="card-title"
                                    ),
                                    dbc.Label(
                                        "SMILES string",
                                        html_for="smiles-input",
                                        className="small",
                                    ),
                                    dbc.Input(
                                        id="smiles-input",
                                        placeholder="e.g. CCO for ethanol",
                                        type="text",
                                        debounce=False,
                                        className="mb-2",
                                    ),
                                    dbc.Label(
                                        "Or choose a sample",
                                        html_for="sample-select",
                                        className="small",
                                    ),
                                    dbc.Select(
                                        id="sample-select",
                                        options=_SAMPLE_MOLECULES,
                                        value="",
                                        className="mb-3",
                                    ),
                                    html.Hr(),
                                    # Display options -------------------------
                                    html.H6("Display Options", className="card-title"),
                                    dbc.Label(
                                        "Visualization mode",
                                        html_for="mode-select",
                                        className="small",
                                    ),
                                    dbc.Select(
                                        id="mode-select",
                                        options=_MODES,
                                        value="ball+stick",
                                        className="mb-2",
                                    ),
                                    dbc.Label(
                                        "Lighting",
                                        html_for="lighting-select",
                                        className="small",
                                    ),
                                    dbc.Select(
                                        id="lighting-select",
                                        options=_LIGHTING_OPTIONS,
                                        value="soft",
                                        className="mb-3",
                                    ),
                                    dbc.Button(
                                        "Visualize",
                                        id="visualize-btn",
                                        color="primary",
                                        n_clicks=0,
                                        className="w-100 mb-3",
                                    ),
                                    html.Div(
                                        id="mol-info", className="text-muted small"
                                    ),
                                    html.Hr(),
                                    dbc.Button(
                                        "Exit",
                                        id="exit-btn",
                                        color="secondary",
                                        outline=True,
                                        n_clicks=0,
                                        className="w-100",
                                        title="Shut down the server and close the app",
                                    ),
                                    dbc.Alert(
                                        id="exit-msg",
                                        is_open=False,
                                        color="success",
                                        className="mt-2 mb-0 small text-center",
                                    ),
                                ],
                                style={"overflowY": "auto", "maxHeight": _VIEWER_H},
                            ),
                            className="h-100",
                        ),
                        md=3,
                        className="mb-2",
                    ),
                    # 3D Viewer -------------------------------------------------
                    dbc.Col(
                        [
                            html.Div(
                                "Search by name, enter a SMILES string, or choose a sample — then click Visualize.",
                                id="graph-placeholder",
                                className="text-muted d-flex align-items-center justify-content-center text-center px-4",
                                style={"height": _VIEWER_H, "fontSize": "1.1rem"},
                            ),
                            html.Div(
                                dcc.Graph(
                                    id="mol-graph",
                                    style={"height": _VIEWER_H},
                                    config={"displayModeBar": True, "scrollZoom": True},
                                    figure=_empty_figure(),
                                ),
                                id="graph-container",
                                style={"display": "none"},
                            ),
                            html.Div(
                                id="error-msg", className="text-danger mt-2 small"
                            ),
                        ],
                        md=9,
                    ),
                ],
                className="g-2",
            ),
        ],
        fluid=True,
        style={"paddingBottom": "0"},
    )


def create_app() -> Dash:
    """Create and return the configured Dash application instance."""
    app = Dash(
        __name__,
        external_stylesheets=[dbc.themes.FLATLY],
        title="plotlyMol3D Viewer",
    )
    app.layout = _build_layout()

    @app.callback(
        Output("smiles-input", "value"),
        Input("sample-select", "value"),
        prevent_initial_call=True,
    )
    def fill_smiles_from_sample(value: str) -> str:
        return value or ""

    @app.callback(
        Output("mol-graph", "figure"),
        Output("mol-info", "children"),
        Output("error-msg", "children"),
        Output("graph-container", "style"),
        Output("graph-placeholder", "style"),
        Input("visualize-btn", "n_clicks"),
        State("smiles-input", "value"),
        State("mode-select", "value"),
        State("lighting-select", "value"),
        prevent_initial_call=True,
    )
    def update_figure(n_clicks: int, smiles: str, mode: str, lighting: str):
        return _render(smiles, mode, lighting)

    @app.callback(
        Output("smiles-input", "value", allow_duplicate=True),
        Output("pubchem-status", "children"),
        Output("pubchem-status", "className"),
        Output("mol-graph", "figure", allow_duplicate=True),
        Output("mol-info", "children", allow_duplicate=True),
        Output("error-msg", "children", allow_duplicate=True),
        Output("graph-container", "style", allow_duplicate=True),
        Output("graph-placeholder", "style", allow_duplicate=True),
        Input("pubchem-btn", "n_clicks"),
        State("pubchem-input", "value"),
        State("mode-select", "value"),
        State("lighting-select", "value"),
        prevent_initial_call=True,
    )
    def search_pubchem(n_clicks: int, name: str, mode: str, lighting: str):
        _placeholder_visible = {"height": _VIEWER_H, "fontSize": "1.1rem"}
        _hide = {"display": "none"}

        def _no_render(smiles_val: str, status: str, cls: str) -> tuple[Any, ...]:
            return (
                smiles_val,
                status,
                cls,
                _empty_figure(),
                "",
                "",
                _hide,
                _placeholder_visible,
            )

        if not name or not name.strip():
            return _no_render(
                "", "Enter a molecule name to search.", "text-warning small mb-3"
            )

        try:
            smiles, cid = _pubchem_name_to_smiles(name.strip())
        except ValueError as exc:
            return _no_render("", str(exc), "text-danger small mb-3")
        except requests.RequestException:
            msg = (
                "Could not reach PubChem. Check your internet connection and try again."
            )
            return _no_render("", msg, "text-danger small mb-3")

        fig, info, err, container_style, placeholder_style = _render(
            smiles, mode, lighting
        )
        status = f"Found: {name.strip().title()} (CID {cid})"
        return (
            smiles,
            status,
            "text-success small mb-3",
            fig,
            info,
            err,
            container_style,
            placeholder_style,
        )

    @app.callback(
        Output("exit-msg", "children"),
        Output("exit-msg", "is_open"),
        Input("exit-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def shutdown(n_clicks: int):
        if n_clicks:
            Timer(0.5, lambda: os.kill(os.getpid(), signal.SIGTERM)).start()
            return "Server stopped. You may close this tab.", True
        return "", False

    return app


def main() -> None:
    """Start the Dash server and open the browser."""
    app = create_app()
    Timer(1.5, lambda: webbrowser.open(f"http://{_HOST}:{_PORT}")).start()
    print(f"\n  plotlyMol3D is running at http://{_HOST}:{_PORT}")
    print("  Press Ctrl+C to stop.\n")
    app.run(debug=False, host=_HOST, port=_PORT)


if __name__ == "__main__":
    main()
