# plotting.py
from __future__ import annotations
from typing import Dict, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt

# Core model utilities
from model.model_core import (
    simulate,
    compute_outputs_along_solution,
    drivers_from_fn_or_const,  # <-- adapter ensures array-based Drivers
)

# Display names & units
NAMES_UNITS = {
    # Environment
    "L":  dict(name="Light (L)",               unit="µmol photons m⁻² s⁻¹"),
    "T":  dict(name="Temperature (T)",         unit="°C"),
    "Nu": dict(name="DIN (Nu)",                unit="mol N L⁻¹"),
    "X":  dict(name="Prey (X)",                unit="mol C L⁻¹"),
    # Outputs (Fig. 5 second row and on)
    "jCP": dict(name="Photosynthesis rate (jCP)",           unit="mol C mol C⁻¹ d⁻¹"),
    "rho_C": dict(name="C shared by symbiont (ρC)",         unit="mol C mol C⁻¹ d⁻¹"),
    "rho_N": dict(name="N shared by host (ρN)",             unit="mol N mol C⁻¹ d⁻¹"),
    "jeL": dict(name="Excess light (jeL)",                  unit="mol photons mol C⁻¹ d⁻¹"),
    "jNPQ": dict(name="Non-photochem. quench. (jNPQ)",      unit="mol photons mol C⁻¹ d⁻¹"),
    "cROS_increase": dict(name="ROS increase (cROS − 1)",   unit="–"),
    "expulsion": dict(name="Expulsion flux b(cROS−1)·jST,0",unit="mol C mol C⁻¹ d⁻¹"),
    "jSG": dict(name="Symbiont growth (jSG)",               unit="mol C mol C⁻¹ d⁻¹"),
    "jHG": dict(name="Host growth (jHG)",                   unit="mol C mol C⁻¹ d⁻¹"),
    "net_S": dict(name="Net symbiont growth Ṡ/S",           unit="d⁻¹"),
    "net_H": dict(name="Net host growth Ḣ/H",               unit="d⁻¹"),
    "S_over_H": dict(name="Symbiont/host ratio S/H",        unit="–"),
}

# Second-row and onward
panel_KEYS = [
    "S_over_H", "jCP", "rho_C", "rho_N",
    "jNPQ", "cROS_increase",
    "net_S", "net_H",
]

def plot_panel(
    params: Dict[str, float],
    drivers,                      # old-style (constants/callables) or new Drivers
    show: bool = True,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    First row: L, T, Nu, X (time series).
    Remaining rows: model outputs in `panel_KEYS`.
    Time runs 0..params['tmax'].
    """
    # Run simulation (simulate() already adapts drivers internally)
    sol = simulate(params=params, drivers=drivers)

    # Compute model output series (also adapts drivers internally)
    series = compute_outputs_along_solution(sol, params, drivers)
    t = series["t"]

    # Adapt drivers *once* to arrays aligned with t for environment plotting
    env = drivers_from_fn_or_const(drivers, params, t)
    L_series, T_series, Nu_series, X_series = env.L, env.T, env.Nu, env.X

    # Grid: 4 columns; 1 env row + enough rows for the panels
    ncols = 4
    rows_for_panels = (len(panel_KEYS) + ncols - 1) // ncols
    nrows = 1 + max(1, rows_for_panels)

    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 3.0 * nrows), sharex=True)
    axes = np.array(axes).ravel()

    # Row 1: environment
    env_keys = [("L", L_series), ("T", T_series), ("Nu", Nu_series), ("X", X_series)]
    for i, (k, arr) in enumerate(env_keys):
        ax = axes[i]
        ax.plot(t, arr)
        ax.set_title(NAMES_UNITS[k]["name"])
        ax.set_ylabel(NAMES_UNITS[k]["unit"])
        ax.grid(True, alpha=0.3)

    # Remaining: model outputs
    start = ncols  # after environment row
    for i, key in enumerate(panel_KEYS):
        ax = axes[start + i]
        ax.plot(t, series[key])
        meta = NAMES_UNITS[key]
        ax.set_title(meta["name"])
        ax.set_ylabel(meta["unit"])
        ax.grid(True, alpha=0.3)

    # Hide unused axes and label x on bottom row
    total_axes = nrows * ncols
    used_axes = start + len(panel_KEYS)
    for j in range(used_axes, total_axes):
        fig.delaxes(axes[j])
    for ax in axes[-ncols:]:
        if ax.has_data():
            ax.set_xlabel("Time (days)")

    fig.tight_layout()
    if show:
        plt.show()
    return fig, axes
