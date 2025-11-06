# fest_plot_panel.py
"""
Quick plot runner for the Option A model (changing environment).
Uses your parameters.py seasonal drivers (AnnualCycle) so L/T/Nu/X vary with time.
"""

from parameters import default_params, drivers_from_params  # ← your real params & seasonal drivers
from plotting import plot_panel


def main():
    # --- parameters and environment (seasonal/annual cycles) ---
    p = default_params()

    # OPTIONAL: make sure there's visible variation in this demo.
    # If your defaults have zero amplitudes, uncomment and tweak:
    # p["L_amp"] = 300.0
    # p["T_amp"] = 4.0
    # p["Nu_amp"] = 1e-7
    # p["X_amp"]  = 1e-7

    drivers = drivers_from_params(p)  # AnnualCycle-based drivers that change over time

    # --- plot ---
    plot_panel(p, drivers, show=True)
    print("Plot generated with changing environment (seasonal drivers).")


if __name__ == "__main__":
    main()
