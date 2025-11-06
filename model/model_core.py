# model_core.py
from __future__ import annotations
import math
import numpy as np
from scipy.integrate import solve_ivp

# ------------------------------------------------------------
# Helpers to accept either constants/callables (old style) or arrays (new style)
# ------------------------------------------------------------
def _as_callable(x):
    """Return a callable f(t) for either a callable or a constant."""
    return x if callable(x) else (lambda t: float(x))

# ------------------------------------------------------------
# Two-substrate Synthesizing-Unit (Kooijman SU; matches your MM)
# ------------------------------------------------------------
def F(jm: float, A: float, B: float) -> float:
    """
    F[max_][A_,B_] := (A B (A+B) max) / (B^2 max + A^2 (B+max) + A B (B+max))
    """
    if jm <= 0.0 or A <= 0.0 or B <= 0.0:
        return 0.0
    num = A * B * (A + B) * jm
    den = (B * B) * jm + (A * A) * (B + jm) + (A * B) * (B + jm)
    if den <= 0.0:
        return 0.0
    return num / den

# ------------------------------------------------------------
# Core RHS with intermediates (no strings, no dict eval)
# y = [S, H, jCP, jSG]
# ------------------------------------------------------------
def compute_rhs_and_internals(
    S: float, H: float, jCP: float, jSG: float,
    L: float, T: float, X: float, Nu: float,
    p: dict[str, float],
):
    # --- unpack frequently used params into locals for speed/readability ---
    yC   = p["yC"];     nNH = p["nNH"];     yCL = p["yCL"]
    KX   = p["KX"];     KN  = p["KN"]
    jXm0 = p["jXm0"];   jNm0 = p["jNm0"];   jHGm0 = p["jHGm0"]
    jCPm0 = p["jCPm0"]; jSGm0 = p["jSGm0"]
    j0HT0 = p["j0HT0"]; j0ST0 = p["j0ST0"]
    sigmaNS = p["sigmaNS"]; nNS = p["nNS"]
    sigmaCH = p["sigmaCH"]; sigmaNH = p["sigmaNH"]; sigmaCS = p["sigmaCS"]
    nNX = p["nNX"]; #rNH_baseline = p.get("rNH", 0.0)
    astar = p["astar"]; kCO2 = p["kCO2"]; kNPQ = p["kNPQ"]; OneOverkROS = p["OneOverkROS"]
    Q10 = p["Q10"]; T0 = p["T0"]; sigmax = p["sigmax"]; steepness = p["steepness"]; ED50 = p["ED50"]
    lam = p.get("lam", 1.0); b = p.get("b", 5.0)

    # --- temperature scalings ---
    alpha = Q10 ** ((T - T0) / 10.0)
    sigT  = sigmax * math.exp(-steepness * (T  - ED50)) / (1.0 + math.exp(-steepness * (T  - ED50)))
    sigT0 = sigmax * math.exp(-steepness * (T0 - ED50)) / (1.0 + math.exp(-steepness * (T0 - ED50)))
    beta  = sigT / sigT0

    jHGm = alpha * jHGm0
    jNm  = alpha * jNm0
    jSGm = alpha * jSGm0
    j0HT = alpha * j0HT0
    j0ST = alpha * j0ST0
    jXm  = alpha * jXm0
    jCPm = alpha * beta * jCPm0

    # --- uptake and pools ---
    jX  = (jXm * X) / (X + KX) if (X + KX) != 0.0 else 0.0
    jN  = (jNm * Nu) / (Nu + KN) if (Nu + KN) != 0.0 else 0.0
    jHT = j0HT
    rNH = sigmaNH * nNH * jHT
    jNH = jN + nNX * jX + rNH #_baseline
    rNS = sigmaNS * nNS * j0ST  # constant in time at given T

    # --- geometry / light ---
    # A = 1.256307 + 1.385969 * exp(-6.479055 * S/H)
    A = 1.256307 + 1.385969 * math.exp(-6.479055 * (S / H)) if H > 0 else 1.256307
    jL = A * astar * L

    # --- partitioning from current states (uses current jCP, jSG) ---
    rhoC = jCP - jSG / yC

    # --- carbon handling up to jHG ---
    jHC = jX + (S / H) * rhoC if H > 0 else jX
    jHG = F(jHGm, yC * jHC, jNH / nNH)

    # --- downstream of jHG ---
    jeC  = jHC - jHG / yC
    jCO2 = jeC * kCO2
    rCH  = (jHT + (jHG * (1 - yC)) / yC) * sigmaCH

    # --- light partitioning and ROS ---
    jeL  = jL - jCP / yCL
    jNPQ = 0.0 if jeL == 0.0 else 1.0 / (1.0 / jeL + 1.0 / kNPQ)
    cROS1 = max(0.0, jeL - jNPQ) * OneOverkROS
    jST   = (1.0 + b * cROS1) * j0ST

    # --- finish partitioning with jHG known ---
    rhoN = jNH - jHG * nNH

    # --- aims/targets ---
    rCS    = sigmaCS * (j0ST + (1 - yC) * jSG / yC)
    jSGaim = F(jSGm, jCP * yC, (rNS + (H / S) * rhoN) / nNS) if S > 0 else F(jSGm, jCP * yC, rNS / nNS)
    jCPaim = F(jCPm, jL * yCL, rCS + (jCO2 + rCH) * (H / S if S > 0 else 0.0)) / (1.0 + cROS1)

    # --- dynamics ---
    dS   = S * (jSG - jST)
    dH   = H * (jHG - jHT)
    djCP = lam * (jCPaim - jCP)
    djSG = lam * (jSGaim - jSG)

    internals = {
        "jST": jST, "jHT": jHT, "jHG": jHG,
        "jCPaim": jCPaim, "jSGaim": jSGaim,
        "jeL": jeL, "jNPQ": jNPQ, "cROS1": cROS1,
        "rhoC": rhoC, "rhoN": rhoN, "A": A, "jL": jL,
        "jX": jX, "jN": jN, "jCO2": jCO2, "rCH": rCH,
    }
    return np.array([dS, dH, djCP, djSG], dtype=float), internals

# ------------------------------------------------------------
# Lightweight drivers with fast interpolants (array-based)
# ------------------------------------------------------------
class Drivers:
    """
    Provide (L, T, X, Nu) at any time t.
    Supply time series as (t_array, value_array) for each; they’ll be linearly interpolated.
    """
    def __init__(self,
                 tL: np.ndarray, L: np.ndarray,
                 tT: np.ndarray, T: np.ndarray,
                 tX: np.ndarray, X: np.ndarray,
                 tN: np.ndarray, Nu: np.ndarray):
        self.tL, self.L = np.asarray(tL), np.asarray(L)
        self.tT, self.T = np.asarray(tT), np.asarray(T)
        self.tX, self.X = np.asarray(tX), np.asarray(X)
        self.tN, self.Nu = np.asarray(tN), np.asarray(Nu)

    def tuple_at(self, t: float) -> tuple[float, float, float, float]:
        L  = float(np.interp(t, self.tL, self.L))
        T  = float(np.interp(t, self.tT, self.T))
        X  = float(np.interp(t, self.tX, self.X))
        Nu = float(np.interp(t, self.tN, self.Nu))
        return (L, T, X, Nu)

# ------------------------------------------------------------
# Adapter: convert old-style drivers (constants/callables) to array Drivers
# ------------------------------------------------------------
def drivers_from_fn_or_const(drivers_like, params: dict, t_eval: np.ndarray | None) -> Drivers:
    """
    Accepts an object with attributes .L, .T, .X, .Nu (each callable(t) or constant),
    samples them on t_eval, and returns an array-interpolated Drivers.
    If t_eval is None, it is built from params.
    """
    if isinstance(drivers_like, Drivers):
        return drivers_like

    # Build t_eval if absent
    if t_eval is None:
        tmax = float(params.get("tmax", 30.0))
        steps_per_day = int(params.get("steps_per_day", 48))
        t_eval = np.linspace(0.0, tmax, int(tmax * steps_per_day) + 1)

    fL  = _as_callable(getattr(drivers_like, "L"))
    fT  = _as_callable(getattr(drivers_like, "T"))
    fX  = _as_callable(getattr(drivers_like, "X"))
    fNu = _as_callable(getattr(drivers_like, "Nu"))

    L  = np.array([fL(t)  for t in t_eval], dtype=float)
    T  = np.array([fT(t)  for t in t_eval], dtype=float)
    X  = np.array([fX(t)  for t in t_eval], dtype=float)
    Nu = np.array([fNu(t) for t in t_eval], dtype=float)

    return Drivers(t_eval, L, t_eval, T, t_eval, X, t_eval, Nu)

# ------------------------------------------------------------
# Factory for RHS + simulator
# ------------------------------------------------------------
def make_rhs(params: dict[str, float], drivers: Drivers):
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        S, H, jCP, jSG = y
        L, T, X, Nu = drivers.tuple_at(t)
        d, _ = compute_rhs_and_internals(S, H, jCP, jSG, L, T, X, Nu, params)
        return d
    return rhs

def simulate(
    params: dict[str, float],
    drivers,  # can be old-style (constants/callables) or new Drivers
    y0: np.ndarray | list[float] = None,
    t_eval: np.ndarray | None = None,
    rtol: float = 1e-3,          # tighter
    atol: float = 1e-6,          # tighter
    method: str = "Radau",       # was DOP853; try Radau (or "BDF")
    max_step: float | None = None,  # cap steps if needed
):
    if y0 is None:
        y0 = params.get("y0", [1.0, 0.1, 0.2, 0.2])
    y0 = np.asarray(y0, dtype=float)

    # Build t_eval if needed
    tmax = float(params.get("tmax", 30.0))
    if t_eval is None:
        steps_per_day = int(params.get("steps_per_day", 48))
        t_eval = np.linspace(0.0, tmax, int(tmax * steps_per_day) + 1)

    # Adapt drivers after t_eval exists
    drivers = drivers_from_fn_or_const(drivers, params, t_eval)

    rhs_fun = make_rhs(params, drivers)

    sol = solve_ivp(
        fun=rhs_fun,
        t_span=(0.0, tmax),
        y0=y0,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        method=method,
        vectorized=False,
    )
    return sol

# ------------------------------------------------------------
# Post-processing: recompute outputs for plotting
# ------------------------------------------------------------
def compute_outputs_along_solution(sol, params: dict, drivers) -> dict[str, np.ndarray]:
    """
    Recompute selected intermediates along the solution for plotting.
    Accepts either old-style drivers (constants/callables) or the new array-based Drivers.
    """
    # 🔌 Ensure we have an array-interpolated Drivers aligned to sol.t
    drivers = drivers_from_fn_or_const(drivers, params, sol.t)

    t = sol.t
    S, H, jCP, jSG = sol.y

    out = {
        "t": t,
        "jCP": np.array(jCP, copy=True),
        "rho_C": np.empty_like(t, dtype=float),
        "rho_N": np.empty_like(t, dtype=float),
        "jeL": np.empty_like(t, dtype=float),
        "jNPQ": np.empty_like(t, dtype=float),
        "cROS_increase": np.empty_like(t, dtype=float),
        "expulsion": np.empty_like(t, dtype=float),
        "jSG": np.array(jSG, copy=True),
        "jHG": np.empty_like(t, dtype=float),
        "net_S": np.empty_like(t, dtype=float),
        "net_H": np.empty_like(t, dtype=float),
        "S_over_H": np.empty_like(t, dtype=float),
    }

    j0ST0 = params.get("j0ST0", 0.03)
    b = params.get("b", 5.0)

    for i, ti in enumerate(t):
        Li, Ti, Xi, Nui = drivers.tuple_at(float(ti))
        _, ints = compute_rhs_and_internals(
            S[i], H[i], jCP[i], jSG[i],
            Li, Ti, Xi, Nui,
            params,
        )
        out["rho_C"][i]        = ints["rhoC"]
        out["rho_N"][i]        = ints["rhoN"]
        out["jeL"][i]          = ints["jeL"]
        out["jNPQ"][i]         = ints["jNPQ"]
        out["cROS_increase"][i]= ints["cROS1"]
        out["expulsion"][i]    = b * ints["cROS1"] * j0ST0
        out["jHG"][i]          = ints["jHG"]
        out["net_S"][i]        = jSG[i] - ints["jST"]
        out["net_H"][i]        = ints["jHG"] - ints["jHT"]
        out["S_over_H"][i]     = (S[i] / H[i]) if H[i] != 0.0 else np.nan

    return out