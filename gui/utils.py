"""Data loading utilities for the SCOPE2 sensitivity GUI."""

from pathlib import Path
import numpy as np
import pandas as pd

LUT_DIR = Path(__file__).parent / "lut_data"

PARAM_META = {
    "Rin":   {"label": "Incoming shortwave radiation",  "units": "W m⁻²",          "default": 600},
    "Ta":    {"label": "Air temperature",               "units": "°C",             "default": 20},
    "LAI":   {"label": "Total leaf area index",         "units": "m² m⁻²",         "default": 3.0},
    "Vcmo":  {"label": "Max carboxylation rate (Vcmax)", "units": "μmol m⁻² s⁻¹", "default": 60},
    "Ca":    {"label": "Atmospheric CO₂ concentration", "units": "ppm",            "default": 410},
    "ea":    {"label": "Vapour pressure",               "units": "hPa",            "default": 15},
    "LIDFa": {"label": "Leaf inclination (LIDFa)",      "units": "—",              "default": -0.35},
}

OUTPUT_META = {
    "Actot": {"label": "Net canopy photosynthesis", "units": "μmol CO₂ m⁻² s⁻¹"},
    "aPAR":  {"label": "Absorbed PAR",              "units": "μmol photons m⁻² s⁻¹"},
    "Rd":    {"label": "Dark respiration",          "units": "μmol CO₂ m⁻² s⁻¹"},
}

PROFILES = ["Uniform", "Top-heavy", "Bottom-heavy"]

PROFILE_COLORS = {
    "Uniform":     "#1f77b4",
    "Top-heavy":   "#2ca02c",
    "Bottom-heavy":"#d62728",
}

# ── Vertical profiles panel ───────────────────────────────────────────────────

ANGLE_SCENARIOS = ["Uniform", "Vertical-on-top", "Horizontal-on-top"]

ANGLE_COLORS = {
    "Uniform":           "#1f77b4",
    "Vertical-on-top":   "#2ca02c",
    "Horizontal-on-top": "#d62728",
}

ANGLE_DESCRIPTIONS = {
    "Uniform":           "Spherical (LIDFa = 0) at every layer",
    "Vertical-on-top":   "LIDFa: +0.4 (top) → −0.4 (bottom) — erectophile near sky",
    "Horizontal-on-top": "LIDFa: −0.4 (top) → +0.4 (bottom) — planophile near sky",
}

PROFILE_METRICS = {
    "A":    {"label": "Net photosynthesis (A)",  "units": "μmol CO₂ m⁻² s⁻¹"},
    "aPAR": {"label": "Absorbed PAR",            "units": "μmol photons m⁻² s⁻¹"},
    "Rd":   {"label": "Dark respiration (Rd)",   "units": "μmol CO₂ m⁻² s⁻¹"},
    "Ps":   {"label": "Sunlit fraction",         "units": "—"},
    "Kn":   {"label": "NPQ rate constant (Kn)",  "units": "s⁻¹"},
}

PROFILE_SWEEP_VALUES = {
    "LAI": [0.5, 1.5, 3.0, 5.0, 8.0],
    "Rin": [100.0, 300.0, 600.0, 900.0, 1200.0],
}


def lut_ready() -> bool:
    return (LUT_DIR / "wavelengths.csv").exists()


def profiles_ready() -> bool:
    return (LUT_DIR / "vertical_profiles.csv").exists()


def load_vertical_profiles(sweep_param: str, sweep_value: float) -> pd.DataFrame:
    """Return per-layer data for all 3 angle scenarios at the nearest sweep value."""
    df = pd.read_csv(LUT_DIR / "vertical_profiles.csv")
    avail = sorted(df[df["sweep_param"] == sweep_param]["sweep_value"].unique())
    closest = min(avail, key=lambda v: abs(v - sweep_value))
    return df[(df["sweep_param"] == sweep_param) & (df["sweep_value"] == closest)]


def load_wavelengths() -> np.ndarray:
    return pd.read_csv(LUT_DIR / "wavelengths.csv", header=None).values.flatten()


def load_sensitivity(param: str, profile: str) -> pd.DataFrame:
    """Return rows for one (param, profile) combination."""
    path = LUT_DIR / f"sensitivity_{param}.csv"
    df = pd.read_csv(path)
    return df[df["profile"] == profile].copy().reset_index(drop=True)


def load_sensitivity_all_profiles(param: str) -> pd.DataFrame:
    """Return all rows for a parameter sweep across all profiles."""
    path = LUT_DIR / f"sensitivity_{param}.csv"
    return pd.read_csv(path)


def load_reflectance(param: str, profile: str, step: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Return (wavelengths, matrix) for spectral overlay plots.

    step: take every `step`-th row to thin the number of spectra shown.
    Returns wavelengths (n_wl,) and reflectance (n_runs, n_wl).
    """
    path = LUT_DIR / f"reflectance_{param}.csv"
    df   = pd.read_csv(path)
    df   = df[df["profile"] == profile].iloc[::step].reset_index(drop=True)
    wl   = load_wavelengths()
    wl_cols = [c for c in df.columns if c.startswith("wl_")]
    refl = df[wl_cols].values
    return wl, refl, df["value"].values
