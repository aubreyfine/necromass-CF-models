"""
model_definitions.py

Core data loading, trait distributions, and CF / fnecC model functions
used in the necromass-CF uncertainty analysis.

Author: A.K. Fine
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

# --------------------------------------------------------------------
# 1. Robust CSV reader
# --------------------------------------------------------------------
def read_csv_robust(path):
    """
    Read CSV with a small set of fallback encodings and delimiter detection.
    """
    for enc in ("utf-8", "utf-8-sig", "cp1252", "latin1"):
        try:
            return pd.read_csv(path, sep=None, engine="python", encoding=enc)
        except UnicodeDecodeError:
            continue
        except Exception:
            continue
    return pd.read_csv(path, sep=None, engine="python",
                       encoding="latin1", encoding_errors="replace")


# --------------------------------------------------------------------
# 2. Load processed data (paths relative to repo root)
# --------------------------------------------------------------------
P_SOILS = "data/processed/soil.data.csv"
P_BACT  = "data/processed/bactMurA.csv"
P_FUNG  = "data/processed/fungGlcN.csv"

if not (os.path.exists(P_SOILS) and os.path.exists(P_BACT) and os.path.exists(P_FUNG)):
    raise FileNotFoundError(
        "One or more processed data files not found. "
        "Expected:\n"
        f"  {P_SOILS}\n"
        f"  {P_BACT}\n"
        f"  {P_FUNG}\n"
        "Make sure processed inputs are present under data/processed/."
    )

soil_raw = read_csv_robust(P_SOILS)
bactMurA = read_csv_robust(P_BACT)
fungGlcN = read_csv_robust(P_FUNG)

# Filter extreme MurA / GlcN, as in CFupdate_v3.py
soil_f = soil_raw.loc[
    (soil_raw["MurA"] <= 412) & (soil_raw["GlcN"] <= 7800)
].copy()

# Compute GS, MS, S (mg g-1)
soils = (
    soil_f.assign(
        GS = soil_f["GlcN"] / 1000.0,
        MS = soil_f["MurA"] / 1000.0,
        S  = soil_f["SOC"]
    )[["MS", "GS", "S"]]
    .replace([np.inf, -np.inf], np.nan)
    .dropna()
)
assert all(c in soils.columns for c in ["MS", "GS", "S"]), "MS/GS/S not found in soils."

# Trait vectors
traits_bact = bactMurA[["Gram", "MurA"]].dropna()
murA_GP_vec = traits_bact.loc[traits_bact["Gram"].str.upper() == "GP", "MurA"].astype(float).values
murA_GN_vec = traits_bact.loc[traits_bact["Gram"].str.upper() != "GP", "MurA"].astype(float).values
glcN_F_vec  = fungGlcN["GlcN"].dropna().astype(float).values


# --------------------------------------------------------------------
# 3. Soil marginals and quantile functions
# --------------------------------------------------------------------
def fit_lognorm_params(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x) & (x > 0)]
    lx = np.log(x)
    return float(lx.mean()), float(lx.std(ddof=1))


ms_meanlog, ms_sdlog = fit_lognorm_params(soils["MS"])
gs_meanlog, gs_sdlog = fit_lognorm_params(soils["GS"])
S_shape, S_loc, S_scale = stats.gamma.fit(soils["S"], floc=0)

def qMS(u):
    u = np.clip(u, 1e-12, 1 - 1e-12)
    return stats.lognorm(s=ms_sdlog, scale=np.exp(ms_meanlog)).ppf(u)

def qGS(u):
    u = np.clip(u, 1e-12, 1 - 1e-12)
    return stats.lognorm(s=gs_sdlog, scale=np.exp(gs_meanlog)).ppf(u)

def qS(u):
    u = np.clip(u, 1e-12, 1 - 1e-12)
    return stats.gamma(a=S_shape, loc=0, scale=S_scale).ppf(u)

def q_emp(u, arr):
    """
    Empirical quantile function with robust clipping and NaN handling.
    """
    u = np.asarray(u, dtype=float)
    u = np.clip(u, 1e-12, 1 - 1e-12)
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.full_like(u, np.nan, dtype=float)
    try:
        return np.quantile(arr, u, method="linear")
    except TypeError:
        # Older NumPy fallback
        return np.quantile(arr, u, interpolation="linear")


# --------------------------------------------------------------------
# 4. CF / necromass model functions
# --------------------------------------------------------------------
EPS      = 1e-9
RB_CONST = 1.16  # profiled r_B for aggregated model (main text value)

def CFB_fun(cB, MGP, MGN, fGP):
    """
    Bacterial conversion factor CF_B given trait inputs.
    """
    den = np.maximum(MGP * fGP + MGN * (1.0 - fGP), EPS)
    return cB / den

def CFF_fun(cF, GF):
    """
    Fungal conversion factor CF_F given trait inputs.
    """
    GF = np.maximum(GF, EPS)
    return cF / GF

def fungNecC_fun(GS, MS, cF, GF, rB):
    """
    Fungal necromass C (mg C g-1 soil) using GS, MS, and a bacterial correction term rB.
    """
    CFF = CFF_fun(cF, GF)
    return np.maximum(GS - rB * MS, 0.0) * CFF


# --------------------------------------------------------------------
# 5. fnecC models: 10-input (integrated) and 5-input (aggregated)
# --------------------------------------------------------------------
def model_fnecC10(U):
    """
    fnecC (%, integrated 10-input model) evaluated at Sobol samples U in [0,1]^10.
    """
    MS  = qMS(U[:, 0])
    GS  = qGS(U[:, 1])
    S   = qS(U[:, 2])
    cB  = 300.0 + 400.0 * U[:, 3]
    MGP = q_emp(U[:, 4], murA_GP_vec)
    MGN = q_emp(U[:, 5], murA_GN_vec)
    fGP = 0.1 + 0.8 * U[:, 6]
    cF  = 300.0 + 400.0 * U[:, 7]
    GF  = q_emp(U[:, 8], glcN_F_vec)
    rB  = 4.0 * U[:, 9]

    CFB = CFB_fun(cB, MGP, MGN, fGP)
    CFF = CFF_fun(cF, GF)

    GS_corr = np.maximum(GS - rB * MS, 0.0)
    S = np.maximum(S, EPS)
    totalNec = MS * CFB + GS_corr * CFF
    return 100.0 * totalNec / S

problem_fnecC10 = {
    "num_vars": 10,
    "names": ["MS", "GS", "S", "cB", "MGP", "MGN", "fGP", "cF", "GF", "rB"],
    "bounds": [[0, 1]] * 10,
}


# Empirical CF_B and CF_F quantiles
def build_empirical_CFs(n=200000, seed=123):
    """
    Monte Carlo sampling to build empirical CF_B and CF_F distributions
    from trait marginals.
    """
    rng = np.random.default_rng(seed)
    U = rng.random((n, 7))  # [cB, MGP, MGN, fGP, cF, GF, dummy]

    cB  = 300.0 + U[:, 0] * 400.0
    MGP = q_emp(U[:, 1], murA_GP_vec)
    MGN = q_emp(U[:, 2], murA_GN_vec)
    fGP = 0.1 + U[:, 3] * 0.8
    CFB = CFB_fun(cB, MGP, MGN, fGP)

    cF  = 300.0 + U[:, 4] * 400.0
    GF  = q_emp(U[:, 5], glcN_F_vec)
    CFF = CFF_fun(cF, GF)

    return np.sort(CFB), np.sort(CFF)

CFB_emp, CFF_emp = build_empirical_CFs()

def qCFB(u):
    u = np.clip(u, 1e-12, 1 - 1e-12)
    try:
        return np.quantile(CFB_emp, u, method="linear")
    except TypeError:
        return np.quantile(CFB_emp, u, interpolation="linear")

def qCFF(u):
    u = np.clip(u, 1e-12, 1 - 1e-12)
    try:
        return np.quantile(CFF_emp, u, method="linear")
    except TypeError:
        return np.quantile(CFF_emp, u, interpolation="linear")


def fnecC5_fun(MS, GS, S, CFB, CFF, rB_const=RB_CONST):
    """
    Aggregated fnecC (%) using MS, GS, S and empirical CF_B, CF_F plus fixed rB.
    """
    S = np.maximum(S, EPS)
    GS_corr = np.maximum(GS - rB_const * MS, 0.0)
    return 100.0 * (MS * CFB + GS_corr * CFF) / S

def model_fnecC5(U, rB_const=RB_CONST):
    MS  = qMS(U[:, 0])
    GS  = qGS(U[:, 1])
    S   = qS(U[:, 2])
    CFB = qCFB(U[:, 3])
    CFF = qCFF(U[:, 4])
    return fnecC5_fun(MS, GS, S, CFB, CFF, rB_const)

problem_fnecC5 = {
    "num_vars": 5,
    "names": ["MS", "GS", "S", "CFB", "CFF"],
    "bounds": [[0, 1]] * 5,
}


# --------------------------------------------------------------------
# 6. Component models (CF_B, CF_F, N_F)
# --------------------------------------------------------------------
def model_CFB(U):
    cB  = 300.0 + U[:, 0] * 400.0
    MGP = q_emp(U[:, 1], murA_GP_vec)
    MGN = q_emp(U[:, 2], murA_GN_vec)
    fGP = 0.1 + U[:, 3] * 0.8
    return CFB_fun(cB, MGP, MGN, fGP)

problem_CFB = {
    "num_vars": 4,
    "names": ["cB", "MGP", "MGN", "fGP"],
    "bounds": [[0, 1]] * 4,
}

def model_CFF(U):
    cF = 300.0 + U[:, 0] * 400.0
    GF = q_emp(U[:, 1], glcN_F_vec)
    return CFF_fun(cF, GF)

problem_CFF = {
    "num_vars": 2,
    "names": ["cF", "GF"],
    "bounds": [[0, 1]] * 2,
}

def model_fungNecC(U):
    GS = qGS(U[:, 0])
    MS = qMS(U[:, 1])
    cF = 300.0 + U[:, 2] * 400.0
    GF = q_emp(U[:, 3], glcN_F_vec)
    rB = U[:, 4] * 4.0
    return fungNecC_fun(GS, MS, cF, GF, rB)

problem_fung = {
    "num_vars": 5,
    "names": ["GS", "MS", "cF", "GF", "rB"],
    "bounds": [[0, 1]] * 5,
}


# --------------------------------------------------------------------
# 7. Pretty parameter labels (for downstream plotting)
# --------------------------------------------------------------------
param_labels = {
    "cB":  r"$\mathit{c}_{B}$", "cF":  r"$\mathit{c}_{F}$",
    "fGP": r"$\mathit{f}_{GP}$", "GF":  r"$\mathit{G}_{F}$",
    "GS":  r"$\mathit{G}_{S}$",  "MGN": r"$\mathit{M}_{GN}$",
    "MGP": r"$\mathit{M}_{GP}$", "MS":  r"$\mathit{M}_{S}$",
    "rB":  r"$\mathit{r}_{B}$",  "S":   r"$\mathit{S}$",
    "CFB": r"$CF_{B}$",          "CFF": r"$CF_{F}$",
}

def pretty_names(names):
    return [param_labels.get(n, n) for n in names]
