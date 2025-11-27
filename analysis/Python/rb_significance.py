"""
rb_significance.py

rB significance analysis:
- ΔNF (mg C g-1 soil)
- ΔfnecC (% of SOC)
- truncation proportion
- R distribution (diagnostic ratio)

Outputs:
- data/derived/gsa/rB_significance_summary.csv
- data/derived/gsa/rB_R_distribution_long.csv

Author: A.K. Fine
"""

import os
import numpy as np
import pandas as pd

from model_definitions import (
    soils, EPS,
    CFB_emp, CFF_emp,
    CFB_fun, CFF_fun,
)


OUTDIR = "data/derived/gsa"
os.makedirs(OUTDIR, exist_ok=True)


def sample_CFF(n=20000, seed=1234):
    rng = np.random.default_rng(seed)
    u = np.clip(rng.random(n), 1e-12, 1 - 1e-12)
    try:
        return np.quantile(CFF_emp, u, method="linear")
    except TypeError:
        return np.quantile(CFF_emp, u, interpolation="linear")


def NF_with_rB(GS, MS, CFF, rB):
    """
    Fungal necromass C (mg C g-1 soil) with bacterial correction.
    """
    GS_corr = np.maximum(GS - rB * MS, 0.0)
    return GS_corr * CFF


def fnecC_with_rB(GS, MS, S, CFF, CFB, rB):
    """
    Total necromass fraction of SOC (%) using aggregated form and given CFB/CFF.
    """
    S = np.maximum(S, EPS)
    GS_corr = np.maximum(GS - rB * MS, 0.0)
    return 100.0 * (MS * CFB + GS_corr * CFF) / S


def bootstrap_ci(x, nboot=2000, subsample=200000, seed=42, q=(2.5, 97.5)):
    """
    Memory-safe bootstrap CI for the median.

    Draws 'subsample' items (with replacement) from x for each replicate,
    computes the median, and returns median(x) and the percentile CI of medians.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    N = x.size
    m = int(min(subsample, N))
    meds = np.empty(nboot, dtype=float)
    for b in range(nboot):
        idx = rng.integers(0, N, size=m)
        meds[b] = np.median(x[idx])
    lo, hi = np.percentile(meds, q)
    return float(np.median(x)), float(lo), float(hi)


def rB_significance_analysis(rB_values=(0.0, 1.0, 1.16, 2.0, 3.0, 4.0),
                             n_cff=20000, seed=1234, outdir=OUTDIR,
                             eps_C=0.001, eps_pct=0.5):
    os.makedirs(outdir, exist_ok=True)

    # soil vectors
    GS = soils["GS"].to_numpy(dtype=float)
    MS = soils["MS"].to_numpy(dtype=float)
    S  = soils["S"].to_numpy(dtype=float)

    # empirical CF samples
    CFF_s = sample_CFF(n=n_cff, seed=seed)
    rng = np.random.default_rng(seed + 1)
    u_cfb = np.clip(rng.random(n_cff), 1e-12, 1 - 1e-12)
    try:
        CFB_s = np.quantile(CFB_emp, u_cfb, method="linear")
    except TypeError:
        CFB_s = np.quantile(CFB_emp, u_cfb, interpolation="linear")

    # broadcast to soil × CF grids
    GS2, MS2, S2 = GS[:, None], MS[:, None], S[:, None]
    CFF2, CFB2   = CFF_s[None, :], CFB_s[None, :]

    # baseline (rB = 0)
    NF0   = NF_with_rB(GS2, MS2, CFF2, rB=0.0)
    fnec0 = fnecC_with_rB(GS2, MS2, S2, CFF2, CFB2, rB=0.0)

    rows = []
    R_summaries = []

    for rB in rB_values:
        # with correction
        NF_rB   = NF_with_rB(GS2, MS2, CFF2, rB=rB)
        fnec_rB = fnecC_with_rB(GS2, MS2, S2, CFF2, CFB2, rB=rB)

        # paired differences
        dNF   = (NF_rB - NF0).ravel()
        dfnec = (fnec_rB - fnec0).ravel()

        # bootstrap medians and 95% CIs
        dNF_med, dNF_lo, dNF_hi       = bootstrap_ci(dNF,   nboot=10000, seed=seed+2)
        dfnec_med, dfnec_lo, dfnec_hi = bootstrap_ci(dfnec, nboot=10000, seed=seed+3)

        # truncation proportion over soils (independent of CF)
        trunc = float(np.mean((GS - rB * MS) <= 0))

        # diagnostic R distribution
        denom  = (CFF2 * np.maximum(GS2, EPS))
        R_vals = (100.0 * (rB * MS2) / denom)
        R_vals = R_vals[np.isfinite(R_vals)].ravel()
        R_med, R_lo, R_hi = bootstrap_ci(R_vals, nboot=10000, seed=seed+4)

        row = {
            "rB": float(rB),
            "dNF_median_mgC_g": dNF_med,
            "dNF_lo95": dNF_lo,
            "dNF_hi95": dNF_hi,
            "dfnec_median_pct": dfnec_med,
            "dfnec_lo95": dfnec_lo,
            "dfnec_hi95": dfnec_hi,
            "truncation_prop": trunc,
            "R_median_pct": R_med,
            "R_lo95": R_lo,
            "R_hi95": R_hi,
            "equiv_NF":    (dNF_lo   >= -eps_C)   and (dNF_hi   <= eps_C),
            "equiv_fnecC": (dfnec_lo >= -eps_pct) and (dfnec_hi <= eps_pct),
        }
        rows.append(row)

        # full R distribution for supplementary figure
        R_summaries.append(pd.DataFrame({"rB": rB, "R_percent": R_vals}))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(outdir, "rB_significance_summary.csv"), index=False)

    R_all = pd.concat(R_summaries, ignore_index=True)
    R_all.to_csv(os.path.join(outdir, "rB_R_distribution_long.csv"), index=False)

    return out


if __name__ == "__main__":
    print("Running rB significance analysis...")
    df = rB_significance_analysis()
    print("Done. Results written to data/derived/gsa/")
