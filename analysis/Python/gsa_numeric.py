"""
gsa_numeric.py

Global sensitivity analysis (Sobol, Delta, Morris) for:
- CF_B
- CF_F
- N_F (fungal necromass C)
- fnecC (aggregated and integrated forms)

Outputs are tidy CSV tables under data/derived/gsa/.

Author: A.K. Fine
"""

import os
import numpy as np
import pandas as pd

from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol
from SALib.analyze.delta import analyze as delta_analyze
from SALib.sample.morris import sample as morris_sample
from SALib.analyze.morris import analyze as morris_analyze

from model_definitions import (
    soils,
    model_CFB,     problem_CFB,
    model_CFF,     problem_CFF,
    model_fungNecC, problem_fung,
    model_fnecC5,  problem_fnecC5,
    model_fnecC10, problem_fnecC10,
)

OUTDIR = "data/derived/gsa"
os.makedirs(OUTDIR, exist_ok=True)


# --------------------------------------------------------------------
# 1. Helper utilities
# --------------------------------------------------------------------
def _next_pow2(n):
    n = int(n)
    return 1 << (n - 1).bit_length()


def run_sobol(problem, model_fun, N=4096, seed=42,
              calc_second_order=False, num_resamples=500, conf_level=0.95):
    """
    Run Sobol analysis with bootstrapped confidence intervals.
    Returns (Si, extra), where extra may hold resampled indices.
    """
    np.random.seed(seed)
    N = _next_pow2(N)
    U = sobol_sample.sample(problem, N, calc_second_order=calc_second_order)
    Y = model_fun(U)

    extra = {}
    try:
        Si = sobol.analyze(
            problem, Y,
            calc_second_order=calc_second_order,
            conf_level=conf_level,
            num_resamples=num_resamples,
            print_to_console=False,
            keep_resamples=True,
        )
        for k in ["S1_resampled", "ST_resampled"]:
            if k in Si and Si[k] is not None:
                extra[k] = np.asarray(Si[k])
    except TypeError:
        Si = sobol.analyze(
            problem, Y,
            calc_second_order=calc_second_order,
            conf_level=conf_level,
            num_resamples=num_resamples,
            print_to_console=False,
        )
    return Si, extra


def summarize_indices(Si, names, boot=None, conf_level=0.95):
    """
    Summarize Sobol indices into median and 95 % interval.
    """
    alpha = 1 - conf_level
    lo_q, hi_q = 100 * (alpha / 2), 100 * (1 - alpha / 2)

    def _summ(one, resampled=None, conf_half=None):
        one = np.asarray(one, dtype=float)
        if resampled is not None and np.size(resampled) > 0:
            med = np.median(resampled, axis=0)
            lo  = np.percentile(resampled, lo_q, axis=0)
            hi  = np.percentile(resampled, hi_q, axis=0)
        else:
            med = one
            lo  = one - (conf_half if conf_half is not None else 0)
            hi  = one + (conf_half if conf_half is not None else 0)
        return med, lo, hi

    S1 = Si["S1"]
    ST = Si["ST"]
    S1_conf = Si.get("S1_conf", np.full_like(S1, np.nan))
    ST_conf = Si.get("ST_conf", np.full_like(ST, np.nan))
    S1_res = boot.get("S1_resampled") if boot else None
    ST_res = boot.get("ST_resampled") if boot else None

    S1_med, S1_lo, S1_hi = _summ(S1, resampled=S1_res, conf_half=S1_conf)
    ST_med, ST_lo, ST_hi = _summ(ST, resampled=ST_res, conf_half=ST_conf)

    return pd.DataFrame({
        "parameter": names,
        "S1_median": S1_med, "S1_lo95": S1_lo, "S1_hi95": S1_hi,
        "ST_median": ST_med, "ST_lo95": ST_lo, "ST_hi95": ST_hi,
    })


def run_delta(problem, model_fun, N=4096, seed=123):
    np.random.seed(seed)
    X = sobol_sample.sample(problem, N, calc_second_order=False)
    Y = model_fun(X)
    return delta_analyze(problem, X, Y, print_to_console=False)


def export_delta(Si, names, fname):
    out = pd.DataFrame({
        "parameter": names,
        "delta": Si["delta"],
        "delta_conf": Si.get("delta_conf", np.full(len(names), np.nan, dtype=float)),
    })
    out.to_csv(fname, index=False)
    return out


def run_morris(problem, model_fun, N_trajectories=1000, num_levels=8,
               optimal_trajectories=10, seed=321):
    np.random.seed(seed)
    X = morris_sample(problem, N=N_trajectories,
                      num_levels=num_levels,
                      optimal_trajectories=optimal_trajectories)
    Y = model_fun(X)
    return morris_analyze(problem, X, Y, conf_level=0.95, print_to_console=False)


def export_morris(Si, names, fname):
    out = pd.DataFrame({
        "parameter": names,
        "mu": Si["mu"],
        "mu_star": Si["mu_star"],
        "mu_star_conf": Si.get("mu_star_conf", np.full(len(names), np.nan)),
        "sigma": Si["sigma"],
    })
    out.to_csv(fname, index=False)
    return out


# --------------------------------------------------------------------
# 2. Main analysis: Sobol, Delta, Morris, CSV outputs
# --------------------------------------------------------------------
def summarize_and_save(Si, boot, names, out_path):
    tab = summarize_indices(Si, names, boot)
    tab.to_csv(out_path, index=False)
    return tab


def main(N=4096):
    print("Running Sobol analyses...")
    Si_CFB,     B_CFB     = run_sobol(problem_CFB,     model_CFB,      N=N)
    Si_CFF,     B_CFF     = run_sobol(problem_CFF,     model_CFF,      N=N)
    Si_fung,    B_fung    = run_sobol(problem_fung,    model_fungNecC, N=N)
    Si_fnecC5,  B_fnecC5  = run_sobol(problem_fnecC5,  model_fnecC5,   N=N)
    Si_fnecC10, B_fnecC10 = run_sobol(problem_fnecC10, model_fnecC10,  N=N)
    print("Sobol analyses completed.")

    print("Running Delta (N_F, fnecC_agg) and Morris (fnecC_int)...")
    Si_delta_fung   = run_delta(problem_fung,   model_fungNecC, N=N)
    tab_delta_fung  = export_delta(Si_delta_fung, problem_fung["names"],
                                   os.path.join(OUTDIR, "delta_NF.csv"))

    Si_delta_fnec5  = run_delta(problem_fnecC5, model_fnecC5,   N=N)
    tab_delta_fnec5 = export_delta(Si_delta_fnec5, problem_fnecC5["names"],
                                   os.path.join(OUTDIR, "delta_fnecC_agg.csv"))

    Si_morris_fnec10  = run_morris(problem_fnecC10, model_fnecC10,
                                   N_trajectories=1500, num_levels=8, optimal_trajectories=10)
    tab_morris_fnec10 = export_morris(Si_morris_fnec10, problem_fnecC10["names"],
                                      os.path.join(OUTDIR, "morris_fnecC_int.csv"))

    print("Exporting tidy Sobol tables with CIs...")
    tab_CFB     = summarize_and_save(Si_CFB,     B_CFB,     problem_CFB["names"],
                                     os.path.join(OUTDIR, "sobol_CFB_bootstrap.csv"))
    tab_CFF     = summarize_and_save(Si_CFF,     B_CFF,     problem_CFF["names"],
                                     os.path.join(OUTDIR, "sobol_CFF_bootstrap.csv"))
    tab_fung    = summarize_and_save(Si_fung,    B_fung,    problem_fung["names"],
                                     os.path.join(OUTDIR, "sobol_NF_bootstrap.csv"))
    tab_fnecC5  = summarize_and_save(Si_fnecC5,  B_fnecC5,  problem_fnecC5["names"],
                                     os.path.join(OUTDIR, "sobol_fnecC_agg_bootstrap.csv"))
    tab_fnecC10 = summarize_and_save(Si_fnecC10, B_fnecC10, problem_fnecC10["names"],
                                     os.path.join(OUTDIR, "sobol_fnecC_int_bootstrap.csv"))

    print("Done. Sobol / Delta / Morris outputs written to:", OUTDIR)


if __name__ == "__main__":
    main()
