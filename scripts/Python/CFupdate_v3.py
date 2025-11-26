"""
CFupdate_v2.py
Global sensitivity analysis (Sobol indices) for CF_B, CF_F, N_F, and f_necC models
(aggregated with profiled r_B; integrated with trait-level expansion)
Author: A.K. Fine
"""

# ============================================================
# 0. Import dependencies
# ============================================================
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless environments

# NEW: use sobol sampler (avoid deprecation)
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol
from SALib.analyze.delta import analyze as delta_analyze
from SALib.sample.morris import sample as morris_sample
from SALib.analyze.morris import analyze as morris_analyze

OUTDIR = "data/derived/gsa"
os.makedirs(OUTDIR, exist_ok=True)

sobol_colors = {"S1": "#336392", "ST": "#a94322"}

# ============================================================
# 1. Load and preprocess data
# ============================================================
def read_csv_robust(path):
    for enc in ("utf-8", "utf-8-sig", "cp1252", "latin1"):
        try:
            return pd.read_csv(path, sep=None, engine="python", encoding=enc)
        except UnicodeDecodeError:
            continue
        except Exception:
            continue
    return pd.read_csv(path, sep=None, engine="python",
                       encoding="latin1", encoding_errors="replace")

# ---- File paths (edit if needed) ----
p_soils = "data/processed/soil.data.csv"
p_bact  = "data/processed/bactMurA.csv"
p_fung  = "data/processed/fungGlcN.csv"

soil_raw = read_csv_robust(p_soils)
bactMurA = read_csv_robust(p_bact)
fungGlcN = read_csv_robust(p_fung)

# ---- Filter outliers (MurA ≤ 412 µg g-1 & GlcN ≤ 7800 µg g-1) ----
soil_f = soil_raw.loc[(soil_raw["MurA"] <= 412) & (soil_raw["GlcN"] <= 7800)].copy()

# ---- Compute GS, MS, and S (mg g-1) ----
soils = (
    soil_f.assign(
        GS = soil_f["GlcN"] / 1000.0,
        MS = soil_f["MurA"] / 1000.0,
        S  = soil_f["SOC"]
    )[["MS", "GS", "S"]]
    .replace([np.inf, -np.inf], np.nan)
    .dropna()
)
assert all(c in soils.columns for c in ["MS","GS","S"]), "MS/GS/S not found."

# ---- Trait vectors ----
traits_bact = bactMurA[["Gram","MurA"]].dropna()
murA_GP_vec = traits_bact.loc[traits_bact["Gram"].str.upper() == "GP", "MurA"].astype(float).values
murA_GN_vec = traits_bact.loc[traits_bact["Gram"].str.upper() != "GP", "MurA"].astype(float).values
glcN_F_vec  = fungGlcN["GlcN"].dropna().astype(float).values

# ============================================================
# 2. Fit marginals for soils
# ============================================================
def fit_lognorm_params(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x) & (x > 0)]
    lx = np.log(x)
    return float(lx.mean()), float(lx.std(ddof=1))

ms_meanlog, ms_sdlog = fit_lognorm_params(soils["MS"])
gs_meanlog, gs_sdlog = fit_lognorm_params(soils["GS"])
S_shape, S_loc, S_scale = stats.gamma.fit(soils["S"], floc=0)

# --- FIX: clip u before ppf to avoid 0/1 ---
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
    u = np.asarray(u, dtype=float)
    u = np.clip(u, 1e-12, 1 - 1e-12)
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.full_like(u, np.nan, dtype=float)
    try:
        return np.quantile(arr, u, method="linear")
    except TypeError:
        return np.quantile(arr, u, interpolation="linear")

def _next_pow2(n):
    n = int(n)
    return 1 << (n - 1).bit_length()

# ============================================================
# 3. Model functions
# ============================================================
EPS = 1e-9
RB_CONST = 1.16  # profiled r_B for aggregated model (sweep below if desired)

def CFB_fun(cB, MGP, MGN, fGP):
    den = np.maximum(MGP * fGP + MGN * (1 - fGP), EPS)
    return cB / den

def CFF_fun(cF, GF):
    GF = np.maximum(GF, EPS)
    return cF / GF

def fungNecC_fun(GS, MS, cF, GF, rB):
    CFF = CFF_fun(cF, GF)
    return np.maximum(GS - rB * MS, 0.0) * CFF

# ======== fnecC (10-input, integrated) ======================
def model_fnecC10(U):
    MS  = qMS(U[:,0])
    GS  = qGS(U[:,1])
    S   = qS(U[:,2])
    cB  = 300.0 + 400.0 * U[:,3]
    MGP = q_emp(U[:,4], murA_GP_vec)
    MGN = q_emp(U[:,5], murA_GN_vec)
    fGP = 0.1 + 0.8 * U[:,6]
    cF  = 300.0 + 400.0 * U[:,7]
    GF  = q_emp(U[:,8], glcN_F_vec)
    rB  = 4.0 * U[:,9]

    CFB = CFB_fun(cB, MGP, MGN, fGP)
    CFF = CFF_fun(cF, GF)

    GS_corr = np.maximum(GS - rB * MS, 0.0)
    S = np.maximum(S, EPS)
    totalNec = MS * CFB + GS_corr * CFF
    return 100.0 * totalNec / S

problem_fnecC10 = {
    "num_vars": 10,
    "names": ["MS","GS","S","cB","MGP","MGN","fGP","cF","GF","rB"],
    "bounds": [[0,1]] * 10
}

# ============================================================
# A. Empirical CF_B and CF_F quantiles
# ============================================================
def build_empirical_CFs(n=200000, seed=123):
    rng = np.random.default_rng(seed)
    U = rng.random((n, 7))  # [cB, MGP, MGN, fGP, cF, GF, dummy]

    cB  = 300 + U[:,0]*400
    MGP = q_emp(U[:,1], murA_GP_vec)
    MGN = q_emp(U[:,2], murA_GN_vec)
    fGP = 0.1 + U[:,3]*0.8
    CFB = CFB_fun(cB, MGP, MGN, fGP)

    cF  = 300 + U[:,4]*400
    GF  = q_emp(U[:,5], glcN_F_vec)
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
    # FIXED typo here
    u = np.clip(u, 1e-12, 1 - 1e-12)
    try:
        return np.quantile(CFF_emp, u, method="linear")
    except TypeError:
        return np.quantile(CFF_emp, u, interpolation="linear")

# ============================================================
# B. Aggregated fnecC with profiled r_B (5 inputs)
# ============================================================
def fnecC5_fun(MS, GS, S, CFB, CFF, rB_const=RB_CONST):
    S = np.maximum(S, EPS)
    GS_corr = np.maximum(GS - rB_const * MS, 0.0)
    return 100.0 * (MS * CFB + GS_corr * CFF) / S

def model_fnecC5(U, rB_const=RB_CONST):
    MS  = qMS(U[:,0])
    GS  = qGS(U[:,1])
    S   = qS(U[:,2])
    CFB = qCFB(U[:,3])
    CFF = qCFF(U[:,4])
    return fnecC5_fun(MS, GS, S, CFB, CFF, rB_const)

problem_fnecC5 = {
    "num_vars": 5,
    "names": ["MS","GS","S","CFB","CFF"],
    "bounds": [[0,1]]*5
}

# ============================================================
# 4. Sobol/Delta/Morris helpers
# ============================================================
def run_sobol(problem, model_fun, N=4096, seed=42, calc_second_order=False,
              num_resamples=500, conf_level=0.95):
    np.random.seed(seed)
    N = _next_pow2(N)
    U = sobol_sample.sample(problem, N, calc_second_order=calc_second_order)
    Y = model_fun(U)
    extra = {}
    try:
        Si = sobol.analyze(problem, Y,
                           calc_second_order=calc_second_order,
                           conf_level=conf_level,
                           num_resamples=num_resamples,
                           print_to_console=False,
                           keep_resamples=True)
        for k in ["S1_resampled", "ST_resampled"]:
            if k in Si and Si[k] is not None:
                extra[k] = np.asarray(Si[k])
    except TypeError:
        Si = sobol.analyze(problem, Y,
                           calc_second_order=calc_second_order,
                           conf_level=conf_level,
                           num_resamples=num_resamples,
                           print_to_console=False)
    return Si, extra

def summarize_indices(Si, names, boot=None, conf_level=0.95):
    alpha = 1 - conf_level
    lo_q, hi_q = 100 * (alpha/2), 100 * (1 - alpha/2)

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

    S1 = Si["S1"];  ST = Si["ST"]
    S1_conf = Si.get("S1_conf", np.full_like(S1, np.nan))
    ST_conf = Si.get("ST_conf", np.full_like(ST, np.nan))
    S1_res = boot.get("S1_resampled") if boot else None
    ST_res = boot.get("ST_resampled") if boot else None

    S1_med, S1_lo, S1_hi = _summ(S1, resampled=S1_res, conf_half=S1_conf)
    ST_med, ST_lo, ST_hi = _summ(ST, resampled=ST_res, conf_half=ST_conf)

    return pd.DataFrame({
        "parameter": names,
        "S1_median": S1_med, "S1_lo95": S1_lo, "S1_hi95": S1_hi,
        "ST_median": ST_med, "ST_lo95": ST_lo, "ST_hi95": ST_hi
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
        "delta_conf": Si.get("delta_conf", np.full(len(names), np.nan, dtype=float))
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
        "sigma": Si["sigma"]
    })
    out.to_csv(fname, index=False)
    return out

# ============================================================
# 5. Define component models (CF_B, CF_F, N_F)
# ============================================================
def model_CFB(U):
    cB  = 300 + U[:,0]*400
    MGP = q_emp(U[:,1], murA_GP_vec)
    MGN = q_emp(U[:,2], murA_GN_vec)
    fGP = 0.1 + U[:,3]*0.8
    return CFB_fun(cB, MGP, MGN, fGP)

problem_CFB = {"num_vars": 4, "names": ["cB","MGP","MGN","fGP"], "bounds": [[0,1]]*4}

def model_CFF(U):
    cF = 300 + U[:,0]*400
    GF = q_emp(U[:,1], glcN_F_vec)
    return CFF_fun(cF, GF)

problem_CFF = {"num_vars": 2, "names": ["cF","GF"], "bounds": [[0,1]]*2}

def model_fungNecC(U):
    GS = qGS(U[:,0])
    MS = qMS(U[:,1])
    cF = 300 + U[:,2]*400
    GF = q_emp(U[:,3], glcN_F_vec)
    rB = U[:,4]*4.0
    return fungNecC_fun(GS, MS, cF, GF, rB)

problem_fung = {"num_vars": 5, "names": ["GS","MS","cF","GF","rB"], "bounds": [[0,1]]*5}

# ============================================================
# 6. Run analyses (base)
# ============================================================
N = 4096
print("Running Sobol analyses...")
Si_CFB,     B_CFB     = run_sobol(problem_CFB,     model_CFB,      N=N)
Si_CFF,     B_CFF     = run_sobol(problem_CFF,     model_CFF,      N=N)
Si_fung,    B_fung    = run_sobol(problem_fung,    model_fungNecC, N=N)
Si_fnecC5,  B_fnecC5  = run_sobol(problem_fnecC5,  model_fnecC5,   N=N)
Si_fnecC10, B_fnecC10 = run_sobol(problem_fnecC10, model_fnecC10,  N=N)
print("All analyses completed.")

print("Running Delta (N_F, fnecC_agg) and Morris (fnecC_int)...")
Si_delta_fung   = run_delta(problem_fung,   model_fungNecC, N=N)
tab_delta_fung  = export_delta(Si_delta_fung, problem_fung["names"], "sobol_outputs/delta_NF.csv")

Si_delta_fnec5  = run_delta(problem_fnecC5, model_fnecC5,   N=N)
tab_delta_fnec5 = export_delta(Si_delta_fnec5, problem_fnecC5["names"], "sobol_outputs/delta_fnecC_agg.csv")

Si_morris_fnec10  = run_morris(problem_fnecC10, model_fnecC10,
                               N_trajectories=1500, num_levels=8, optimal_trajectories=10)
tab_morris_fnec10 = export_morris(Si_morris_fnec10, problem_fnecC10["names"],
                                  "sobol_outputs/morris_fnecC_int.csv")

# ============================================================
# 7. Export tidy Sobol tables with CIs
# ============================================================
def summarize_and_save(Si, boot, names, out_path):
    tab = summarize_indices(Si, names, boot)
    tab.to_csv(out_path, index=False)
    return tab

tab_CFB     = summarize_and_save(Si_CFB,     B_CFB,     problem_CFB["names"],    f"{OUTDIR}/sobol_CFB_bootstrap.csv")
tab_CFF     = summarize_and_save(Si_CFF,     B_CFF,     problem_CFF["names"],    f"{OUTDIR}/sobol_CFF_bootstrap.csv")
tab_fung    = summarize_and_save(Si_fung,    B_fung,    problem_fung["names"],   f"{OUTDIR}/sobol_NF_bootstrap.csv")
tab_fnecC5  = summarize_and_save(Si_fnecC5,  B_fnecC5,  problem_fnecC5["names"], f"{OUTDIR}/sobol_fnecC_agg_bootstrap.csv")
tab_fnecC10 = summarize_and_save(Si_fnecC10, B_fnecC10, problem_fnecC10["names"],f"{OUTDIR}/sobol_fnecC_int_bootstrap.csv")

# ============================================================
# 8. Plotting utilities
# ============================================================
param_labels = {
    "cB":  r"$\mathit{c}_{B}$", "cF":  r"$\mathit{c}_{F}$",
    "fGP": r"$\mathit{f}_{GP}$", "GF":  r"$\mathit{G}_{F}$",
    "GS":  r"$\mathit{G}_{S}$",  "MGN": r"$\mathit{M}_{GN}$",
    "MGP": r"$\mathit{M}_{GP}$", "MS":  r"$\mathit{M}_{S}$",
    "rB":  r"$\mathit{r}_{B}$",  "S":   r"$\mathit{S}$",
    "CFB": r"$CF_{B}$",          "CFF": r"$CF_{F}$",
}
def pretty_names(names): return [param_labels.get(n, n) for n in names]

def plot_sensitivity(df, title, fname, kind=None):
    cols = set(df.columns.str.lower())
    if kind is None:
        if {"s1_median","st_median"}.issubset(cols): kind = "sobol"
        elif {"delta"}.issubset(cols):               kind = "delta"
        elif {"mu_star","sigma"}.issubset(cols):     kind = "morris"
        else: raise ValueError("Cannot infer plot kind from columns.")

    plt.rcParams.update({
        "mathtext.fontset": "dejavusans",
        "font.size": 10, "axes.labelsize": 10,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
    })

    if kind == "sobol":
        d = (df.set_index("parameter")
               .loc[df.sort_values("ST_median")["parameter"]]
               .reset_index())
        y = np.arange(len(d))
        fig, ax = plt.subplots(figsize=(6.8, 5.0))
        # S1
        x   = d["S1_median"].to_numpy()
        xlo = d["S1_median"] - d["S1_lo95"]
        xhi = d["S1_hi95"]   - d["S1_median"]
        ax.errorbar(x, y, xerr=[xlo, xhi], fmt="o", color=sobol_colors["S1"],
                    label="First-order (median, 95% CI)", capsize=2, markersize=5)
        # ST
        xt   = d["ST_median"].to_numpy()
        xtlo = d["ST_median"] - d["ST_lo95"]
        xthi = d["ST_hi95"]   - d["ST_median"]
        ax.errorbar(xt, y+0.15, xerr=[xtlo, xthi], fmt="s", color=sobol_colors["ST"],
                    label="Total-order (median, 95% CI)", capsize=2, markersize=5)

        ax.set_yticks(y, pretty_names(d["parameter"]))
        ax.set_xlabel("Sobol index")
        ax.set_title(title)
        ax.legend(loc="lower right")
        fig.savefig(fname, dpi=300, bbox_inches="tight")
        plt.close(fig)

    elif kind == "delta":
        d = df.sort_values("delta")
        y = np.arange(len(d))
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        ax.errorbar(d["delta"], y, xerr=d.get("delta_conf", np.zeros_like(y)),
                    fmt="o", capsize=3, label="Δ ± CI")
        ax.set_yticks(y, pretty_names(d["parameter"]))
        ax.set_xlabel("Delta (moment-independent sensitivity)")
        ax.set_title(title)
        ax.legend().remove()
        fig.savefig(fname, dpi=300, bbox_inches="tight")
        plt.close(fig)

    elif kind == "morris":
        d = df.copy()
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        s = 400 * (d["mu_star"] / (d["mu_star"].max() + 1e-12)) + 30
        ax.scatter(d["mu_star"], d["sigma"], s=s, alpha=0.65)
        for _, r in d.iterrows():
            ax.text(r["mu_star"], r["sigma"],
                    param_labels.get(r["parameter"], r["parameter"]),
                    fontsize=8, ha="left", va="center")
        ax.set_xlabel(r"Morris $\mu^{\ast}$ (importance)")
        ax.set_ylabel(r"Morris $\sigma$ (nonlinearity / interactions)")
        ax.set_title(title)
        fig.savefig(fname, dpi=300, bbox_inches="tight")
        plt.close(fig)

# Per-model plots (serve as Supplementary)
plot_sensitivity(tab_CFB, "CF_B — Sobol (median ±95% CI)", f"{OUTDIR}/sobol_CFB_bootstrap.png", kind="sobol")
plot_sensitivity(tab_CFF, "CF_F — Sobol (median ±95% CI)", f"{OUTDIR}/sobol_CFF_bootstrap.png", kind="sobol")
plot_sensitivity(tab_fung, "N_F — Sobol (median ±95% CI)", f"{OUTDIR}/sobol_NF_bootstrap.png", kind="sobol")
plot_sensitivity(tab_CFB, "fₙₑcC (aggregated, r_B profiled) — Sobol", f"{OUTDIR}/sobol_fnecC_agg_bootstrap.png", kind="sobol")
plot_sensitivity(tab_CFB, "fₙₑcC (integrated) — Sobol", f"{OUTDIR}/sobol_fnecC_int_bootstrap.png", kind="sobol")

# ============================================================
# 9. Figure 2 (main): 4 panels; integrated = Top-5 by ST
# ============================================================
from matplotlib.gridspec import GridSpec

def _topk_by_ST(df, k=5):
    d = df.sort_values("ST_median", ascending=False)
    return d.head(k).sort_values("ST_median")

def _panel(ax, df, title):
    d = df.set_index("parameter").loc[df.sort_values("ST_median")["parameter"]].reset_index()
    y = np.arange(len(d))
    x   = d["S1_median"].to_numpy()
    xlo = (d["S1_median"] - d["S1_lo95"]).to_numpy()
    xhi = (d["S1_hi95"]   - d["S1_median"]).to_numpy()
    ax.errorbar(x, y, xerr=[xlo, xhi], fmt="o", color=sobol_colors["S1"], capsize=2, markersize=4, label="S₁")
    xt   = d["ST_median"].to_numpy()
    xtlo = (d["ST_median"] - d["ST_lo95"]).to_numpy()
    xthi = (d["ST_hi95"]   - d["ST_median"]).to_numpy()
    ax.errorbar(xt, y+0.15, xerr=[xtlo, xthi], fmt="s", color=sobol_colors["ST"], capsize=2, markersize=4, label="S_T")
    ax.set_yticks(y, pretty_names(d["parameter"]))
    ax.set_xlabel("Sobol index")
    ax.set_title(title, pad=4)
    ax.grid(axis="x", alpha=0.15)
    return ax

def figure2_panels(tab_CFB, tab_NF, tab_fnec_agg, tab_fnec_int,
                   outpng=f"{OUTDIR}/Figure2_CFpipeline_Sobol.png"):
    plt.rcParams.update({"mathtext.fontset": "dejavusans", "font.size": 9})
    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)  # FIX: constrained layout
    gs = GridSpec(2, 2, figure=fig)

    axA = fig.add_subplot(gs[0,0]); _panel(axA, tab_CFB, "A  CF_B")
    axB = fig.add_subplot(gs[0,1]); _panel(axB, tab_NF,  "B  N_F")
    axC = fig.add_subplot(gs[1,0]); _panel(axC, tab_fnec_agg, "C  fₙₑcC (aggregated)")
    axD = fig.add_subplot(gs[1,1]); _panel(axD, _topk_by_ST(tab_fnec_int, k=5), "D  fₙₑcC (integrated, Top-5)")

    # Figure-level legend
    handles, labels = axD.get_legend_handles_labels()
    fig.legend(handles, ["S₁ (median, 95% CI)", "S_T (median, 95% CI)"],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)

# Build main Figure 2
figure2_panels(
    tab_CFB=tab_CFB, tab_NF=tab_fung,
    tab_fnec_agg=tab_fnecC5, tab_fnec_int=tab_fnecC10,
    outpng="sobol_outputs/Figure2_CFpipeline_Sobol.png"
)
tab_fnecC10_top5 = _topk_by_ST(tab_fnecC10, k=5)
tab_fnecC10_top5.to_csv("sobol_outputs/sobol_fnecC_int_top5_for_Fig2.csv", index=False)

# ============================================================
# Save raw model outputs for verification
# ============================================================

# Reuse the same sampling inputs from your Sobol analyses
N = 4096
U5  = sobol_sample.sample(problem_fnecC5, N, calc_second_order=False)
U10 = sobol_sample.sample(problem_fnecC10, N, calc_second_order=False)

# Evaluate both models
Y5  = model_fnecC5(U5)
Y10 = model_fnecC10(U10)

# Save as CSVs for verification / plotting
import pandas as pd
pd.DataFrame({
    "fnecC_agg": Y5
}).to_csv("sobol_outputs/fnecC5_raw.csv", index=False)

pd.DataFrame({
    "fnecC_int": Y10
}).to_csv("sobol_outputs/fnecC10_raw.csv", index=False)

# ============================================================
# Save raw model outputs for CF_B, CF_F, and N_F (fungal term)
# ============================================================

print("Saving CF_B, CF_F, and N_F raw outputs...")

# 1. CF_B (Eq. 1) — 4 inputs: cB, MGP, MGN, fGP
N = 4096
U_b = sobol_sample.sample(problem_CFB, N, calc_second_order=False)
Y_CFB = model_CFB(U_b)
pd.DataFrame({"CF_B": Y_CFB}).to_csv("sobol_outputs/CFB_raw.csv", index=False)

# 2. CF_F (Eq. 2a) — 2 inputs: cF, GF
U_f = sobol_sample.sample(problem_CFF, N, calc_second_order=False)
Y_CFF = model_CFF(U_f)
pd.DataFrame({"CF_F": Y_CFF}).to_csv("sobol_outputs/CFF_raw.csv", index=False)

# 3. N_F (Eq. 2b) — 5 inputs: GS, MS, cF, GF, rB
U_nf = sobol_sample.sample(problem_fung, N, calc_second_order=False)
Y_NF = model_fungNecC(U_nf)
pd.DataFrame({"N_F": Y_NF}).to_csv("sobol_outputs/NF_raw.csv", index=False)

print("All raw model outputs saved to sobol_outputs/")

# ============================================================
# 10. rB profile sweep for aggregated model + stacked figure
# ============================================================
def run_sobol_fnecC5_profiles(rB_list=(0, 1.0, 1.16, 2.0, 3.0, 4.0),
                              N=4096, seed=42, conf_level=0.95, outdir="sobol_outputs"):
    os.makedirs(outdir, exist_ok=True)
    tabs = []
    for rB_val in rB_list:
        Si, boot = run_sobol(problem_fnecC5, lambda U: model_fnecC5(U, rB_const=rB_val),
                             N=N, seed=seed, calc_second_order=False,
                             num_resamples=500, conf_level=conf_level)
        tab = summarize_indices(Si, problem_fnecC5["names"], boot, conf_level=conf_level)
        tab.insert(0, "rB", rB_val)
        tab.to_csv(os.path.join(outdir, f"sobol_fnecC_agg_rB{str(rB_val).replace('.','p')}.csv"), index=False)
        tabs.append(tab)
    comb = pd.concat(tabs, ignore_index=True)
    comb.to_csv(os.path.join(outdir, "sobol_fnecC_agg_profiles_long.csv"), index=False)
    return comb

def figure2C_profiles_fnecC_agg(long_tab, rB_list=None,
                                outpng="sobol_outputs/Figure2C_fnecC_agg_profiles.png"):
    from matplotlib.gridspec import GridSpec
    if rB_list is None:
        rB_list = sorted(long_tab["rB"].unique(), key=lambda x: float(x))
    plt.rcParams.update({"mathtext.fontset": "dejavusans", "font.size": 9})
    n = len(rB_list)
    height = 1.8 * n + 1.0
    fig = plt.figure(figsize=(7.2, height), constrained_layout=True)  # FIX: constrained layout
    gs = GridSpec(n, 1, figure=fig)

    for i, rB_val in enumerate(rB_list):
        ax = fig.add_subplot(gs[i, 0])
        df = long_tab.loc[long_tab["rB"] == rB_val,
                          ["parameter","S1_median","S1_lo95","S1_hi95","ST_median","ST_lo95","ST_hi95"]].copy()
        d = df.set_index("parameter").loc[df.sort_values("ST_median")["parameter"]].reset_index()
        y = np.arange(len(d))

        # S1
        x   = d["S1_median"].to_numpy()
        xlo = (d["S1_median"] - d["S1_lo95"]).to_numpy()
        xhi = (d["S1_hi95"]   - d["S1_median"]).to_numpy()
        ax.errorbar(x, y, xerr=[xlo, xhi], fmt="o", color=sobol_colors["S1"], capsize=2, markersize=4, label="S₁")

        # ST
        xt   = d["ST_median"].to_numpy()
        xtlo = (d["ST_median"] - d["ST_lo95"]).to_numpy()
        xthi = (d["ST_hi95"]   - d["ST_median"]).to_numpy()
        ax.errorbar(xt, y + 0.15, xerr=[xtlo, xthi], fmt="s", color=sobol_colors["ST"], capsize=2, markersize=4, label="S_T")

        ax.set_yticks(y, pretty_names(d["parameter"]))
        ax.set_xlabel("Sobol index")
        ax.set_title(f"fₙₑcC (aggregated), profiled r_B = {rB_val:g}", pad=2)
        ax.grid(axis="x", alpha=0.15)
        if i == 0:
            ax.legend(loc="lower right", frameon=False, fontsize=8)

    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---- RUN the rB profile sweep + stacked figure ----
rB_profiles = (0.0, 1.0, 1.16, 2.0, 3.0, 4.0)
agg_profiles_long = run_sobol_fnecC5_profiles(rB_list=rB_profiles, N=4096, seed=42)
figure2C_profiles_fnecC_agg(agg_profiles_long, rB_list=rB_profiles,
                            outpng="sobol_outputs/Figure2C_fnecC_agg_profiles.png")
# ============================================================
# 11. rB significance: ΔNF, ΔfnecC, truncation, and R distribution
# ============================================================
import numpy as np
import pandas as pd

def sample_CFF(n=20000, seed=1234):
    rng = np.random.default_rng(seed)
    u = np.clip(rng.random(n), 1e-12, 1 - 1e-12)
    try:
        return np.quantile(CFF_emp, u, method="linear")
    except TypeError:
        return np.quantile(CFF_emp, u, interpolation="linear")

def NF_with_rB(GS, MS, CFF, rB):
    """Fungal necromass C (mg C g-1 soil) with bacterial correction."""
    GS_corr = np.maximum(GS - rB * MS, 0.0)
    return GS_corr * CFF

def fnecC_with_rB(GS, MS, S, CFF, CFB, rB):
    """Total necromass fraction of SOC (%) using aggregated form and given CFB/CFF."""
    S = np.maximum(S, EPS)
    GS_corr = np.maximum(GS - rB * MS, 0.0)
    return 100.0 * (MS * CFB + GS_corr * CFF) / S

# --- replace your bootstrap_ci with a streaming subsample bootstrap ---
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
    import os
    os.makedirs(outdir, exist_ok=True)

    # soil vectors
    GS = soils["GS"].to_numpy(dtype=float)
    MS = soils["MS"].to_numpy(dtype=float)
    S  = soils["S"].to_numpy(dtype=float)

    # empirical CF samples
    CFF_s = sample_CFF(n=n_cff, seed=seed)
    rng = np.random.default_rng(seed+1)
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

        # build the row dict in one go (no bare strings!)
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
            "equiv_NF":    (dNF_lo   >= -eps_C)  and (dNF_hi   <= eps_C),
            "equiv_fnecC": (dfnec_lo >= -eps_pct) and (dfnec_hi <= eps_pct),
        }
        rows.append(row)

        # optional: keep full R distribution (for supplementary figure)
        R_summaries.append(pd.DataFrame({"rB": rB, "R_percent": R_vals}))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(outdir, "rB_significance_summary.csv"), index=False)

    R_all = pd.concat(R_summaries, ignore_index=True)
    R_all.to_csv(os.path.join(outdir, "rB_R_distribution_long.csv"), index=False)

    return out
# ---- run and save
rB_sig = rB_significance_analysis(
    rB_values=(0.0, 1.0, 1.16, 2.0, 3.0, 4.0),
    n_cff=2000,
    seed=1234,
    outdir=OUTDIR,
    eps_C=0.001,
    eps_pct=0.5
)
print("rB significance summary:\n", rB_sig.head())
# ============================================================
# 12. Plot: rB effects — ΔNF, ΔfnecC, and R vs rB
# ============================================================
import matplotlib.pyplot as plt

def plot_rB_effects(summary_csv="sobol_outputs/rB_significance_summary.csv",
                    outpng="sobol_outputs/Fig4C_rB_Delta_and_R.png"):
    df = pd.read_csv(summary_csv)

    # Order by rB
    df = df.sort_values("rB").reset_index(drop=True)

    # Figure
    plt.rcParams.update({
        "mathtext.fontset": "dejavusans",
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    })
    fig, axes = plt.subplots(2, 1, figsize=(6.2, 5.2), constrained_layout=True)

    # ---------- Top: ΔNF and ΔfnecC ----------
    ax = axes[0]
    # ΔNF (mg C g-1)
    ax.errorbar(
        df["rB"], df["dNF_median_mgC_g"],
        yerr=[df["dNF_median_mgC_g"] - df["dNF_lo95"],
              df["dNF_hi95"] - df["dNF_median_mgC_g"]],
        fmt="o-", markersize=6, color="black", mfc="black", mec="black",
        capsize=3, label=r"$\Delta N_F$ (mg C g$^{-1}$)", zorder=3
)

    # twin axis for ΔfnecC (%-SOC)
    ax2 = ax.twinx()
    ax2.errorbar(df["rB"], df["dfnec_median_pct"],
                 yerr=[df["dfnec_median_pct"] - df["dfnec_lo95"],
                       df["dfnec_hi95"] - df["dfnec_median_pct"]],
                 fmt="s--", capsize=3, label=r"$\Delta f_{\mathrm{necC}}$ (pp)")

    ax.axhline(0, color="0.7", linestyle=":", linewidth=1)
    ax.set_xlabel(None)
    ax.set_ylabel(r"$\Delta N_F$ (mg C g$^{-1}$)")
    ax2.set_ylabel(r"$\Delta f_{\mathrm{necC}}$ (percentage points)")

    # single legend combining both axes
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left", frameon=False)

    # ---------- Bottom: R (median ± 95% CI) ----------
    ax = axes[1]
    ax.errorbar(df["rB"], df["R_median_pct"],
                yerr=[df["R_median_pct"] - df["R_lo95"],
                      df["R_hi95"] - df["R_median_pct"]],
                fmt="o-", capsize=3, label=r"$R$ median (95% CI)")

    # annotate truncation proportion
    for x, p in zip(df["rB"], df["truncation_prop"]):
        if p > 0:
            ax.text(x, 0, f"trunc={p*100:.1f}%", ha="center", va="bottom", fontsize=8)

    ax.set_xlabel(r"$r_B$")
    ax.set_ylabel(r"$R$ (% of fungal GlcN subtracted)")
    ax.axhline(0, color="0.7", linestyle=":", linewidth=1)
    ax.legend(loc="upper left", frameon=False)

    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)

# create the panel
plot_rB_effects()
