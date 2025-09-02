import os, sys, argparse, pickle, multiprocessing as mp
from functools import partial

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from tqdm.auto import tqdm


# ---------- Small utilities ----------

def zscore_series(x: pd.Series):
    x = pd.to_numeric(x, errors="coerce")
    mu = x.mean()
    sd = x.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return None
    return (x - mu) / sd


def ks_stat_uniform(pvals):
    p = np.asarray(pvals, dtype=float)
    m = np.isfinite(p)
    p = p[m]
    if p.size == 0:
        return np.nan
    return stats.kstest(p, 'uniform').statistic


def signed_logp_vec(pvals, effects):
    """Return signed -log10(p) with sign from effects; NaNs filtered out."""
    pvals = np.asarray(pvals, dtype=float)
    effects = np.asarray(effects, dtype=float)
    m = np.isfinite(pvals) & np.isfinite(effects)
    if not np.any(m):
        return np.array([]), m
    s = -np.log10(np.clip(pvals[m], 1e-300, 1.0)) * np.sign(effects[m])
    return s, m


def lambda_gc_from_p(p):
    # two-sided z via p; chi2=z^2; lambda = median(chi2)/median(chi2_1)
    p = np.asarray(p, dtype=float)
    m = np.isfinite(p)
    p = p[m]
    if p.size == 0:
        return np.nan
    z = stats.norm.isf(p / 2.0)  # two-sided
    chi2 = z ** 2
    return np.median(chi2) / 0.455936  # median(chi2_1)


# ---------- Core per-gene fit on current dataframe ----------

def fit_full_and_contrast(df2, use_shet=True, robust=None):
    """
    df2 columns required:
      single-predictor: post_mean, perturb_beta_cond1, (optional) shet
      two-predictor   : post_mean, perturb_beta_cond1, perturb_beta_cond2, (optional) shet
    Returns: (beta_or_beta_sum, p_value, n, k, R2)
    """
    y = pd.to_numeric(df2["post_mean"], errors="coerce")
    y_z = zscore_series(y)
    if y_z is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    has_c2 = "perturb_beta_cond2" in df2.columns

    x1 = zscore_series(df2["perturb_beta_cond1"]) if "perturb_beta_cond1" in df2.columns else None
    x2 = zscore_series(df2["perturb_beta_cond2"]) if has_c2 else None

    if x1 is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    X = {}
    if has_c2 and x2 is not None:
        # Two-predictor mode
        X["perturb_beta_cond1"] = x1
        X["perturb_beta_cond2"] = x2
        test_sum = True
    else:
        # Single-predictor mode
        X["perturb_beta"] = x1
        test_sum = False

    if use_shet and ("shet" in df2.columns):
        svec = pd.to_numeric(df2["shet"], errors="coerce")
        if svec.nunique(dropna=True) > 1:
            smean = svec.mean()
            sstd = svec.std(ddof=0)
            X["shet"] = (svec - smean) / sstd if (sstd and np.isfinite(sstd) and sstd > 0) else svec

    X = pd.DataFrame(X)
    X = sm.add_constant(X, has_constant="add")

    try:
        res = sm.OLS(y_z, X, missing="drop").fit() if robust is None else sm.OLS(y_z, X, missing="drop").fit(cov_type=robust)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    V = res.cov_params()
    n = int(res.nobs)
    k = int(X.shape[1] - 1)
    R2 = float(res.rsquared) if res.rsquared is not None else np.nan

    if not test_sum:
        # Single predictor: report β and its p
        beta = float(res.params.get("perturb_beta", np.nan))
        if "perturb_beta" in V.index:
            var = float(V.loc["perturb_beta", "perturb_beta"])
            if np.isfinite(var) and var < 0 and var > -1e-12:
                var = 0.0
            se = float(np.sqrt(var)) if (np.isfinite(var) and var >= 0) else np.nan
        else:
            se = np.nan
        if np.isfinite(beta) and np.isfinite(se) and se > 0 and np.isfinite(res.df_resid) and res.df_resid > 0:
            t = beta / se
            p = float(2 * stats.t.sf(abs(t), df=int(res.df_resid)))
        else:
            p = np.nan
        return beta, p, n, k, R2

    # Two predictors: test β1 + β2
    beta_c1 = float(res.params.get("perturb_beta_cond1", np.nan))
    beta_c2 = float(res.params.get("perturb_beta_cond2", np.nan))
    beta_sum = (beta_c1 + beta_c2) if np.isfinite(beta_c1) and np.isfinite(beta_c2) else np.nan

    if {"perturb_beta_cond1", "perturb_beta_cond2"}.issubset(V.index):
        v11 = float(V.loc["perturb_beta_cond1", "perturb_beta_cond1"])
        v22 = float(V.loc["perturb_beta_cond2", "perturb_beta_cond2"])
        v12 = float(V.loc["perturb_beta_cond1", "perturb_beta_cond2"])
        var_sum = v11 + v22 + 2.0 * v12
        if np.isfinite(var_sum) and var_sum < 0 and var_sum > -1e-12:
            var_sum = 0.0
        se_sum = float(np.sqrt(var_sum)) if (np.isfinite(var_sum) and var_sum >= 0) else np.nan
    else:
        se_sum = np.nan

    if np.isfinite(beta_sum) and np.isfinite(se_sum) and se_sum > 0 and np.isfinite(res.df_resid) and res.df_resid > 0:
        t_sum = beta_sum / se_sum
        p_sum = float(2 * stats.t.sf(abs(t_sum), df=int(res.df_resid)))
    else:
        p_sum = np.nan
    return beta_sum, p_sum, n, k, R2


# ---------- Freedman–Lane residual permutation processing for one response gene ----------
# Note: Uses the shared permutation from the current iteration

def freedman_lane_one(df2, base_order_invpos, perm_invpos, use_shet=True, robust=None):
    """
    Process one response gene using the shared permutation from the current iteration.
    
    df2: long table for a single response gene (rows=targets)
    base_order_invpos: dict ensg -> 0..M-1 positions in the common base target list
    perm_invpos: array(len=M) mapping base OLD index -> permuted NEW rank (inverse permutation)
                 This is the SAME permutation used for ALL response genes in this iteration
    Returns function that yields (beta_or_beta_sum, p_value) for this gene using the shared permutation
    """
    # Reduced model: y ~ const + (shet if varying); NO perturb betas here
    y = pd.to_numeric(df2["post_mean"], errors="coerce")
    y_z = zscore_series(y)
    if y_z is None:
        return None

    Z = pd.DataFrame({"const": 1.0}, index=df2.index)
    use_shet_here = False
    if use_shet and ("shet" in df2.columns):
        svec = pd.to_numeric(df2["shet"], errors="coerce")
        if svec.nunique(dropna=True) > 1:
            Z["shet"] = svec
            use_shet_here = True

    # Fit reduced
    try:
        res_red = sm.OLS(y_z, Z, missing="drop").fit()
    except Exception:
        return None

    yhat_red = pd.Series(res_red.fittedvalues, index=Z.index).reindex(df2.index)
    r = (y_z - yhat_red).astype(float).values  # residuals in df2 row order

    # Map df2 rows to base indices (if any missing, fall back to within-df2 permutation)
    ensg_subset = df2["ensg"].values
    try:
        base_idx = np.array([base_order_invpos[e] for e in ensg_subset], dtype=int)
        perm_ranks = perm_invpos[base_idx]
    except KeyError:
        base_idx = np.arange(len(df2), dtype=int)
        perm_ranks = np.arange(len(df2), dtype=int)

    perm_order = np.argsort(perm_ranks)

    def refit_on_permuted_residuals():
        r_pi = r[perm_order]
        y_pi = yhat_red.values + r_pi
        df2_pi = df2.copy()
        df2_pi["post_mean"] = y_pi
        beta_sum, p_sum, *_ = fit_full_and_contrast(df2_pi, use_shet=use_shet_here, robust=robust)
        return beta_sum, p_sum

    return refit_on_permuted_residuals


# ---------- Pair design ----------

def build_pair_design(dfA: pd.DataFrame, dfB: pd.DataFrame = None):
    """
    If dfB is None, single-condition mode: returns dfA restricted to itself (columns unchanged).
    Otherwise, restrict both to common targets and return pair.
    """
    if dfB is None:
        return dfA.copy(), None
    common_targets = sorted(set(dfA.columns).intersection(dfB.columns))
    A2 = dfA.loc[:, common_targets]
    B2 = dfB.loc[:, common_targets]
    return A2, B2


# ---------- Single iteration analysis ----------

def run_single_iteration(
    cond_mats,
    ensg2sym,
    sym2ensg,
    lof,
    shet,
    condA, condB=None, seed=1, robust=None, min_targets=20, use_shet=True,
    iteration_idx=0
):
    """
    Run a single permutation iteration and return results.
    This is designed to be run independently for each iteration.
    """
    rng = np.random.RandomState(seed + iteration_idx)  # Different seed for each iteration

    # Build design(s)
    if condB is None:
        dfA2, dfB2 = build_pair_design(cond_mats[condA], None)
        common_resp = dfA2.index
        dfA2 = dfA2.loc[common_resp, :]
    else:
        dfA2, dfB2 = build_pair_design(cond_mats[condA], cond_mats[condB])
        common_resp = dfA2.index.intersection(dfB2.index)
        dfA2 = dfA2.loc[common_resp, :]
        dfB2 = dfB2.loc[common_resp, :]

    # Base target universe (for consistent permutations across genes)
    base_targets_sym = list(dfA2.columns)
    base_targets_ensg = pd.Series(base_targets_sym).map(sym2ensg).values
    keep = pd.notna(base_targets_ensg)
    base_targets_ensg = base_targets_ensg[keep]
    base_order_invpos = {e: i for i, e in enumerate(base_targets_ensg)}
    M = len(base_targets_ensg)

    # Generate single permutation
    perm = rng.permutation(M)                 # perm[new] = old
    inv = np.empty(M, dtype=int)
    inv[perm] = np.arange(M, dtype=int)      # inv[old]=new
    perm_invpos = inv

    # --- Observed fits ---
    rows = []
    per_gene_df2_cache = {}
    label = f"{condA}+{condB}" if condB is not None else condA

    for resp_ensg in common_resp:
        if condB is None:
            # Single predictor
            tmp = pd.DataFrame({
                "target_symbol": dfA2.columns,
                "perturb_beta_cond1": dfA2.loc[resp_ensg, :].values
            })
        else:
            # Two predictors
            tmp = pd.DataFrame({
                "target_symbol": dfA2.columns,
                "perturb_beta_cond1": dfA2.loc[resp_ensg, :].values,
                "perturb_beta_cond2": dfB2.loc[resp_ensg, :].values
            })

        tmp["ensg"] = tmp["target_symbol"].map(sym2ensg)
        tmp = tmp.dropna(subset=["ensg"])
        tmp = tmp[tmp["ensg"] != resp_ensg]

        df2 = tmp.merge(lof, on="ensg", how="inner")
        if shet is not None:
            df2 = df2.merge(shet, on="ensg", how="inner")
        if df2.shape[0] <= min_targets:
            continue

        beta_sum, p_sum, n, k, R2 = fit_full_and_contrast(df2, use_shet=use_shet, robust=robust)
        rows.append({
            "gene": ensg2sym.get(resp_ensg, resp_ensg),
            "ensg": resp_ensg,
            "beta": beta_sum,
            "p_beta": p_sum,
            "n": n, "k": k, "R2": R2
        })
        per_gene_df2_cache[resp_ensg] = df2

    obs_df = pd.DataFrame(rows).reset_index(drop=True)
    resp_list = list(obs_df["ensg"].values)
    obs_p = obs_df["p_beta"].values
    obs_beta = obs_df["beta"].values

    # --- Build per-gene reduced-model refitters (closures) ---
    refitters = []
    for resp in resp_list:
        df2 = per_gene_df2_cache.get(resp, None)
        refitters.append((resp, df2) if df2 is not None else None)

    # --- Single permutation ---
    p_vec = np.full(len(resp_list), np.nan, dtype=float)
    b_vec = np.full(len(resp_list), np.nan, dtype=float)

    for i, item in enumerate(refitters):
        if item is None:
            continue
        _, df2 = item
        refit = freedman_lane_one(df2, base_order_invpos, perm_invpos, use_shet=use_shet, robust=robust)
        if refit is None:
            continue
        beta_pi, p_pi = refit()
        p_vec[i] = p_pi
        b_vec[i] = beta_pi

    # Return results for this iteration
    return {
        "iteration": iteration_idx,
        "label": label,
        "obs_df": obs_df,
        "perm_pvals": p_vec,
        "perm_betas": b_vec,
        "base_order_invpos": base_order_invpos,
        "M": M
    }


# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond_mats_file", required=True)
    ap.add_argument("--ensg2sym_file", required=True)
    ap.add_argument("--sym2ensg_file", required=True)
    ap.add_argument("--lof_file", required=True)
    ap.add_argument("--shet_file", required=True)
    ap.add_argument("--conds", required=True, help="Comma-separated list; entries may be 'A' or 'A+B'")
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--min_targets", required=True)
    ap.add_argument("--iteration_idx", type=int, required=True, help="Current iteration index (0-based)")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Loading cond_mats from {args.cond_mats_file}")
    with open(args.cond_mats_file, "rb") as f:
        cond_mats = pickle.load(f)
    print(f"Loading ensg2sym from {args.ensg2sym_file}")
    with open(args.ensg2sym_file, "rb") as f:
        ensg2sym = pickle.load(f)
    print(f"Loading sym2ensg from {args.sym2ensg_file}")
    with open(args.sym2ensg_file, "rb") as f:
        sym2ensg = pickle.load(f)
    print(f"Loading lof from {args.lof_file}")
    lof = pd.read_csv(args.lof_file, sep="\t")
    print(f"Loading shet from {args.shet_file}")
    shet = pd.read_csv(args.shet_file, sep="\t")

    # Ensure lof contains 'ensg' and post_mean; merge shet if not already included
    if "shet" not in lof.columns and "shet" in shet.columns and "ensg" in shet.columns:
        lof = lof.merge(shet[["ensg", "shet"]], on="ensg", how="left")

    conds = [c.strip() for c in args.conds.split(",") if c.strip()]
    min_targets = int(args.min_targets)
    seed = int(args.seed)
    iteration_idx = int(args.iteration_idx)

    # Run single iteration for each condition
    results = {}
    for cond in conds:
        if "+" in cond:
            condA, condB = [s.strip() for s in cond.split("+", 1)]
        else:
            condA, condB = cond, None  # single-predictor mode

        print(f"Running iteration {iteration_idx} for condition: {cond}")
        result = run_single_iteration(
            cond_mats,
            ensg2sym,
            sym2ensg,
            lof,
            shet,
            condA, condB,
            seed=seed, robust=None,
            min_targets=min_targets, use_shet=True,
            iteration_idx=iteration_idx
        )
        results[result["label"]] = result

    # Save results for this iteration
    output_file = os.path.join(args.output_dir, f"iteration_{iteration_idx:04d}.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    print(f"Saved iteration {iteration_idx} results to {output_file}")


if __name__ == "__main__":
    main()
