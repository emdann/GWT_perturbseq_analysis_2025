import os, sys, argparse, pickle, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


# ---------- Small utilities ----------

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


# ---------- Integration functions ----------

def load_iteration_results(output_dir, expected_iterations):
    """
    Load all iteration results from the output directory.
    Returns: dict {label -> (obs_df, perm_pvals, perm_betas)}
    """
    print(f"Loading iteration results from {output_dir}")
    
    # Find all iteration files
    iteration_files = glob.glob(os.path.join(output_dir, "iteration_*.pkl"))
    iteration_files.sort()
    
    print(f"Found {len(iteration_files)} iteration files")
    
    if len(iteration_files) != expected_iterations:
        print(f"Warning: Expected {expected_iterations} iterations but found {len(iteration_files)} files")
    
    # Load first file to get structure
    with open(iteration_files[0], "rb") as f:
        first_result = pickle.load(f)
    
    # Initialize results structure
    integrated_results = {}
    for label in first_result.keys():
        # Get observed data from first iteration (should be same for all)
        obs_df = first_result[label]["obs_df"]
        
        # Initialize arrays for permutation results
        n_genes = len(obs_df)
        n_iterations = len(iteration_files)
        
        perm_pvals = np.full((n_iterations, n_genes), np.nan, dtype=float)
        perm_betas = np.full((n_iterations, n_genes), np.nan, dtype=float)
        
        integrated_results[label] = {
            "obs_df": obs_df,
            "perm_pvals": perm_pvals,
            "perm_betas": perm_betas
        }
    
    # Load all iteration results
    for i, iteration_file in enumerate(iteration_files):
        print(f"Loading iteration {i+1}/{len(iteration_files)}: {os.path.basename(iteration_file)}")
        
        with open(iteration_file, "rb") as f:
            iteration_result = pickle.load(f)
        
        for label, result in iteration_result.items():
            if label in integrated_results:
                integrated_results[label]["perm_pvals"][i, :] = result["perm_pvals"]
                integrated_results[label]["perm_betas"][i, :] = result["perm_betas"]
    
    return integrated_results


def calculate_global_calibration(obs_df, perm_pvals, perm_betas):
    """
    Calculate global calibration metrics from observed and permutation results.
    """
    obs_p = obs_df["p_beta"].values
    obs_beta = obs_df["beta"].values
    
    # KS test for uniformity
    ks_obs = ks_stat_uniform(obs_p)
    ks_perm = np.array([ks_stat_uniform(pp) for pp in perm_pvals])
    global_p = (1 + np.sum(ks_perm >= ks_obs)) / (len(ks_perm) + 1)
    
    # Lambda GC
    lam_obs = lambda_gc_from_p(obs_p)
    lam_perm = np.array([lambda_gc_from_p(pp) for pp in perm_pvals if np.isfinite(pp).any()])
    lam_ratio = lam_obs / np.nanmedian(lam_perm) if (np.isfinite(lam_obs) and lam_perm.size) else np.nan
    
    return {
        "KS_stat_observed": ks_obs,
        "KS_global_p": global_p,
        "lambda_GC_observed": lam_obs,
        "lambda_GC_perm_median": np.nanmedian(lam_perm) if lam_perm.size else np.nan,
        "lambda_GC_ratio": lam_ratio
    }


def qq_bands_from_perms(obs_p, obs_beta, perm_pvals, perm_betas, alpha=0.05):
    """
    Build 1-alpha permutation bands for signed -log10(p).
    Returns y_obs_sorted, lower, upper; all length K (number of finite obs points).
    """
    y_obs, mask = signed_logp_vec(obs_p, obs_beta)
    if y_obs.size == 0:
        return np.array([]), np.array([]), np.array([])

    y_obs_sorted = np.sort(y_obs)
    K = len(y_obs_sorted)

    # Build permutation signed -log10 p's (sorted)
    mats = []
    for p_vec, b_vec in zip(perm_pvals, perm_betas):
        if not np.isfinite(p_vec).any():
            continue
        y, m = signed_logp_vec(p_vec, b_vec)
        if y.size == 0:
            continue
        y = np.sort(y)
        if len(y) >= K:
            mats.append(y[:K])
        else:
            pad = np.full(K - len(y), y[-1])
            mats.append(np.concatenate([y, pad]))

    if len(mats) == 0:
        lower = upper = np.full(K, np.nan)
    else:
        perm_mat = np.vstack(mats)
        lower = np.quantile(perm_mat, alpha / 2, axis=0)
        upper = np.quantile(perm_mat, 1 - alpha / 2, axis=0)

    return y_obs_sorted, lower, upper


def qq_signed_logp(
    pvals, effects, ax, label=None, genes=None, top_n=0,
    gene_set=None, gene_set_color=None, base_color=None
):
    """
    Scatter a signed QQ plot and return (ax, top_genes_list, theor_x).
    NaN p/effects are dropped consistently.
    """
    pvals = np.asarray(pvals, dtype=float)
    effects = np.asarray(effects, dtype=float)
    mask = np.isfinite(pvals) & np.isfinite(effects)
    if genes is not None:
        genes = np.asarray(genes)[mask] if hasattr(genes, "__len__") else None
    pvals = pvals[mask]
    effects = effects[mask]

    if pvals.size == 0:
        return ax, [], np.array([])

    signed_logp = -np.log10(np.clip(pvals, 1e-300, 1.0)) * np.sign(effects)
    order = np.argsort(signed_logp)
    obs = signed_logp[order]
    n = len(obs)
    q = (np.arange(1, n + 1) - 0.5) / n

    theor = np.empty(n)
    lower = q <= 0.5
    theor[lower] = np.log10(2 * q[lower])
    theor[~lower] = -np.log10(2 * (1 - q[~lower]))

    scatter_color = base_color if base_color is not None else None
    ax.scatter(theor, obs, s=15, alpha=0.7, label=label, color=scatter_color)

    # Optional highlight of a gene set
    if gene_set is not None and genes is not None:
        gene_set = set(gene_set)
        name2idx = {g: i for i, g in enumerate(genes)}
        present = [name2idx[g] for g in gene_set if g in name2idx]
        if present:
            inv = np.empty_like(order)
            inv[order] = np.arange(n)
            sorted_pos = [int(inv[idx]) for idx in present]
            ax.scatter(theor[sorted_pos], obs[sorted_pos], s=50, alpha=0.8,
                       color=gene_set_color, edgecolors='black', linewidth=1,
                       label=f'{label}_gene_set' if label else 'gene_set')

    top_genes_list = []
    if genes is not None and top_n > 0:
        top_idxs = np.argsort(np.abs(signed_logp))[-top_n:]
        top_genes_list = [(genes[idx], signed_logp[idx], pvals[idx], effects[idx]) for idx in top_idxs]

    return ax, top_genes_list, theor


def parse_color_map(s):
    """
    Parse --base_color_map like:
      Stim8hr:#1f77b4,Stim48hr:#ff7f0e,Stim8hr+Stim48hr:#2ca02c
    Returns dict {cond_string: color}
    """
    cm = {}
    if not s:
        return cm
    for part in s.split(","):
        if ":" not in part:
            continue
        k, v = part.split(":", 1)
        cm[k.strip()] = v.strip()
    return cm


# ---------- Main integration ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output_dir", required=True, help="Directory containing iteration results")
    ap.add_argument("--base_color_map", required=True, help='Comma-separated "cond:color" entries')
    ap.add_argument("--permutations", type=int, required=True, help="Total number of permutations run")
    ap.add_argument("--final_output_dir", required=True, help="Directory for final integrated results")
    args = ap.parse_args()

    os.makedirs(args.final_output_dir, exist_ok=True)

    # Load all iteration results
    integrated_results = load_iteration_results(args.output_dir, args.permutations)
    
    # Calculate global calibration for each condition
    gc_rows = []
    for label, result in integrated_results.items():
        print(f"\nProcessing {label}...")
        
        obs_df = result["obs_df"]
        perm_pvals = result["perm_pvals"]
        perm_betas = result["perm_betas"]
        
        # Calculate global calibration
        global_calib = calculate_global_calibration(obs_df, perm_pvals, perm_betas)
        
        print(f"  KS global p = {global_calib['KS_global_p']:.3g}")
        print(f"  λGC = {global_calib['lambda_GC_observed']:.3f} (perm median {global_calib['lambda_GC_perm_median']:.3f})")
        
        # Save observed results
        safe_label = label.replace("+", "_")
        obs_csv = os.path.join(args.final_output_dir, f"{safe_label}_observed.csv")
        obs_df.to_csv(obs_csv, index=False)
        print(f"  wrote {obs_csv}")
        
        # Store for global summary
        row = {"cond": label, **global_calib}
        gc_rows.append(row)
        
        # Store in integrated results for plotting
        integrated_results[label]["global_calib"] = global_calib
    
    # Save global calibration summary
    gc_df = pd.DataFrame(gc_rows)
    gc_csv = os.path.join(args.final_output_dir, "freedman_lane_global_calibration.csv")
    gc_df.to_csv(gc_csv, index=False)
    print(f"\nWrote {gc_csv}")
    
    # Save full integrated results
    results_pkl = os.path.join(args.final_output_dir, "freedman_lane_integrated_results.pkl")
    with open(results_pkl, "wb") as f:
        pickle.dump(integrated_results, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Wrote {results_pkl}")
    
    # Create QQ plot with permutation bands
    print("\nCreating QQ plot...")
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    color_map = parse_color_map(args.base_color_map)
    
    for label, result in integrated_results.items():
        obs_df = result["obs_df"]
        perm_pvals = result["perm_pvals"]
        perm_betas = result["perm_betas"]
        
        # Calculate bands
        y_obs_sorted, band_lo, band_hi = qq_bands_from_perms(
            obs_df["p_beta"].values, obs_df["beta"].values, perm_pvals, perm_betas, alpha=0.05
        )
        
        # Scatter plot
        base_color = color_map.get(label, None)
        ax, _top, theor_x = qq_signed_logp(
            obs_df["p_beta"].values,
            obs_df["beta"].values,
            ax,
            label=f"{label}",
            genes=obs_df["gene"].values,
            top_n=0,
            gene_set=None,
            gene_set_color=None,
            base_color=base_color
        )
        
        # Fill bands if shapes match
        K = min(len(theor_x), len(band_lo), len(band_hi), len(y_obs_sorted))
        if K > 0:
            ax.fill_between(theor_x[:K], band_lo[:K], band_hi[:K],
                            color=base_color if base_color is not None else None,
                            alpha=0.12, linewidth=0)
    
    # Diagonal guide & labels
    curx = ax.get_xlim()
    cury = ax.get_ylim()
    lim = max(abs(curx[0]), abs(curx[1]), abs(cury[0]), abs(cury[1]))
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=1, linestyle='--', alpha=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    
    ax.set_xlabel("Theoretical signed $-\\log_{10}(P)$")
    ax.set_ylabel("Observed signed $-\\log_{10}(P)$")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles[::-1], labels[::-1], title="Condition")
    plt.tight_layout()
    
    png_path = os.path.join(args.final_output_dir, "freedman_lane_permutation.png")
    plt.savefig(png_path, dpi=200)
    plt.close()
    print(f"Wrote {png_path}")
    
    print(f"\nIntegration complete! Results saved to {args.final_output_dir}")


if __name__ == "__main__":
    main()
