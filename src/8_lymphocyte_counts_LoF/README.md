# Regulator-burden correlation
Regulator-burden correlation analysis for lymphocyte count. Aims to identify core genes for different conditions (e.g., Stim8hr, Stime48hr) and their regulators.

---
## QQ-plots of regulator-burden correlations
`lymph_reg_burden_correlation.ipynb` does the following:
1. plot gene-level regulator-burden QQ-plots for a given trait
2. compares signed -log10P values across conditions (Stim8hr vs Stim48hr) and performs GO annotations for top genes
3. evaluates whether combining perturbation effects for different conditions gives us more signal
4. in a separate script, uses permutation methods outlined below to evaluate whether the distribution of signed -log10P values is calibrated under the null

### Freedman-Lane permutation to test calibration of burden-regulator correlation QQ plot
Uses Freedman-Lane residual permutation to test if the distribution of signed -log10P values is calibrated under the null. 

Goal: preserve the joint structure of $s_{het}$ with burden effect and with regulator effect.

1.	Fit the reduced model without $\text{Regulator effect}$:
$$\text{Burden effect} = s_{het}\alpha + r, \quad \hat r = \text{Burden effect} - s_{het}\hat\alpha.$$
2.	For each permutation $\pi$, permute the residuals $\hat r$ across samples: $\hat r^{(\pi)} = P_\pi \hat r$.
3.	Construct a permuted response that keeps $s_{het}$ fixed:
$$\text{Burden effect}^{(\pi)} = s_{het}\hat\alpha + \hat r^{(\pi)}.$$
4.	Fit the full model on ($\text{Burden effect}^{(\pi)}$, $\text{Regulator effect}$, $s_{het}$) and extract the test statistic for $\text{Regulator effect}$ for the QQ-plot.

This yields a valid null even when $\text{Regulator effect}$ and $s_{het}$ are correlated, and it preserves the correlation/dependence structure of $\text{Burden effect}$.

To run this in parallel across permutations:

```bash
chmod +x /oak/stanford/groups/pritch/users/rgoto/bin/GWT_causal_gene_discovery/src/5_reg_burden_correlation/freedman_lane_array.sh
/oak/stanford/groups/pritch/users/rgoto/bin/GWT_causal_gene_discovery/src/5_reg_burden_correlation/freedman_lane_array.sh
```

### An alternative way to test the calibration of burden-regulator correlation QQ plot
An alternative test if the distribution of signed -log10P values is calibrated under the null. This method randomly reverses the sign of the residual instead of permuting them across genes.

The goal remains the same: preserve the joint structure of $s_{het}$ with burden effect and with regulator effect.

1.	Fit the reduced model without $\text{Regulator effect}$:
$$\text{Burden effect} = s_{het}\alpha + r, \quad \hat r = \text{Burden effect} - s_{het}\hat\alpha.$$
2.	For each permutation $\pi$, permute the residuals $\hat r$ across samples: $\hat r^{(\pi)} = P_\pi \hat r$.
2. For each permutation￼$\pi$, randomly reverse the sign of residual $\hat r$: \hat r^{\pi}_i = s_i \hat r_i, where $s_i \sim \text{Bernoulli}(\tfrac{1}{2}$ independently takes values $\pm 1$. 
3.	Construct a permuted response that keeps $s_{het}$ fixed:
$$\text{Burden effect}^{(\pi)} = s_{het}\hat\alpha + \hat r^{(\pi)}.$$
4.	Fit the full model on ($\text{Burden effect}^{(\pi)}$, $\text{Regulator effect}$, $s_{het}$) and extract the test statistic for $\text{Regulator effect}$ for the QQ-plot.

This yields a valid null even when $\text{Regulator effect}$ and $s_{het}$ are correlated, and it preserves the correlation/dependence structure of $\text{Burden effect}$.

To run this in parallel:

```bash
chmod +x /oak/stanford/groups/pritch/users/rgoto/bin/GWT_causal_gene_discovery/src/5_reg_burden_correlation/sign_freedman_lane_array.sh
/oak/stanford/groups/pritch/users/rgoto/bin/GWT_causal_gene_discovery/src/5_reg_burden_correlation/sign_freedman_lane_array.sh
```

---
## Obtaining regulators of core gene(s) of interest
`core_gene_regulators.ipynb` obtains the regulators of a given (core) gene for Stim8hr and Stim48hr as well as their LoF burden test effects. Useful for identifying regulators unique to each condition-core gene pair. 