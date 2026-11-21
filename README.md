# Monetary Tightening or Generative AI? Isolating Post-2022 U.S. Labor Shocks via Conditional Forecasting

**Authors:** Eda Kal Arıca & Esra Şekerci  
**Affiliation:** Department of Statistics, Middle East Technical University (STAT 563 — Track A Empirical Project)  
**Project Website:** https://esrasekerci.github.io/generative_ai_labor_shock/  

---

## Abstract

The public release of generative artificial intelligence (AI) in late 2022 coincided directly with the sharpest Federal Reserve monetary tightening cycle in forty years, creating an acute macroeconomic identification challenge: did post-2022 labor cooling reflect policy-induced tightening, omitted macroeconomic confounders, or an exogenous technological shock? Standard dummy regressions cannot separate simultaneous disturbances that propagate through credit markets and real activity with distributed lags.

We formulate a two-part multivariate forecasting study over a six-variable monthly macro-financial system from FRED-MD (2010–2026), comprising policy rate (`FEDFUNDS`), broad money (`M2SL`), credit spread (`BAA10YM`), industrial output (`INDPRO`), and dual labor targets: **JOLTS Total Nonfarm Job Openings** (`HWI`, directly reporting BLS `JTSJOL` since the 2017-09 FRED-MD vintage) and nonfarm payrolls (`PAYEMS`). To neutralize the 2020 COVID-19 outlier, we pair ten monthly **2020:03–2020:12** pulse dummies in frequentist models (`COVID-VAR(2)` and formal Case 3 `VECM(r = 2)` following Johansen rank $r = 2$) with exact marginal-likelihood **Lenza–Primiceri (2022)** volatility scaling in Bayesian VARs (Log Bayes Factor $= +414.1$).

* **Part A (Single-Origin Counterfactual & Falsification, 2022M11–2026M8):** Conditioning on realized post-2022 monetary and financial paths via exact Waggoner–Zha (1999) projections reveals a persistent vacancy-versus-payroll asymmetry: realized JOLTS openings (`HWI`) fell $-44.3\%$ below the cointegration-preserving `LP-BVAR(3)` counterfactual ($-20.5\%$ below `VECM(r = 2)`), whereas payrolls (`PAYEMS`) expanded by $+5.08$ million net jobs ($-3.8\%$ under `LP-BVAR(3)`, $+0.5\%, p = 0.509$ under `VECM(r = 2)`—zero break). Because conditional residuals also absorb fiscal, immigration, and reallocation confounders, we audit attribution via 13 pre-2022 quarterly pseudo-breaks (0/13 exceedances) and sectoral AI-exposure placebos (`USGOVT` $-0.2\%, p = 0.925$; high-AI `USFIRE` $-9.0\%$; rate-sensitive `USCONS` $-13.3\%$).
* **Part B (Rolling-Origin Multi-Horizon Evaluation, $N = 73$ origins, $h = 1\text{–}12$):** Across 38 clean pre-COVID and 35 post-2022 expanding origins (zero post-2022 tuning), **`VECM(r = 2)` ranks #1 on both targets** (`HWI` MASE $= 2.68$, CRPS $= 335.8$; `PAYEMS` MASE $= 1.55$, CRPS $= 214.7$, cutting CRPS by $-54.3\%$ vs. `COVID-VAR(2)`, Diebold–Mariano HLN $t = -2.54, p = 0.013$), followed by `LP-BVAR(3)` on payrolls (CRPS $= 375.0$), while `Unscaled-BVAR(2)` finishes last on payrolls (CRPS $= 539.5$).

Full data documentation, interactive figures, econometric specifications, and references are published on the [project website](https://esrasekerci.github.io/generative_ai_labor_shock/).

---

## Repository Structure

```text
.
├── fred_data.csv                                                # FRED-MD raw monthly dataset (St. Louis Fed)
├── requirements.txt                                             # Python environment dependencies
├── proposal/project_proposal_edakalarica_esrasekerci.pdf
├── src/
│   ├── run_all.R                                                # Master reproduction script (executes full 5-step pipeline)
│   ├── data/
│   │   └── prepare_web_data.py                                  # Compiles fred_data.csv -> docs/assets/data/fred_macro.js
│   └── econometrics/
│       ├── theme_publication.R                                  # Unified IEEE print (3.5 in) & web widescreen (8.4 in) theme
│       ├── 01_pretesting_and_cointegration.R                    # Unit-root battery, Wei outlier audit, Johansen (r = 2)
│       ├── 02_var_estimation_and_structural_irf.R               # Common-sample LR, COVID-VAR(2), VECM(r = 2), IRFs, WZ, HAC
│       ├── 03_bvar_minnesota_covid_scaling.R                    # Exact ML Lenza-Primiceri BVAR(2)/BVAR(3) & fan charts
│       ├── 04_falsification_and_placebo_tests.R                 # 13 pre-2022 pseudo-breaks & 6 sectoral AI-exposure placebos
│       └── 05_rolling_origin_evaluation.R                       # Expanding-window (N = 73, h = 1..12) MASE, CRPS & DM tests
├── outputs/
│   ├── econometrics/                                            # Replicated empirical tables & test diagnostics (CSV & LaTeX)
│   │   ├── table_ur_levels_train.csv
│   │   ├── table_ur_trans_train.csv
│   │   ├── table_ur_trans_full.csv
│   │   ├── table_multivariate_outliers.csv
│   │   ├── table_johansen_cointegration.csv
│   │   ├── table_vecm_cointegrating_vectors.csv
│   │   ├── table_var_lag_selection.csv
│   │   ├── table_var_residual_diagnostics.csv
│   │   ├── table_residual_correlation_matrix.csv
│   │   ├── table_counterfactual_trajectories.csv
│   │   ├── table_hac_break_tests.csv / table4_hac_break_tests.tex
│   │   ├── table_bvar_comparison.csv / table5_bvar_comparison.tex
│   │   ├── table_bvar_hyperparameters.csv
│   │   ├── table_falsification_placebo.csv / table_falsification_sectoral.csv / table6_falsification_tests.tex
│   │   └── table_rolling_origin_eval.csv / table_rolling_origin_by_horizon.csv / table7_rolling_origin_dm.tex
│   └── figures/                                                 # High-resolution (300 DPI) IEEE print & web figures
│       ├── fig1_covid_outlier_audit.png (.pdf)                  # & fig1_web_outlier_audit.png
│       ├── fig2_structural_irfs.png (.pdf)                      # & fig2_web_structural_irfs.png
│       ├── fig3_ordering_robustness_irf.png (.pdf)              # & fig3_web_ordering_robustness.png
│       ├── fig4_conditional_counterfactual_paths.png (.pdf)     # & fig4_web_counterfactuals.png
│       ├── fig5_bvar_posterior_bands.png (.pdf)                 # & fig5_web_bvar_bands.png
│       ├── fig6_bvar_vs_var_counterfactual.png (.pdf)           # & fig6_web_bvar_comparison.png
│       ├── fig7_falsification_placebo.png (.pdf)                # & fig7_web_falsification_placebo.png
│       └── fig8_rolling_origin_horizons.png (.pdf)              # & fig8_web_rolling_origin_horizons.png
└── docs/                                                        # GitHub Pages static site (main -> /docs)
    ├── .nojekyll
    ├── index.html
    └── assets/
        ├── css/style.css
        ├── js/app.js
        ├── data/fred_macro.js
        └── images/                                              # Synchronized 300 DPI widescreen web figures
```

---

## Reproducing the Empirical Econometric Pipeline

All econometric procedures are fully deterministic (`set.seed(563)`), implemented in R, and adhere strictly to course reference textbooks and literature:
* **Enders, W. (2015)**: *Applied Econometric Time Series* (4th ed.), Ch. 4 (Unit Roots), Ch. 5 (Common-Sample LR Lag Selection, VAR & IRF Ordering Audit), Ch. 6 (Johansen Cointegration & VECM).
* **Wei, W. W. S. (2019)**: *Multivariate Time Series Analysis and Applications*, Ch. 2 (Secs. 2.5–2.7: VAR Diagnostics, MIO/MAO Outlier Detection & Pulse Intervention).
* **Cryer, J. D. & Chan, K.-S. (2008)** & **Shumway, R. H. & Stoffer, D. S. (2011)**: Time-series stationarity and residual diagnostics.
* **Waggoner, D. F. & Zha, T. (1999)**: *Conditional Forecasts in Dynamic Multivariate Models*, Review of Economics and Statistics (exact minimum-MSE projection & conditional predictive density draws).
* **Giannone, D., Lenza, M. & Primiceri, G. E. (2015)** & **Lenza, M. & Primiceri, G. E. (2022)**: Conjugate Minnesota prior marginal-likelihood selection & COVID-19 innovation volatility scaling.
* **Diebold, F. X. & Mariano, R. S. (1995)**, **Harvey, D., Leybourne, S. & Newbold, P. (1997)**, **Hyndman, R. J. & Koehler, A. B. (2006)**, **Gneiting, T. & Raftery, A. E. (2007)**, and **Newey, W. K. & West, K. D. (1987)**: Multi-horizon forecast accuracy (MASE, CRPS, DM-HLN tests) and HAC inference.

### Required R Packages

```r
install.packages(c("vars", "urca", "MASS", "sandwich", "lmtest", "ggplot2", "patchwork", "scales"))
```

### Complete Replication (Single Command)

To reproduce **every** table, diagnostic CSV, LaTeX table, IEEE print figure, and GitHub Pages figure from raw `fred_data.csv` (~25 seconds):

```bash
Rscript src/run_all.R
```

### Individual Replication Steps

1. **Unit-Root Battery, Multivariate Outlier Audit & Johansen Cointegration**:
   ```bash
   Rscript src/econometrics/01_pretesting_and_cointegration.R
   ```
2. **Common-Sample Lag Selection, COVID-Controlled VAR(2), Formal VECM(r = 2), IRFs, Counterfactuals & HAC Tests**:
   ```bash
   Rscript src/econometrics/02_var_estimation_and_structural_irf.R
   ```
3. **Exact Marginal-Likelihood Minnesota BVAR with Lenza–Primiceri COVID Scaling & Predictive Fan Charts**:
   ```bash
   Rscript src/econometrics/03_bvar_minnesota_covid_scaling.R
   ```
4. **Attribution Falsification Battery (13 Pre-2022 Pseudo-Breaks & 6 Sectoral AI-Exposure Placebos)**:
   ```bash
   Rscript src/econometrics/04_falsification_and_placebo_tests.R
   ```
5. **Rolling-Origin Multi-Horizon Forecast Evaluation ($N = 73$ Origins, $h = 1\text{–}12$, MASE, CRPS & Diebold–Mariano)**:
   ```bash
   Rscript src/econometrics/05_rolling_origin_evaluation.R
   ```

---

## Web Site Data Compilation & Local Preview

To rebuild the interactive dataset for GitHub Pages:

```bash
python3 src/data/prepare_web_data.py
```

To preview the academic website locally:

```bash
python3 -m http.server 8000 --directory docs
```
