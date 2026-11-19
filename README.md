# Monetary Tightening or Generative AI? Isolating Post-2022 U.S. Labor Shocks via Conditional Forecasting

**Authors:** Eda Kal Arıca & Esra Şekerci  
**Affiliation:** Department of Statistics, Middle East Technical University (STAT 563 — Track A Empirical Project)  
**Project Website:** https://esrasekerci.github.io/generative_ai_labor_shock/  

---

## Abstract

The public release of generative artificial intelligence (AI) in late 2022 coincided directly with the sharpest Federal Reserve monetary tightening cycle in forty years, creating an acute macroeconomic identification challenge: did post-2022 labor cooling reflect policy-induced tightening or an exogenous technological shock? Standard dummy regressions cannot separate simultaneous disturbances that propagate through credit markets and real activity with distributed lags.

We formulate a multivariate conditional forecasting framework over a tightly coupled six-variable monthly macro-financial system from FRED-MD (2010–2026), comprising policy rate (`FEDFUNDS`), broad money (`M2SL`), credit spread (`BAA10YM`), real output (`INDPRO`), and dual labor targets: job openings (`HWI`) and nonfarm payrolls (`PAYEMS`). Grounded in textbook econometric methodology, we execute a confirmatory unit-root battery, Johansen cointegration rank tests ($r = 2$), Wei (2019) multivariate outlier audits, COVID-controlled VAR(2) estimation, and exact marginal-likelihood Minnesota Bayesian VARs with Lenza–Primiceri (2022) COVID-19 volatility scaling (Log Bayes Factor $= +414.1$).

Conditioning on realized post-2022 monetary and financial paths via exact Waggoner–Zha (1999) minimum-MSE projections and 500-draw predictive fan charts, we uncover a stark **12-to-1 proportional bifurcation**: under the cointegration-preserving Levels LP-BVAR(3), realized job openings (`HWI`) plunged **$-44.3\%$** below policy-conditioned expectations ($\bar{e} = -122.8$ thousand/month, $P(\bar{e} < 0) = 1.000$), breaking completely below the 95% predictive fan chart, whereas aggregate payrolls (`PAYEMS`) expanded by $+5.08$ million net jobs and tracked inside the 2022–2024 predictive bands ($-3.8\%$ end gap). This sharp asymmetry isolates an idiosyncratic hiring-freeze shock orthogonal to monetary policy during the generative AI era.

Full data documentation, interactive figures, econometric specifications, and references are published on the [project website](https://esrasekerci.github.io/generative_ai_labor_shock/).

---

## Repository Structure (Completed Through Phase 2)

```text
.
├── fred_data.csv                                                # FRED-MD raw monthly dataset (St. Louis Fed)
├── requirements.txt                                             # Python environment dependencies
├── proposal/project_proposal_edakalarica_esrasekerci.pdf
├── src/
│   ├── run_all.R                                                # Master reproduction script (executes full 3-step pipeline)
│   ├── data/
│   │   └── prepare_web_data.py                                  # Compiles fred_data.csv -> docs/assets/data/fred_macro.js
│   └── econometrics/
│       ├── theme_publication.R                                  # Unified IEEE print (3.5 in) & web widescreen (8.4 in) theme
│       ├── 01_pretesting_and_cointegration.R                    # Phase 1: Unit-root battery, Wei outlier audit, Johansen tests
│       ├── 02_var_estimation_and_structural_irf.R               # Phase 1: Common-sample LR, COVID-VAR(2), IRFs, WZ, HAC tests
│       └── 03_bvar_minnesota_covid_scaling.R                    # Phase 2: Exact ML Lenza-Primiceri BVAR(2)/BVAR(3) & fan charts
├── outputs/
│   ├── econometrics/                                            # Replicated empirical tables & test diagnostics (CSV & LaTeX)
│   │   ├── table_ur_levels_train.csv
│   │   ├── table_ur_trans_train.csv
│   │   ├── table_ur_trans_full.csv
│   │   ├── table_multivariate_outliers.csv
│   │   ├── table_johansen_cointegration.csv
│   │   ├── table_var_lag_selection.csv
│   │   ├── table_var_residual_diagnostics.csv
│   │   ├── table_residual_correlation_matrix.csv
│   │   ├── table_counterfactual_trajectories.csv
│   │   ├── table_hac_break_tests.csv
│   │   ├── table4_hac_break_tests.tex
│   │   ├── table_bvar_comparison.csv
│   │   ├── table_bvar_hyperparameters.csv
│   │   └── table5_bvar_comparison.tex
│   └── figures/                                                 # High-resolution (300 DPI) IEEE print & web figures
│       ├── fig1_covid_outlier_audit.png (.pdf)                  # & fig1_web_outlier_audit.png
│       ├── fig2_structural_irfs.png (.pdf)                      # & fig2_web_structural_irfs.png
│       ├── fig3_ordering_robustness_irf.png (.pdf)              # & fig3_web_ordering_robustness.png
│       ├── fig4_conditional_counterfactual_paths.png (.pdf)     # & fig4_web_counterfactuals.png
│       ├── fig5_bvar_posterior_bands.png (.pdf)                 # & fig5_web_bvar_bands.png
│       └── fig6_bvar_vs_var_counterfactual.png (.pdf)           # & fig6_web_bvar_comparison.png
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
* **Enders, W. (2015)**: *Applied Econometric Time Series* (4th ed.), Ch. 4 (Unit Roots), Ch. 5 (Common-Sample LR Lag Selection, VAR & IRF Ordering Audit), Ch. 6 (Johansen Cointegration).
* **Wei, W. W. S. (2019)**: *Multivariate Time Series Analysis and Applications*, Ch. 2 (Secs. 2.5–2.7: VAR Diagnostics, MIO/MAO Outlier Detection & Pulse Intervention).
* **Cryer, J. D. & Chan, K.-S. (2008)** & **Shumway, R. H. & Stoffer, D. S. (2011)**: Time-series stationarity and residual diagnostics.
* **Waggoner, D. F. & Zha, T. (1999)**: *Conditional Forecasts in Dynamic Multivariate Models*, Review of Economics and Statistics (exact minimum-MSE projection & conditional predictive density draws).
* **Giannone, D., Lenza, M. & Primiceri, G. E. (2015)** & **Lenza, M. & Primiceri, G. E. (2022)**: Conjugate Minnesota prior marginal-likelihood selection & COVID-19 innovation volatility scaling.
* **Newey, W. K. & West, K. D. (1987)**: *HAC Covariance Matrix Estimation*, Econometrica.

### Required R Packages

```r
install.packages(c("vars", "urca", "MASS", "sandwich", "lmtest", "ggplot2", "patchwork", "scales"))
```

### Complete Replication (Single Command)

To reproduce **every** table, diagnostic CSV, LaTeX table, IEEE print figure, and GitHub Pages figure across Phases 1 and 2 from raw `fred_data.csv` (~21 seconds):

```bash
Rscript src/run_all.R
```

### Individual Replication Steps

1. **Unit-Root Battery, Multivariate Outlier Audit & Johansen Cointegration**:
   ```bash
   Rscript src/econometrics/01_pretesting_and_cointegration.R
   ```
2. **Common-Sample Lag Selection, COVID-Controlled VAR(2), Cholesky IRFs, Counterfactuals & HAC Tests**:
   ```bash
   Rscript src/econometrics/02_var_estimation_and_structural_irf.R
   ```
3. **Exact Marginal-Likelihood Minnesota BVAR with Lenza–Primiceri COVID Scaling & Predictive Fan Charts**:
   ```bash
   Rscript src/econometrics/03_bvar_minnesota_covid_scaling.R
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
