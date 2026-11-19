# Monetary Tightening or Generative AI? Isolating Post-2022 U.S. Labor Shocks via Conditional Forecasting

**Authors:** Eda Kal Arıca & Esra Şekerci  
**Affiliation:** Department of Statistics, Middle East Technical University (STAT 563 — Track A Empirical Project)  
**Project Website:** https://esrasekerci.github.io/generative_ai_labor_shock/  

---

## Abstract

The public release of generative artificial intelligence (AI) in late 2022 coincided directly with the sharpest Federal Reserve monetary tightening cycle in forty years, creating an acute macroeconomic identification challenge: did post-2022 labor cooling reflect policy-induced tightening or an exogenous technological shock? Standard dummy regressions cannot separate simultaneous disturbances that propagate through credit markets and real activity with distributed lags. 

We formulate a multivariate conditional forecasting framework over a tightly coupled six-variable monthly macro-financial system from FRED-MD (2010–2026), comprising policy rate (`FEDFUNDS`), broad money (`M2SL`), credit spread (`BAA10YM`), real output (`INDPRO`), and dual labor targets: job openings (`HWI`) and nonfarm payrolls (`PAYEMS`). Conditioning on realized post-2022 policy and financial paths via Waggoner–Zha (1999) exact structural projections, we test for systematic conditional forecast errors with Newey–West HAC standard errors. 

Our empirical results establish that while aggregate payroll employment tracks policy-conditioned expectations with statistical precision ($t = +0.085, p = 0.933$), job openings experienced an acute, highly significant structural collapse ($\bar{e} = -186.27$ thousand openings/month, $t = -5.729, p = 7.87 \times 10^{-7}$). This sharp bifurcation confirms that post-2022 monetary policy fully accounts for aggregate employment growth, but cannot explain the unprecedented contraction in job openings, isolating a distinct, idiosyncratic labor demand shock during the generative AI era.

Full data documentation, interactive figures, econometric specifications, and references are published on the [project website](https://esrasekerci.github.io/generative_ai_labor_shock/).

---

## Repository Structure

```text
.
├── fred_data.csv                                                # FRED-MD raw monthly dataset (St. Louis Fed)
├── requirements.txt                                             # Environment dependencies
├── proposal/project_proposal_edakalarica_esrasekerci.pdf
├── src/
│   ├── run_all.R                                                # Master reproduction script (executes full pipeline)
│   ├── data/
│   │   └── prepare_web_data.py                                  # Compiles fred_data.csv -> docs/assets/data/fred_macro.js
│   └── econometrics/
│       ├── 01_pretesting_and_cointegration.R                    # Unit-root battery, Johansen tests, outlier audit
│       └── 02_var_estimation_and_structural_irf.R             # VAR lag selection, Cholesky IRFs, counterfactuals, HAC tests
├── outputs/
│   ├── econometrics/                                            # Replicated empirical tables & test diagnostics (CSV)
│   │   ├── table_ur_levels_train.csv
│   │   ├── table_ur_trans_train.csv
│   │   ├── table_multivariate_outliers.csv
│   │   ├── table_johansen_cointegration.csv
│   │   ├── table_var_lag_selection.csv
│   │   ├── table_var_residual_diagnostics.csv
│   │   ├── table_residual_correlation_matrix.csv
│   │   └── table_hac_break_tests.csv
│   └── figures/                                                 # High-resolution (300 DPI) publication figures
│       ├── fig1_covid_outlier_audit.png (.pdf)
│       ├── fig2_structural_irfs.png (.pdf)
│       ├── fig3_ordering_robustness_irf.png (.pdf)
│       └── fig4_conditional_counterfactual_paths.png (.pdf)
└── docs/                                                        # GitHub Pages static site (main -> /docs)
    ├── .nojekyll
    ├── index.html
    └── assets/
        ├── css/style.css
        ├── js/app.js
        └── data/fred_macro.js
```

---

## Reproducing the Empirical Econometric Pipeline

All econometric procedures are implemented in R and adhere to course reference textbooks:
* **Enders, W. (2015)**: *Applied Econometric Time Series* (4th ed.), Ch. 4 (Unit Roots), Ch. 5 (VAR & IRF Identification), Ch. 6 (Cointegration).
* **Wei, W. W. S. (2019)**: *Multivariate Time Series Analysis and Applications*, Ch. 2 (Outlier Diagnostics & VAR Modeling).
* **Waggoner, D. F. & Zha, T. (1999)**: *Conditional Forecasts in Dynamic Multivariate Models*, Review of Economics and Statistics.
* **Newey, W. K. & West, K. D. (1987)**: *HAC Covariance Matrix Estimation*, Econometrica.

### Required R Packages

```r
install.packages(c("vars", "urca", "sandwich", "lmtest", "ggplot2", "patchwork", "scales"))
```

### Complete Replication (Single Command)

To run the full econometric estimation, pre-testing battery, outlier audit, structural IRF bootstrap simulations, conditional projections, and HAC break tests from raw data:

```bash
Rscript src/run_all.R
```

### Individual Replication Steps

1. **Unit-Root Battery, Cointegration & Outlier Diagnostics**:
   ```bash
   Rscript src/econometrics/01_pretesting_and_cointegration.R
   ```
2. **VAR Estimation, Structural IRFs, Counterfactuals & HAC Break Tests**:
   ```bash
   Rscript src/econometrics/02_var_estimation_and_structural_irf.R
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
