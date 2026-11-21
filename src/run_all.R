#!/usr/bin/env Rscript
# ==============================================================================
# run_all.R - Master Reproduction Script for Empirical Econometric Pipeline
# Project: Monetary Tightening or Generative AI? Isolating Post-2022 U.S. Labor Shocks
# Authors: Eda Kal Arıca & Esra Şekerci (METU Department of Statistics)
#
# Methodological references:
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.). Wiley.
#   - Wei, W. W. S. (2019). Multivariate Time Series Analysis and Applications. Wiley.
#   - Lenza, M. & Primiceri, G. E. (2022). Journal of Applied Econometrics.
#   - Giannone, D., Lenza, M. & Primiceri, G. E. (2015). Review of Economics and Statistics.
#   - Waggoner, D. F. & Zha, T. (1999). Review of Economics and Statistics.
#   - Newey, W. K. & West, K. D. (1987). Econometrica.
#   - Diebold, F. X. & Mariano, R. S. (1995). Journal of Business & Economic Statistics.
# ==============================================================================

cat("======================================================================\n")
cat("Starting Complete Empirical Econometric Pipeline Replication (5 Steps)\n")
cat("======================================================================\n\n")

t0 <- Sys.time()

cat(">>> Step 1/5: Unit-Root Battery, Wei Outlier Audit & Johansen Cointegration...\n")
source("src/econometrics/01_pretesting_and_cointegration.R")

cat("\n>>> Step 2/5: Common-Sample LR, COVID-VAR(2), VECM(r=2), IRFs & HAC Tests...\n")
source("src/econometrics/02_var_estimation_and_structural_irf.R")

cat("\n>>> Step 3/5: Exact Marginal-Likelihood Lenza-Primiceri BVAR & Fan Charts...\n")
source("src/econometrics/03_bvar_minnesota_covid_scaling.R")

cat("\n>>> Step 4/5: Falsification Battery (Pre-2022 Pseudo-Breaks & Sectoral Placebos)...\n")
source("src/econometrics/04_falsification_and_placebo_tests.R")

cat("\n>>> Step 5/5: Rolling-Origin Evaluation (h=1..12: MASE, CRPS & Diebold-Mariano)...\n")
source("src/econometrics/05_rolling_origin_evaluation.R")

t_elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat("\n======================================================================\n")
cat(sprintf("Complete 5-step pipeline finished cleanly in %s seconds.\n", t_elapsed))
cat("Generated Outputs:\n")
cat("  - outputs/econometrics/ (CSV tables, LaTeX tables 4-7, & diagnostics)\n")
cat("  - outputs/figures/      (Unified-punto IEEE Print & Web Figs 1-8)\n")
cat("======================================================================\n")
