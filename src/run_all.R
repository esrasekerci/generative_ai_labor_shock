#!/usr/bin/env Rscript
# ==============================================================================
# run_all.R - Master Reproduction Script for Empirical Econometric Pipeline
# Project: Monetary Tightening or Generative AI? Isolating Post-2022 U.S. Labor Shocks
# Authors: Eda Kal Arıca & Esra Şekerci (METU Department of Statistics)
#
# Methodological references:
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.). Wiley.
#   - Wei, W. W. S. (2019). Multivariate Time Series Analysis and Applications. Wiley.
#   - Waggoner, D. F. & Zha, T. (1999). Review of Economics and Statistics.
#   - Newey, W. K. & West, K. D. (1987). Econometrica.
# ==============================================================================

cat("======================================================================\n")
cat("Starting Complete Empirical Econometric Pipeline Replication\n")
cat("======================================================================\n\n")

t0 <- Sys.time()

cat(">>> Step 1/2: Executing Unit-Root Battery, Cointegration & Outlier Diagnostics...\n")
source("src/econometrics/01_pretesting_and_cointegration.R")

cat("\n>>> Step 2/2: Executing VAR Lag Selection, Structural IRFs, & Conditional Forecasts...\n")
source("src/econometrics/02_var_estimation_and_structural_irf.R")

t_elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat("\n======================================================================\n")
cat(sprintf("Pipeline completed successfully in %s seconds.\n", t_elapsed))
cat("Generated Outputs:\n")
cat("  - outputs/econometrics/ (CSV tables & diagnostics)\n")
cat("  - outputs/figures/      (Publication figures: Fig 1-4)\n")
cat("======================================================================\n")
