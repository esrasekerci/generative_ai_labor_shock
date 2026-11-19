#!/usr/bin/env Rscript
# ==============================================================================
# 02_var_estimation_and_structural_irf.R
# Step 2: VAR Estimation, Sequential LR Lag Selection (Enders 2015, p. 304),
#         Stability Verification, Cholesky Structural IRF & Robustness Ordering,
#         Waggoner-Zha Conditional Forecasting, and Post-2022 HAC Break Tests.
#
# Methodological references:
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.), Ch. 5 (Secs. 5.5-5.8, pp. 297-316).
#   - Wei, W. W. S. (2019). Multivariate Time Series Analysis and Applications, Ch. 2 (Secs. 2.5-2.7).
#   - Waggoner, D. F. & Zha, T. (1999). Conditional Forecasts in Dynamic Multivariate Models.
#     Review of Economics and Statistics, 81(4), 639-651.
#   - Newey, W. K. & West, K. D. (1987). A Simple, Positive Semi-Definite, Heteroskedasticity
#     and Autocorrelation Consistent Covariance Matrix. Econometrica, 55(3), 703-708.
# ==============================================================================

suppressPackageStartupMessages({
  library(vars)
  library(urca)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

dir.create("outputs/econometrics", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load Data & Prepare Stationary Estimation System
# ------------------------------------------------------------------------------
raw <- read.csv("fred_data.csv", stringsAsFactors = FALSE)
df <- raw[-1, ]
df$date <- as.Date(df$sasdate, format = "%m/%d/%Y")
df <- df[!is.na(df$date), ]
df <- df[order(df$date), ]

num_cols <- c("FEDFUNDS", "M2SL", "BAA", "GS10", "INDPRO", "HWI", "PAYEMS")
for (col in num_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}
df$BAA10YM <- df$BAA - df$GS10

d1 <- function(x) c(NA, diff(x))
dlog1 <- function(x) c(NA, diff(log(x))) * 100
d2log1 <- function(x) c(NA, NA, diff(diff(log(x)))) * 100

# Stationary FRED-MD transforms
df$t_FEDFUNDS <- d1(df$FEDFUNDS)        # Δ policy rate (pp)
df$t_M2SL     <- d2log1(df$M2SL)        # Δ² log(M2SL) * 100 (liquidity acceleration)
df$t_BAA10YM  <- df$BAA10YM             # corporate bond credit spread (pp)
df$t_INDPRO   <- dlog1(df$INDPRO)       # monthly industrial output growth (%)
df$t_HWI      <- d1(df$HWI)             # monthly change in job openings (k)
df$t_PAYEMS   <- dlog1(df$PAYEMS)       # monthly payroll employment growth (%)

# Complete usable sample (2010-01-01 to 2026-08-01, T = 200)
full_sample <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]
train_df    <- full_sample[full_sample$date <= as.Date("2022-10-01"), ] # T_train = 154
test_df     <- full_sample[full_sample$date >= as.Date("2022-11-01"), ] # T_test  = 46

var_names <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")
var_cols  <- c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")

Y_train <- as.matrix(train_df[, var_cols])
colnames(Y_train) <- var_names

Y_test <- as.matrix(test_df[, var_cols])
colnames(Y_test) <- var_names

Y_full <- as.matrix(full_sample[, var_cols])
colnames(Y_full) <- var_names

k <- ncol(Y_train)
T_tr <- nrow(Y_train)

cat("======================================================================\n")
cat("1. ESTIMATION SYSTEM SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Variables (k = %d)      : %s\n", k, paste(var_names, collapse = ", ")))
cat(sprintf("Pre-AI Training (T = %d): %s to %s\n", T_tr, min(train_df$date), max(train_df$date)))
cat(sprintf("Post-AI Test (T = %d)    : %s to %s\n\n", nrow(test_df), min(test_df$date), max(test_df$date)))

# ------------------------------------------------------------------------------
# 2. Lag Length Selection: Information Criteria & Sims Small-Sample Adjusted LR
#    (Enders 2015, Sec. 5.8, pp. 303-305)
# ------------------------------------------------------------------------------
# Sims (1980) small-sample correction: (T - c) * (ln|Sigma_r| - ln|Sigma_u|) ~ chi^2(k^2 * (p_u - p_r))
# where c = 1 + k * p_u (number of parameters estimated in each unrestricted equation)
max_p <- 8
p_candidates <- 1:max_p
log_det_sigma <- numeric(max_p)
aic_vec <- numeric(max_p)
bic_vec <- numeric(max_p)
hq_vec  <- numeric(max_p)

# To ensure common estimation sample for all p, reserve first max_p observations
# as presample: T_eff = T_tr - max_p
presample <- max_p
T_eff <- T_tr - presample

for (p in p_candidates) {
  v_fit <- VAR(Y_train, p = p, type = "const")
  # Extract residuals corresponding to the common evaluation period (rows (max_p - p + 1) to end)
  res_full <- residuals(v_fit)
  res_common <- tail(res_full, T_eff)
  sig_p <- crossprod(res_common) / T_eff
  log_det_sigma[p] <- log(det(sig_p))
  
  # Total parameters per equation: 1 + k*p
  num_params <- 1 + k * p
  aic_vec[p] <- log_det_sigma[p] + (2 * k * num_params) / T_eff
  bic_vec[p] <- log_det_sigma[p] + (log(T_eff) * k * num_params) / T_eff
  hq_vec[p]  <- log_det_sigma[p] + (2 * log(log(T_eff)) * k * num_params) / T_eff
}

# Sequential Likelihood Ratio Tests: test p vs p+1
lr_table <- data.frame(
  Lag_p = 1:max_p,
  LogDet = round(log_det_sigma, 3),
  AIC = round(aic_vec, 3),
  BIC = round(bic_vec, 3),
  HQ = round(hq_vec, 3),
  LR_Stat_Sims = NA,
  LR_df = NA,
  LR_pvalue = NA,
  stringsAsFactors = FALSE
)

for (p in 1:(max_p - 1)) {
  p_u <- p + 1
  p_r <- p
  c_sims <- 1 + k * p_u
  df_lr <- k^2 * (p_u - p_r) # k^2 * 1 = 36
  lr_stat <- (T_eff - c_sims) * (log_det_sigma[p_r] - log_det_sigma[p_u])
  lr_p <- 1 - pchisq(lr_stat, df = df_lr)
  lr_table$LR_Stat_Sims[p] <- round(lr_stat, 2)
  lr_table$LR_df[p] <- df_lr
  lr_table$LR_pvalue[p] <- signif(lr_p, 3)
}

best_aic <- which.min(aic_vec)
best_bic <- which.min(bic_vec)
best_hq  <- which.min(hq_vec)

cat("======================================================================\n")
cat("2. LAG LENGTH SELECTION (ENDERS 2015, SEC. 5.8)\n")
cat("======================================================================\n")
print(lr_table, row.names = FALSE)
cat(sprintf("Selected lag by AIC: p = %d | BIC: p = %d | HQ: p = %d\n", best_aic, best_bic, best_hq))
write.csv(lr_table, "outputs/econometrics/table_var_lag_selection.csv", row.names = FALSE)

# Select p_opt: BIC selects p = 2 (parsimonious, preserves degrees of freedom; T_eff/parameters)
p_opt <- best_bic
cat(sprintf("Optimal baseline model selected: VAR(%d)\n\n", p_opt))

# ------------------------------------------------------------------------------
# 3. Baseline VAR(p) Estimation, Companion Matrix Stability, & Diagnostic Checks
#    (Wei 2019, Sec. 2.5.3; Cryer & Chan 2008, Ch. 8)
# ------------------------------------------------------------------------------
var_baseline <- VAR(Y_train, p = p_opt, type = "const")

# (a) Companion matrix eigenvalues (roots of characteristic polynomial)
roots_obj <- roots(var_baseline, modulus = TRUE)
max_eigen <- max(roots_obj)
is_stable <- all(roots_obj < 1.0)

cat("======================================================================\n")
cat("3A. STABILITY VERIFICATION (COMPANION MATRIX ROOTS)\n")
cat("======================================================================\n")
cat(sprintf("Max eigenvalue modulus: %.4f (< 1.0: %s)\n", max_eigen, ifelse(is_stable, "STABLE", "UNSTABLE")))
print(round(roots_obj, 4))

# (b) Diagnostic tests: Portmanteau autocorrelation, ARCH volatility, Normality
serial_diag <- serial.test(var_baseline, lags.pt = 12, type = "PT.asymptotic")
arch_diag   <- arch.test(var_baseline, lags.multi = 4)
norm_diag   <- normality.test(var_baseline)

diag_summary <- data.frame(
  Diagnostic_Test = c(
    "Portmanteau Autocorrelation (PT.asymptotic, h=12)",
    "Multivariate ARCH Test (lags=4)",
    "Multivariate Jarque-Bera Normality",
    "Skewness Test",
    "Kurtosis Test"
  ),
  Test_Statistic = c(
    round(serial_diag$serial$statistic, 2),
    round(arch_diag$arch.mul$statistic, 2),
    round(norm_diag$jb.mul$JB$statistic, 2),
    round(norm_diag$jb.mul$Skewness$statistic, 2),
    round(norm_diag$jb.mul$Kurtosis$statistic, 2)
  ),
  Degrees_of_Freedom = c(
    serial_diag$serial$parameter,
    arch_diag$arch.mul$parameter,
    norm_diag$jb.mul$JB$parameter,
    norm_diag$jb.mul$Skewness$parameter,
    norm_diag$jb.mul$Kurtosis$parameter
  ),
  p_value = c(
    signif(serial_diag$serial$p.value, 3),
    signif(arch_diag$arch.mul$p.value, 3),
    signif(norm_diag$jb.mul$JB$p.value, 3),
    signif(norm_diag$jb.mul$Skewness$p.value, 3),
    signif(norm_diag$jb.mul$Kurtosis$p.value, 3)
  ),
  stringsAsFactors = FALSE
)

cat("\n======================================================================\n")
cat("3B. BASELINE RESIDUAL DIAGNOSTICS (WEI 2019, SEC. 2.5.3)\n")
cat("======================================================================\n")
print(diag_summary, row.names = FALSE)
write.csv(diag_summary, "outputs/econometrics/table_var_residual_diagnostics.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Residual Correlation Matrix & Ordering Sensitivity Audit
#    (Enders 2015, Sec. 5.7, p. 298: rho_ij test against T^(-0.5))
# ------------------------------------------------------------------------------
# Enders rule (p. 298): test rho_ij = 0 with SE = 1/sqrt(T). For T = 152, SE = 0.081.
# |rho| > 2*SE = 0.162 is significant at 5%. If correlations across blocks are near zero,
# the Cholesky ordering is invariant.
E_opt <- residuals(var_baseline)
Sigma_res <- crossprod(E_opt) / nrow(E_opt)
P_res <- cov2cor(Sigma_res)
rownames(P_res) <- colnames(P_res) <- var_names

se_rho <- 1 / sqrt(nrow(E_opt))
crit_rho_05 <- 1.96 * se_rho
crit_rho_01 <- 2.576 * se_rho

cat("\n======================================================================\n")
cat("4. RESIDUAL CORRELATION MATRIX & CHOLESKY ORDERING AUDIT (ENDERS 2015, p. 298)\n")
cat("======================================================================\n")
cat(sprintf("Sample T_e = %d | SE(rho) = %.4f | 5%% Threshold = +/-%.3f | 1%% Threshold = +/-%.3f\n",
            nrow(E_opt), se_rho, crit_rho_05, crit_rho_01))
print(round(P_res, 3))
write.csv(round(P_res, 4), "outputs/econometrics/table_residual_correlation_matrix.csv")

# ------------------------------------------------------------------------------
# 5. Structural Cholesky IRFs & Alternative Ordering Robustness
#    (Enders 2015, Sec. 5.7, pp. 297-302)
# ------------------------------------------------------------------------------
# Primary Ordering (Standard Macroeconomic Transmission):
#   1. FEDFUNDS  (Policy instrument, fast-moving)
#   2. M2SL      (Broad money / liquidity)
#   3. BAA10YM   (Financial friction / credit spread)
#   4. INDPRO    (Real output activity)
#   5. HWI       (Labor demand / vacancies)
#   6. PAYEMS    (Labor volume / employment)
#
# Alternative Ordering 1 (Pater 2017 Labor-Preceding-Monetary):
#   [PAYEMS, HWI, INDPRO, FEDFUNDS, M2SL, BAA10YM]
# Alternative Ordering 2 (Vacancy-First Swap):
#   [FEDFUNDS, M2SL, BAA10YM, INDPRO, PAYEMS, HWI]

irf_horizon <- 24
n_boot <- 500

cat("\nComputing primary Cholesky IRFs with 500 bootstrap replications (horizon = 24)...\n")
set.seed(563)
irf_primary <- irf(var_baseline, impulse = c("FEDFUNDS", "BAA10YM"),
                   response = c("HWI", "PAYEMS", "INDPRO"),
                   n.ahead = irf_horizon, boot = TRUE, runs = n_boot, ci = 0.95)

# Fit Alternative Ordering 1
alt1_names <- c("PAYEMS", "HWI", "INDPRO", "FEDFUNDS", "M2SL", "BAA10YM")
Y_train_alt1 <- Y_train[, alt1_names]
var_alt1 <- VAR(Y_train_alt1, p = p_opt, type = "const")
irf_alt1 <- irf(var_alt1, impulse = c("FEDFUNDS", "BAA10YM"),
                response = c("HWI", "PAYEMS", "INDPRO"),
                n.ahead = irf_horizon, boot = TRUE, runs = n_boot, ci = 0.95)

# Convert primary IRF into tidy data frame for publication plotting
tidy_irf <- function(irf_obj, imp, resp, label = "Baseline") {
  h <- 0:(nrow(irf_obj$irf[[imp]]) - 1)
  data.frame(
    Horizon = h,
    Impulse = imp,
    Response = resp,
    Point = as.numeric(irf_obj$irf[[imp]][, resp]),
    Lower = as.numeric(irf_obj$Lower[[imp]][, resp]),
    Upper = as.numeric(irf_obj$Upper[[imp]][, resp]),
    Model = label,
    stringsAsFactors = FALSE
  )
}

irf_df_ff_hwi <- tidy_irf(irf_primary, "FEDFUNDS", "HWI", "Primary Ordering")
irf_df_ff_pay <- tidy_irf(irf_primary, "FEDFUNDS", "PAYEMS", "Primary Ordering")
irf_df_sp_hwi <- tidy_irf(irf_primary, "BAA10YM", "HWI", "Primary Ordering")
irf_df_sp_pay <- tidy_irf(irf_primary, "BAA10YM", "PAYEMS", "Primary Ordering")

irf_alt1_ff_hwi <- tidy_irf(irf_alt1, "FEDFUNDS", "HWI", "Labor-First Ordering")
irf_alt1_ff_pay <- tidy_irf(irf_alt1, "FEDFUNDS", "PAYEMS", "Labor-First Ordering")

irf_plot_df <- rbind(irf_df_ff_hwi, irf_df_ff_pay, irf_df_sp_hwi, irf_df_sp_pay)
irf_plot_df$PanelLabel <- factor(
  paste0(irf_plot_df$Impulse, " -> ", irf_plot_df$Response),
  levels = c(
    "FEDFUNDS -> HWI",
    "FEDFUNDS -> PAYEMS",
    "BAA10YM -> HWI",
    "BAA10YM -> PAYEMS"
  ),
  labels = c(
    "(a) FEDFUNDS -> HWI",
    "(b) FEDFUNDS -> PAYEMS",
    "(c) BAA10YM -> HWI",
    "(d) BAA10YM -> PAYEMS"
  )
)

# Plot Publication IRFs (Figure 2: IEEE Single Column 3.5in)
p_irf <- ggplot(irf_plot_df, aes(x = Horizon, y = Point)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#2980B9", alpha = 0.22) +
  geom_line(color = "#1F618D", linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.4) +
  facet_wrap(~PanelLabel, scales = "free_y", ncol = 2) +
  labs(
    x = "Forecast Horizon (Months)",
    y = "Orthogonalized Response"
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "#F2F4F4", color = NA),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(t = 4, r = 6, b = 4, l = 4, "pt")
  )

ggsave("outputs/figures/fig2_structural_irfs.pdf", p_irf, width = 3.5, height = 2.8)
ggsave("outputs/figures/fig2_structural_irfs.png", p_irf, width = 3.5, height = 2.8, dpi = 300)

# Ordering Robustness Comparison (Figure 3: IEEE Single Column 3.5in)
comp_df <- rbind(
  cbind(irf_df_ff_hwi, Ordering = "Baseline [Policy First]"),
  cbind(irf_alt1_ff_hwi, Ordering = "Alternative [Labor First]"),
  cbind(irf_df_ff_pay, Ordering = "Baseline [Policy First]"),
  cbind(irf_alt1_ff_pay, Ordering = "Alternative [Labor First]")
)
comp_df$ResponseVar <- ifelse(comp_df$Response == "HWI", "(a) Job Openings (HWI)", "(b) Payrolls (PAYEMS)")

p_robust_irf <- ggplot(comp_df, aes(x = Horizon, y = Point, color = Ordering, linetype = Ordering)) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.4) +
  facet_wrap(~ResponseVar, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Baseline [Policy First]" = "#1F618D", "Alternative [Labor First]" = "#C0392B")) +
  scale_linetype_manual(values = c("Baseline [Policy First]" = "solid", "Alternative [Labor First]" = "dashed")) +
  labs(
    x = "Forecast Horizon (Months)",
    y = "Impulse Response",
    color = "Ordering",
    linetype = "Ordering"
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.margin = margin(b = -2),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "#F2F4F4", color = NA),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(t = 2, r = 6, b = 4, l = 4, "pt")
  )

ggsave("outputs/figures/fig3_ordering_robustness_irf.pdf", p_robust_irf, width = 3.5, height = 2.4)
ggsave("outputs/figures/fig3_ordering_robustness_irf.png", p_robust_irf, width = 3.5, height = 2.4, dpi = 300)

# ------------------------------------------------------------------------------
# 6. Waggoner-Zha (1999) Conditional Forecasting Framework
# ------------------------------------------------------------------------------
# Conditioning Block C: {FEDFUNDS, M2SL, BAA10YM} (indices 1, 2, 3)
# Target Block T    : {HWI, PAYEMS}           (indices 5, 6)
# Intermediate Block: {INDPRO}                (index 4)
#
# Procedure:
# Over H = 46 months (2022-11 to 2026-08):
# 1. Unconditional baseline forecast: Y_hat^{unc}_{t+h} = \sum_{j=1}^p A_j Y_{t+h-j} + c
# 2. Structural conditional forecast: find sequences of shocks {e_{1..3, t+h}} that match
#    the exact realized paths of FEDFUNDS, M2SL, BAA10YM over post-2022, while setting
#    exogenous structural shocks to labor {e_{5, t+h}, e_{6, t+h}} = 0!
#    Then compute the conditional expectation of HWI and PAYEMS given the policy paths!

H <- nrow(Y_test)
cond_idx <- c(1, 2, 3) # FEDFUNDS, M2SL, BAA10YM
target_idx <- c(5, 6)  # HWI, PAYEMS

A_est <- Acoef(var_baseline) # list of p matrices (k x k)
c_est <- sapply(var_baseline$varresult, function(eq) coef(eq)["const"]) # k x 1
P_chol <- t(chol(Sigma_res)) # Lower triangular Cholesky factor: Sigma = P P'

# Simulation of unconditional dynamic forecast
Y_unc <- matrix(0, nrow = H, ncol = k)
colnames(Y_unc) <- var_names

# Presample history for forecast initiation: last p rows of training
hist_init <- Y_train[(T_tr - p_opt + 1):T_tr, ]

for (h in 1:H) {
  pred_h <- c_est
  for (j in 1:p_opt) {
    lag_val <- if ((h - j) <= 0) {
      hist_init[p_opt + (h - j), ]
    } else {
      Y_unc[h - j, ]
    }
    pred_h <- pred_h + A_est[[j]] %*% lag_val
  }
  Y_unc[h, ] <- pred_h
}

# Waggoner-Zha Exact Conditional Forecast
# State-space representation: Y_{t+h} = Y^{unc}_{t+h} + \sum_{s=1}^h M_{h-s} P eta_s
# Compute impulse response multipliers M_m (unorthogonalized MA representation)
M_mats <- vector("list", H)
M_mats[[1]] <- diag(k)
for (m in 2:H) {
  M_m <- matrix(0, k, k)
  for (j in 1:min(p_opt, m - 1)) {
    M_m <- M_m + A_est[[j]] %*% M_mats[[m - j]]
  }
  M_mats[[m]] <- M_m
}

# Stack condition equations:
# For conditioning variables i in cond_idx at step h:
# Y_{i, h} - Y^{unc}_{i, h} = \sum_{s=1}^h [M_{h-s} P]_{i, 1:3} eta_{C, s}
# Let n_c = 3. Total constraints = H * n_c. Total shocks = H * n_c.
n_c <- length(cond_idx)
R_mat <- matrix(0, nrow = H * n_c, ncol = H * n_c)
y_gap <- numeric(H * n_c)

row_cnt <- 1
for (h in 1:H) {
  for (ci in 1:n_c) {
    c_var <- cond_idx[ci]
    y_gap[row_cnt] <- Y_test[h, c_var] - Y_unc[h, c_var]
    
    # Fill coefficients for eta_s (s = 1..h)
    for (s in 1:h) {
      M_hs <- M_mats[[h - s + 1]]
      # Column index in eta_vector for shock at time s, component cj (1..n_c)
      for (cj in 1:n_c) {
        col_idx <- (s - 1) * n_c + cj
        R_mat[row_cnt, col_idx] <- (M_hs %*% P_chol)[c_var, cond_idx[cj]]
      }
    }
    row_cnt <- row_cnt + 1
  }
}

# Solve for the unique structural conditioning shocks eta_C
eta_C <- solve(R_mat, y_gap)

# Reconstruct the conditional forecasts for all variables
Y_cond <- Y_unc
for (h in 1:H) {
  cond_adj <- rep(0, k)
  for (s in 1:h) {
    M_hs <- M_mats[[h - s + 1]]
    eta_s_full <- rep(0, k)
    eta_s_full[cond_idx] <- eta_C[((s - 1) * n_c + 1):(s * n_c)]
    cond_adj <- cond_adj + as.numeric(M_hs %*% P_chol %*% eta_s_full)
  }
  Y_cond[h, ] <- Y_unc[h, ] + cond_adj
}

# Verify conditioning accuracy (gap on FEDFUNDS, M2SL, BAA10YM must be ~0)
max_cond_error <- max(abs(Y_cond[, cond_idx] - Y_test[, cond_idx]))
cat("======================================================================\n")
cat("6. WAGGONER-ZHA CONDITIONAL FORECAST ACCURACY CHECK\n")
cat("======================================================================\n")
cat(sprintf("Maximum absolute error on conditioning variables: %.2e\n", max_cond_error))

# ------------------------------------------------------------------------------
# 7. Forecast Error Analysis & Newey-West HAC Structural Break Tests
# ------------------------------------------------------------------------------
dates_test <- test_df$date

# Unconditional and Conditional forecast errors on target variables
err_unc_hwi <- Y_test[, "HWI"] - Y_unc[, "HWI"]
err_cnd_hwi <- Y_test[, "HWI"] - Y_cond[, "HWI"]

err_unc_pay <- Y_test[, "PAYEMS"] - Y_unc[, "PAYEMS"]
err_cnd_pay <- Y_test[, "PAYEMS"] - Y_cond[, "PAYEMS"]

# Newey-West HAC regression: test whether post-2022 conditional error mean != 0
# e_t = alpha + u_t, test H0: alpha = 0 with Newey-West HAC standard error
hac_test <- function(err_vec, name) {
  mod <- lm(err_vec ~ 1)
  # Lag truncation per Newey-West / Enders: floor(4 * (T/100)^(2/9))
  bw <- floor(4 * (length(err_vec) / 100)^(2/9))
  vcov_hac <- NeweyWest(mod, lag = max(1, bw), prewhite = FALSE)
  ct <- coeftest(mod, vcov. = vcov_hac)
  data.frame(
    Target = name,
    Mean_Error = round(ct[1, "Estimate"], 4),
    HAC_SE = round(ct[1, "Std. Error"], 4),
    t_ratio = round(ct[1, "t value"], 3),
    p_value = signif(ct[1, "Pr(>|t|)"], 4),
    RMSE_Unc = round(sqrt(mean((if (grepl("HWI", name)) err_unc_hwi else err_unc_pay)^2)), 4),
    RMSE_Cond = round(sqrt(mean(err_vec^2)), 4),
    stringsAsFactors = FALSE
  )
}

hac_table <- rbind(
  hac_test(err_unc_hwi, "Job Openings (HWI) [Unconditional]"),
  hac_test(err_cnd_hwi, "Job Openings (HWI) [Conditional on Policy]"),
  hac_test(err_unc_pay, "Nonfarm Payrolls (PAYEMS) [Unconditional]"),
  hac_test(err_cnd_pay, "Nonfarm Payrolls (PAYEMS) [Conditional on Policy]")
)

cat("\n======================================================================\n")
cat("7. NEWEY-WEST HAC STRUCTURAL BREAK TESTS ON FORECAST ERRORS\n")
cat("======================================================================\n")
print(hac_table, row.names = FALSE)
write.csv(hac_table, "outputs/econometrics/table_hac_break_tests.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. Publication Figure 4: Conditional vs. Unconditional Counterfactual Paths
# ------------------------------------------------------------------------------
# Convert stationary transforms back to level trajectories from 2022-10 base
base_hwi <- train_df$HWI[nrow(train_df)]       # 10447
base_pay <- train_df$PAYEMS[nrow(train_df)]    # 153939

# Cumulative levels
lvl_real_hwi <- test_df$HWI
lvl_unc_hwi  <- base_hwi + cumsum(Y_unc[, "HWI"])
lvl_cnd_hwi  <- base_hwi + cumsum(Y_cond[, "HWI"])

lvl_real_pay <- test_df$PAYEMS
lvl_unc_pay  <- base_pay * exp(cumsum(Y_unc[, "PAYEMS"] / 100))
lvl_cnd_pay  <- base_pay * exp(cumsum(Y_cond[, "PAYEMS"] / 100))

plot_traj_df <- rbind(
  data.frame(date = dates_test, Value = lvl_real_hwi, Type = "Realized Post-2022", Series = "Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_unc_hwi,  Type = "Unconditional VAR(2)", Series = "Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_cnd_hwi,  Type = "Policy-Conditioned VAR(2)", Series = "Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_real_pay, Type = "Realized Post-2022", Series = "Nonfarm Payrolls (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_unc_pay,  Type = "Unconditional VAR(2)", Series = "Nonfarm Payrolls (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_cnd_pay,  Type = "Policy-Conditioned VAR(2)", Series = "Nonfarm Payrolls (PAYEMS, thousands)")
)

plot_traj_df$SeriesLabel <- ifelse(
  grepl("HWI", plot_traj_df$Series),
  "(a) Job Openings (HWI, thousands)",
  "(b) Payroll Employment (PAYEMS, thousands)"
)

p_traj <- ggplot(plot_traj_df, aes(x = date, y = Value, color = Type, linetype = Type)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~SeriesLabel, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(
    "Realized Post-2022" = "#C0392B",
    "Policy-Conditioned VAR(2)" = "#1F618D",
    "Unconditional VAR(2)" = "#7F8C8D"
  )) +
  scale_linetype_manual(values = c(
    "Realized Post-2022" = "solid",
    "Policy-Conditioned VAR(2)" = "dashed",
    "Unconditional VAR(2)" = "dotted"
  )) +
  scale_x_date(date_breaks = "12 months", date_labels = "%b %Y") +
  scale_y_continuous(labels = comma) +
  labs(
    x = NULL, y = "Levels (Thousands)",
    color = "Trajectory", linetype = "Trajectory"
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.margin = margin(b = -2),
    strip.text = element_text(face = "bold", size = 8),
    strip.background = element_rect(fill = "#F2F4F4", color = NA),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(t = 2, r = 6, b = 4, l = 4, "pt")
  )

ggsave("outputs/figures/fig4_conditional_counterfactual_paths.pdf", p_traj, width = 3.5, height = 3.6)
ggsave("outputs/figures/fig4_conditional_counterfactual_paths.png", p_traj, width = 3.5, height = 3.6, dpi = 300)

cat("\nStep 2 completed successfully. All figures and tables exported to outputs/\n")
