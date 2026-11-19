#!/usr/bin/env Rscript
# ==============================================================================
# 02_var_estimation_and_structural_irf.R
# Step 2: Common-Sample Lag Selection (Enders 2015, p. 304), Wei (2019, Sec. 2.7)
#         Outlier-Adjusted VAR(p) Estimation, Stability Verification, Cholesky
#         Structural IRF & Ordering Invariance Audit (Enders 2015, p. 298),
#         Waggoner-Zha (1999) Conditional Forecasting, & Newey-West HAC Tests.
#
# Methodological references:
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.), Ch. 5 (Secs. 5.5-5.8).
#   - Wei, W. W. S. (2019). Multivariate Time Series Analysis and Applications, Ch. 2 (Secs. 2.5-2.7).
#   - Waggoner, D. F. & Zha, T. (1999). Conditional Forecasts in Dynamic Multivariate Models.
#     Review of Economics and Statistics, 81(4), 639-651.
#   - Newey, W. K. & West, K. D. (1987). Econometrica, 55(3), 703-708.
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

source("src/econometrics/theme_publication.R")

dir.create("outputs/econometrics", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/assets/images", recursive = TRUE, showWarnings = FALSE)

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

df$t_FEDFUNDS <- d1(df$FEDFUNDS)        # Δ policy rate (pp)
df$t_M2SL     <- d2log1(df$M2SL)        # Δ² log(M2SL) * 100
df$t_BAA10YM  <- df$BAA10YM             # corporate bond credit spread (pp)
df$t_INDPRO   <- dlog1(df$INDPRO)       # monthly industrial output growth (%)
df$t_HWI      <- d1(df$HWI)             # monthly change in job openings (k)
df$t_PAYEMS   <- dlog1(df$PAYEMS)       # monthly payroll employment growth (%)

full_sample <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]
train_df    <- full_sample[full_sample$date <= as.Date("2022-10-01"), ] # T_train = 154
test_df     <- full_sample[full_sample$date >= as.Date("2022-11-01"), ] # T_test  = 46

var_names <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")
var_cols  <- c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")

Y_train <- as.matrix(train_df[, var_cols]); colnames(Y_train) <- var_names
Y_test  <- as.matrix(test_df[, var_cols]);  colnames(Y_test)  <- var_names

k    <- ncol(Y_train)
T_tr <- nrow(Y_train)
H    <- nrow(Y_test)

# Wei (2019, Sec. 2.7) Pulse Dummies for Detected COVID-19 Outliers (2020M3-2020M5)
covid_dum <- cbind(
  D2020M3 = as.numeric(train_df$date == as.Date("2020-03-01")),
  D2020M4 = as.numeric(train_df$date == as.Date("2020-04-01")),
  D2020M5 = as.numeric(train_df$date == as.Date("2020-05-01"))
)
n_dum <- ncol(covid_dum)

cat("======================================================================\n")
cat("1. ESTIMATION SYSTEM SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Variables (k = %d)      : %s\n", k, paste(var_names, collapse = ", ")))
cat(sprintf("Pre-AI Training (T = %d): %s to %s\n", T_tr, min(train_df$date), max(train_df$date)))
cat(sprintf("Post-AI Test (H = %d)   : %s to %s\n", H, min(test_df$date), max(test_df$date)))
cat("Outlier Intervention    : Wei (2019, Sec. 2.7) pulse dummies for 2020M3-2020M5\n\n")

# ------------------------------------------------------------------------------
# 2. Lag Length Selection on Exact Common Effective Sample (Enders 2015, p. 304)
#    Every candidate VAR(p) is estimated over the IDENTICAL T_eff = T_tr - max_p
#    observations by slicing Y_train[(max_p - p + 1):T_tr, ]!
# ------------------------------------------------------------------------------
max_p <- 8
T_eff_lr <- T_tr - max_p

log_det_sigma <- numeric(max_p)
aic_vec <- numeric(max_p)
bic_vec <- numeric(max_p)
hq_vec  <- numeric(max_p)

for (p in 1:max_p) {
  # Exact common estimation window: dependent variable always rows (max_p + 1):T_tr
  idx_slice <- (max_p - p + 1):T_tr
  v_fit <- VAR(Y_train[idx_slice, ], p = p, type = "const", exogen = covid_dum[idx_slice, ])
  res_p <- residuals(v_fit) # nrow = T_eff_lr = 146 identically for all p
  sig_p <- crossprod(res_p) / T_eff_lr
  log_det_sigma[p] <- log(det(sig_p))

  # Total parameters per equation: 1 (const) + n_dum + k*p
  num_params <- 1 + n_dum + k * p
  aic_vec[p] <- log_det_sigma[p] + (2 * k * num_params) / T_eff_lr
  bic_vec[p] <- log_det_sigma[p] + (log(T_eff_lr) * k * num_params) / T_eff_lr
  hq_vec[p]  <- log_det_sigma[p] + (2 * log(log(T_eff_lr)) * k * num_params) / T_eff_lr
}

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

# Sims (1980) small-sample LR correction (Enders 2015, p. 304, Eq. 5.56):
# (T - c) * (ln|Sigma_r| - ln|Sigma_u|) ~ chi^2(k^2), c = 1 + n_dum + k * p_u
for (p in 1:(max_p - 1)) {
  p_u <- p + 1
  c_sims <- 1 + n_dum + k * p_u
  df_lr <- k^2
  lr_stat <- (T_eff_lr - c_sims) * (log_det_sigma[p] - log_det_sigma[p_u])
  lr_table$LR_Stat_Sims[p] <- round(lr_stat, 2)
  lr_table$LR_df[p] <- df_lr
  lr_table$LR_pvalue[p] <- signif(1 - pchisq(lr_stat, df = df_lr), 3)
}

best_aic <- which.min(aic_vec)
best_bic <- which.min(bic_vec)
best_hq  <- which.min(hq_vec)

cat("======================================================================\n")
cat("2. COMMON-SAMPLE LAG LENGTH SELECTION (ENDERS 2015, p. 304)\n")
cat("======================================================================\n")
print(lr_table, row.names = FALSE)
cat(sprintf("Selected lag (common T_eff = %d): AIC = %d | BIC = %d | HQ = %d\n",
            T_eff_lr, best_aic, best_bic, best_hq))
write.csv(lr_table, "outputs/econometrics/table_var_lag_selection.csv", row.names = FALSE)

# Select p_opt = 2 (parsimonious BIC/HQ-aligned benchmark for monthly dynamics)
p_opt <- max(2, best_bic)
cat(sprintf("Selected baseline specification: COVID-Controlled VAR(%d)\n\n", p_opt))

# ------------------------------------------------------------------------------
# 3. Baseline Outlier-Adjusted VAR(2) Estimation, Stability & Diagnostics
#    (Wei 2019, Secs. 2.5.3 & 2.7; Enders 2015, Ch. 5)
# ------------------------------------------------------------------------------
var_baseline <- VAR(Y_train, p = p_opt, type = "const", exogen = covid_dum)
var_uncontrolled <- VAR(Y_train, p = p_opt, type = "const") # Retained for contrast

roots_obj <- roots(var_baseline, modulus = TRUE)
max_eigen <- max(roots_obj)
is_stable <- all(roots_obj < 1.0)

cat("======================================================================\n")
cat("3A. STABILITY VERIFICATION (COMPANION MATRIX ROOTS)\n")
cat("======================================================================\n")
cat(sprintf("Max eigenvalue modulus: %.4f (< 1.0: %s)\n", max_eigen, ifelse(is_stable, "STABLE", "UNSTABLE")))
print(round(roots_obj, 4))

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
# 4. Residual Correlation Matrix & Cholesky Ordering Audit (Enders 2015, p. 298)
# ------------------------------------------------------------------------------
E_opt <- residuals(var_baseline)
T_e   <- nrow(E_opt)
df_eq <- T_e - (1 + n_dum + k * p_opt)
Sigma_res <- crossprod(E_opt) / df_eq # Exact OLS finite-sample Sigma matching vars::irf
P_res <- cov2cor(Sigma_res)
rownames(P_res) <- colnames(P_res) <- var_names

se_rho <- 1 / sqrt(T_e)
crit_rho_05 <- 1.96 * se_rho
crit_rho_01 <- 2.576 * se_rho

cat("\n======================================================================\n")
cat("4. RESIDUAL CORRELATION MATRIX & CHOLESKY ORDERING AUDIT (ENDERS 2015, p. 298)\n")
cat("======================================================================\n")
cat(sprintf("Sample T_e = %d | SE(rho) = %.4f | 5%% Threshold = +/-%.3f | 1%% Threshold = +/-%.3f\n",
            T_e, se_rho, crit_rho_05, crit_rho_01))
print(round(P_res, 3))
write.csv(round(P_res, 4), "outputs/econometrics/table_residual_correlation_matrix.csv")

# ------------------------------------------------------------------------------
# 5. Structural Cholesky IRFs & Ordering Invariance Robustness
#    (Enders 2015, Sec. 5.7, pp. 297-302)
# ------------------------------------------------------------------------------
irf_horizon <- 24
n_boot <- 500

cat("\nComputing COVID-controlled Cholesky IRFs (500 bootstrap reps, horizon = 24)...\n")
set.seed(563)
irf_primary <- irf(var_baseline, impulse = c("FEDFUNDS", "BAA10YM"),
                   response = c("HWI", "PAYEMS", "INDPRO"),
                   n.ahead = irf_horizon, boot = TRUE, runs = n_boot, ci = 0.95)

# Alternative Ordering (Labor-First: PAYEMS, HWI, INDPRO, FEDFUNDS, M2SL, BAA10YM)
alt1_names <- c("PAYEMS", "HWI", "INDPRO", "FEDFUNDS", "M2SL", "BAA10YM")
var_alt1 <- VAR(Y_train[, alt1_names], p = p_opt, type = "const", exogen = covid_dum)
set.seed(563)
irf_alt1 <- irf(var_alt1, impulse = c("FEDFUNDS", "BAA10YM"),
                response = c("HWI", "PAYEMS", "INDPRO"),
                n.ahead = irf_horizon, boot = TRUE, runs = n_boot, ci = 0.95)

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

irf_plot_df <- rbind(
  tidy_irf(irf_primary, "FEDFUNDS", "HWI",    "Primary Ordering"),
  tidy_irf(irf_primary, "FEDFUNDS", "PAYEMS", "Primary Ordering"),
  tidy_irf(irf_primary, "BAA10YM",  "HWI",    "Primary Ordering"),
  tidy_irf(irf_primary, "BAA10YM",  "PAYEMS", "Primary Ordering")
)
irf_plot_df$PanelLabel <- factor(
  paste0(irf_plot_df$Impulse, " -> ", irf_plot_df$Response),
  levels = c("FEDFUNDS -> HWI", "FEDFUNDS -> PAYEMS", "BAA10YM -> HWI", "BAA10YM -> PAYEMS"),
  labels = c("(a) FEDFUNDS -> HWI", "(b) FEDFUNDS -> PAYEMS", "(c) BAA10YM -> HWI", "(d) BAA10YM -> PAYEMS")
)

# Unified Figure 2 builder (Print 3.5 in & Web 4.4 in — ZERO strip clipping!)
build_fig2 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  ggplot(irf_plot_df, aes(x = Horizon, y = Point)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = PAL$band68, alpha = 0.28) +
    geom_line(color = PAL$navy, linewidth = if (is_web) 0.85 else 0.65) +
    geom_hline(yintercept = 0, linetype = "dashed", color = PAL$gray, linewidth = 0.35) +
    facet_wrap(~PanelLabel, scales = "free_y", ncol = 2) +
    labs(x = "Forecast Horizon (Months)", y = "Orthogonalized Response") +
    th +
    theme(
      strip.text = element_text(face = "bold", size = if (is_web) 9.0 else 6.2,
                                margin = margin(t = 2, b = 2, l = 1, r = 1))
    )
}

p_irf_print <- build_fig2(is_web = FALSE)
ggsave("outputs/figures/fig2_structural_irfs.pdf", p_irf_print, width = 3.5, height = 2.8)
ggsave("outputs/figures/fig2_structural_irfs.png", p_irf_print, width = 3.5, height = 2.8, dpi = 300)

p_irf_web <- build_fig2(is_web = TRUE)
ggsave("outputs/figures/fig2_web_structural_irfs.png", p_irf_web, width = 4.6, height = 3.6, dpi = 300)
file.copy("outputs/figures/fig2_web_structural_irfs.png", "docs/assets/images/fig2_structural_irfs.png", overwrite = TRUE)

# Unified Figure 3 builder (Ordering Robustness — Print 3.5 in & Web 4.6 in)
comp_df <- rbind(
  cbind(tidy_irf(irf_primary, "FEDFUNDS", "HWI"),    Ordering = "Baseline (Policy First)"),
  cbind(tidy_irf(irf_alt1,    "FEDFUNDS", "HWI"),    Ordering = "Alternative (Labor First)"),
  cbind(tidy_irf(irf_primary, "FEDFUNDS", "PAYEMS"), Ordering = "Baseline (Policy First)"),
  cbind(tidy_irf(irf_alt1,    "FEDFUNDS", "PAYEMS"), Ordering = "Alternative (Labor First)")
)
comp_df$ResponseVar <- ifelse(comp_df$Response == "HWI", "(a) Job Openings (HWI)", "(b) Payrolls (PAYEMS)")

build_fig3 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  ggplot(comp_df, aes(x = Horizon, y = Point, color = Ordering, linetype = Ordering)) +
    geom_line(linewidth = if (is_web) 0.85 else 0.65) +
    geom_hline(yintercept = 0, linetype = "dashed", color = PAL$gray, linewidth = 0.35) +
    facet_wrap(~ResponseVar, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Baseline (Policy First)" = PAL$navy, "Alternative (Labor First)" = PAL$red)) +
    scale_linetype_manual(values = c("Baseline (Policy First)" = "solid", "Alternative (Labor First)" = "dashed")) +
    guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
    labs(x = "Forecast Horizon (Months)", y = "Impulse Response", color = NULL, linetype = NULL) +
    th +
    theme(
      strip.text = element_text(face = "bold", size = if (is_web) 9.0 else 6.4),
      legend.text = element_text(size = if (is_web) 8.8 else 6.2)
    )
}

p_rob_print <- build_fig3(is_web = FALSE)
ggsave("outputs/figures/fig3_ordering_robustness_irf.pdf", p_rob_print, width = 3.5, height = 2.4)
ggsave("outputs/figures/fig3_ordering_robustness_irf.png", p_rob_print, width = 3.5, height = 2.4, dpi = 300)

p_rob_web <- build_fig3(is_web = TRUE)
ggsave("outputs/figures/fig3_web_ordering_robustness.png", p_rob_web, width = 4.6, height = 3.6, dpi = 300)
file.copy("outputs/figures/fig3_web_ordering_robustness.png", "docs/assets/images/fig3_ordering_robustness_irf.png", overwrite = TRUE)

# ------------------------------------------------------------------------------
# 6. Waggoner-Zha (1999) Conditional Forecasting Engine
#    Implements BOTH:
#      (a) Exact Minimum-MSE All-Shocks Projection (Waggoner & Zha 1999, Eq. 11;
#          Banbura et al. 2015) — strictly invariant to Cholesky ordering!
#      (b) Structural-Policy-Shocks-Only Projection (setting non-policy shocks = 0)
# ------------------------------------------------------------------------------
cond_idx <- c(1, 2, 3) # FEDFUNDS, M2SL, BAA10YM
hist_init <- Y_train[(T_tr - p_opt + 1):T_tr, ]

compute_wz <- function(v_obj, Y_init, Y_eval, cond_idx, mode = c("min_mse", "structural")) {
  mode <- match.arg(mode)
  A_est <- Acoef(v_obj)
  c_est <- sapply(v_obj$varresult, function(eq) coef(eq)["const"])
  E_v   <- residuals(v_obj)
  n_par <- length(coef(v_obj$varresult[[1]]))
  Sig_v <- crossprod(E_v) / (nrow(E_v) - n_par)
  P_c   <- t(chol(Sig_v))
  p_v   <- length(A_est)
  k_v   <- ncol(Y_eval)
  H_v   <- nrow(Y_eval)
  n_c   <- length(cond_idx)

  Y_unc <- matrix(0, H_v, k_v); colnames(Y_unc) <- colnames(Y_eval)
  for (h in 1:H_v) {
    pr <- c_est
    for (j in 1:p_v) {
      lv <- if ((h - j) <= 0) Y_init[p_v + (h - j), ] else Y_unc[h - j, ]
      pr <- pr + A_est[[j]] %*% lv
    }
    Y_unc[h, ] <- pr
  }

  M_mats <- vector("list", H_v); M_mats[[1]] <- diag(k_v)
  for (m in 2:H_v) {
    Mm <- matrix(0, k_v, k_v)
    for (j in 1:min(p_v, m - 1)) Mm <- Mm + A_est[[j]] %*% M_mats[[m - j]]
    M_mats[[m]] <- Mm
  }

  if (mode == "min_mse") {
    R <- matrix(0, H_v * n_c, H_v * k_v); r <- numeric(H_v * n_c); rc <- 1
    for (h in 1:H_v) {
      for (ci in 1:n_c) {
        cv <- cond_idx[ci]
        r[rc] <- Y_eval[h, cv] - Y_unc[h, cv]
        for (s in 1:h) R[rc, ((s - 1) * k_v + 1):(s * k_v)] <- (M_mats[[h - s + 1]] %*% P_c)[cv, ]
        rc <- rc + 1
      }
    }
    eta_all <- as.numeric(t(R) %*% solve(R %*% t(R), r))
    Y_cond <- Y_unc
    for (h in 1:H_v) {
      adj <- rep(0, k_v)
      for (s in 1:h) adj <- adj + as.numeric(M_mats[[h - s + 1]] %*% P_c %*% eta_all[((s - 1) * k_v + 1):(s * k_v)])
      Y_cond[h, ] <- Y_unc[h, ] + adj
    }
  } else {
    R <- matrix(0, H_v * n_c, H_v * n_c); r <- numeric(H_v * n_c); rc <- 1
    for (h in 1:H_v) {
      for (ci in 1:n_c) {
        cv <- cond_idx[ci]
        r[rc] <- Y_eval[h, cv] - Y_unc[h, cv]
        for (s in 1:h) {
          MP <- M_mats[[h - s + 1]] %*% P_c
          for (cj in 1:n_c) R[rc, (s - 1) * n_c + cj] <- MP[cv, cond_idx[cj]]
        }
        rc <- rc + 1
      }
    }
    eta_c <- solve(R, r)
    Y_cond <- Y_unc
    for (h in 1:H_v) {
      adj <- rep(0, k_v)
      for (s in 1:h) {
        es <- rep(0, k_v); es[cond_idx] <- eta_c[((s - 1) * n_c + 1):(s * n_c)]
        adj <- adj + as.numeric(M_mats[[h - s + 1]] %*% P_c %*% es)
      }
      Y_cond[h, ] <- Y_unc[h, ] + adj
    }
  }
  colnames(Y_cond) <- colnames(Y_eval)
  list(Y_unc = Y_unc, Y_cond = Y_cond)
}

wz_base_mmse   <- compute_wz(var_baseline,     hist_init, Y_test, cond_idx, mode = "min_mse")
wz_base_struct <- compute_wz(var_baseline,     hist_init, Y_test, cond_idx, mode = "structural")
wz_uncov_mmse  <- compute_wz(var_uncontrolled, hist_init, Y_test, cond_idx, mode = "min_mse")

Y_unc  <- wz_base_mmse$Y_unc
Y_cond <- wz_base_mmse$Y_cond

max_cond_error <- max(abs(Y_cond[, cond_idx] - Y_test[, cond_idx]))
cat("======================================================================\n")
cat("6. WAGGONER-ZHA CONDITIONAL FORECAST ACCURACY CHECK\n")
cat("======================================================================\n")
cat(sprintf("Maximum absolute error on conditioning variables: %.2e\n", max_cond_error))

# ------------------------------------------------------------------------------
# 7. Newey-West HAC Structural Break Tests (Newey & West 1987; bandwidth b = 3)
# ------------------------------------------------------------------------------
dates_test <- test_df$date

err_unc_hwi <- Y_test[, "HWI"]    - Y_unc[, "HWI"]
err_cnd_hwi <- Y_test[, "HWI"]    - Y_cond[, "HWI"]
err_unc_pay <- Y_test[, "PAYEMS"] - Y_unc[, "PAYEMS"]
err_cnd_pay <- Y_test[, "PAYEMS"] - Y_cond[, "PAYEMS"]

hac_test <- function(err_vec, ref_unc, name) {
  mod <- lm(err_vec ~ 1)
  bw <- floor(4 * (length(err_vec) / 100)^(2/9)) # b = 3 for H = 46
  vcov_hac <- NeweyWest(mod, lag = max(1, bw), prewhite = FALSE)
  ct <- coeftest(mod, vcov. = vcov_hac)
  data.frame(
    Target     = name,
    Mean_Error = round(ct[1, "Estimate"], 4),
    HAC_SE     = round(ct[1, "Std. Error"], 4),
    t_ratio    = round(ct[1, "t value"], 3),
    p_value    = signif(ct[1, "Pr(>|t|)"], 4),
    RMSE_Unc   = round(sqrt(mean(ref_unc^2)), 4),
    RMSE_Cond  = round(sqrt(mean(err_vec^2)), 4),
    stringsAsFactors = FALSE
  )
}

hac_table <- rbind(
  hac_test(err_unc_hwi, err_unc_hwi, "Job Openings (HWI) [Unconditional]"),
  hac_test(err_cnd_hwi, err_unc_hwi, "Job Openings (HWI) [Policy-Conditioned]"),
  hac_test(err_unc_pay, err_unc_pay, "Nonfarm Payrolls (PAYEMS) [Unconditional]"),
  hac_test(err_cnd_pay, err_unc_pay, "Nonfarm Payrolls (PAYEMS) [Policy-Conditioned]")
)

cat("\n======================================================================\n")
cat("7. NEWEY-WEST HAC STRUCTURAL BREAK TESTS (COVID-CONTROLLED VAR(2))\n")
cat("======================================================================\n")
print(hac_table, row.names = FALSE)
write.csv(hac_table, "outputs/econometrics/table_hac_break_tests.csv", row.names = FALSE)

# Format p-values cleanly for LaTeX
fmt_p <- function(p) {
  if (p < 0.001) {
    ex <- floor(log10(p))
    sprintf("%.2f \\times 10^{%d}", p / 10^ex, ex)
  } else {
    sprintf("%.3f", p)
  }
}

tex_tab4 <- sprintf(
"\\begin{table}[t]
\\centering
\\caption{Newey--West HAC Structural Break Tests on Post-2022 Errors}
\\label{tab:hac_tests}
\\begin{minipage}{\\columnwidth}
\\centering
\\fontsize{7.0pt}{8.4pt}\\selectfont
\\setlength{\\tabcolsep}{2.3pt}
\\begin{tabular}{lccccc}
\\toprule
\\textbf{Target \\& Regime} & $\\bar{e}$ & \\textbf{HAC SE} & $t$-\\textbf{ratio} & $p$-\\textbf{value} & \\textbf{RMSE} \\\\
\\midrule
\\multicolumn{6}{l}{\\textit{Job Openings (\\texttt{HWI}, Thousands of Openings/Month)}} \\\\
Unconditional VAR(2) & $%.2f$ & $%.2f$ & $%.2f$ & $%s$ & $%.2f$ \\\\
Policy-Conditioned & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%s}$ & $%.2f$ \\\\
\\midrule
\\multicolumn{6}{l}{\\textit{Nonfarm Payrolls (\\texttt{PAYEMS}, Monthly \\%% Growth)}} \\\\
Unconditional VAR(2) & $%+.3f$ & $%.3f$ & $%.2f$ & $%s$ & $%.3f$ \\\\
Policy-Conditioned & $\\mathbf{%+.3f}$ & $\\mathbf{%.3f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%s}$ & $%.3f$ \\\\
\\bottomrule
\\end{tabular}
\\vskip 2pt
{\\scriptsize \\textit{Notes}: COVID-controlled stationary VAR(2) (2020M3--2020M5 pulse dummies per Wei 2019, Sec.~2.7) evaluated over 2022M11--2026M8 ($H = 46$). Conditioning block $\\mathcal{C} = \\{\\text{FEDFUNDS}, \\text{M2SL}, \\text{BAA10YM}\\}$ via exact Waggoner--Zha (1999) projection. HAC SEs use Newey--West (1987) Bartlett kernel (bandwidth $b = 3$).\\par}
\\end{minipage}
\\end{table}
",
hac_table$Mean_Error[1], hac_table$HAC_SE[1], hac_table$t_ratio[1], fmt_p(hac_table$p_value[1]), hac_table$RMSE_Cond[1],
hac_table$Mean_Error[2], hac_table$HAC_SE[2], hac_table$t_ratio[2], fmt_p(hac_table$p_value[2]), hac_table$RMSE_Cond[2],
hac_table$Mean_Error[3], hac_table$HAC_SE[3], hac_table$t_ratio[3], fmt_p(hac_table$p_value[3]), hac_table$RMSE_Cond[3],
hac_table$Mean_Error[4], hac_table$HAC_SE[4], hac_table$t_ratio[4], fmt_p(hac_table$p_value[4]), hac_table$RMSE_Cond[4]
)
cat(tex_tab4, file = "outputs/econometrics/table4_hac_break_tests.tex")

# ------------------------------------------------------------------------------
# 8. Publication Figure 4: Conditional Counterfactual Levels (Unified Punto!)
#    Zero legend clipping: 2-row wrapped legend on 3.5 in print; 2-col widescreen web
# ------------------------------------------------------------------------------
base_hwi <- train_df$HWI[T_tr]
base_pay <- train_df$PAYEMS[T_tr]

lvl_real_hwi <- test_df$HWI
lvl_unc_hwi  <- base_hwi + cumsum(Y_unc[, "HWI"])
lvl_cnd_hwi  <- base_hwi + cumsum(Y_cond[, "HWI"])

lvl_real_pay <- test_df$PAYEMS
lvl_unc_pay  <- base_pay * exp(cumsum(Y_unc[, "PAYEMS"] / 100))
lvl_cnd_pay  <- base_pay * exp(cumsum(Y_cond[, "PAYEMS"] / 100))

cat(sprintf("\nCumulative Level Outcomes (2026M8):\n"))
cat(sprintf("  HWI   : Realized = %.0fk | Conditional = %.0fk (Gap = %.0fk, %.1f%%) | Unconditional = %.0fk\n",
            lvl_real_hwi[H], lvl_cnd_hwi[H], lvl_real_hwi[H] - lvl_cnd_hwi[H],
            100 * (lvl_real_hwi[H] / lvl_cnd_hwi[H] - 1), lvl_unc_hwi[H]))
cat(sprintf("  PAYEMS: Realized = %.0fk | Conditional = %.0fk (Gap = %.0fk, %.1f%%) | Unconditional = %.0fk\n",
            lvl_real_pay[H], lvl_cnd_pay[H], lvl_real_pay[H] - lvl_cnd_pay[H],
            100 * (lvl_real_pay[H] / lvl_cnd_pay[H] - 1), lvl_unc_pay[H]))

plot_traj_df <- rbind(
  data.frame(date = dates_test, Value = lvl_real_hwi, Type = "Realized Post-2022",       SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_cnd_hwi,  Type = "Policy-Conditioned VAR(2)", SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_unc_hwi,  Type = "Unconditional VAR(2)",     SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_real_pay, Type = "Realized Post-2022",       SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_cnd_pay,  Type = "Policy-Conditioned VAR(2)", SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_unc_pay,  Type = "Unconditional VAR(2)",     SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)")
)
plot_traj_df$Type <- factor(plot_traj_df$Type, levels = c("Realized Post-2022", "Policy-Conditioned VAR(2)", "Unconditional VAR(2)"))

build_fig4 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  ggplot(plot_traj_df, aes(x = date, y = Value, color = Type, linetype = Type)) +
    geom_line(linewidth = if (is_web) 0.95 else 0.70) +
    facet_wrap(~SeriesLabel, scales = "free_y", ncol = if (is_web) 2 else 1) +
    scale_color_manual(values = c(
      "Realized Post-2022"        = PAL$red,
      "Policy-Conditioned VAR(2)" = PAL$navy,
      "Unconditional VAR(2)"      = PAL$gray
    )) +
    scale_linetype_manual(values = c(
      "Realized Post-2022"        = "solid",
      "Policy-Conditioned VAR(2)" = "dashed",
      "Unconditional VAR(2)"      = "dotted"
    )) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    guides(
      color    = guide_legend(nrow = if (is_web) 1 else 2, byrow = TRUE),
      linetype = guide_legend(nrow = if (is_web) 1 else 2, byrow = TRUE)
    ) +
    labs(x = NULL, y = "Levels (Thousands)", color = NULL, linetype = NULL) +
    th
}

p_traj_print <- build_fig4(is_web = FALSE)
ggsave("outputs/figures/fig4_conditional_counterfactual_paths.pdf", p_traj_print, width = 3.5, height = 3.5)
ggsave("outputs/figures/fig4_conditional_counterfactual_paths.png", p_traj_print, width = 3.5, height = 3.5, dpi = 300)

p_traj_web <- build_fig4(is_web = TRUE)
ggsave("outputs/figures/fig4_web_counterfactuals.png", p_traj_web, width = 8.4, height = 3.9, dpi = 300)
file.copy("outputs/figures/fig4_web_counterfactuals.png", "docs/assets/images/fig4_web_counterfactuals.png", overwrite = TRUE)

write.csv(plot_traj_df, "outputs/econometrics/table_counterfactual_trajectories.csv", row.names = FALSE)
cat("\nStep 2 completed successfully. All unified figures and tables exported.\n")
