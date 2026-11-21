#!/usr/bin/env Rscript
# ==============================================================================
# 02_var_estimation_and_structural_irf.R
# Step 2: Common-Sample Lag Selection (Enders 2015, p. 304), COVID-Controlled
#         VAR(2) & Formal Cointegrated VECM(r=2) with 2020:03-2020:12 Dummies
#         (Wei 2019 Sec. 2.7; Enders 2015 Ch. 5 & 6), Cholesky Structural IRFs,
#         Waggoner-Zha (1999) Conditional Counterfactuals, & Newey-West HAC Tests.
#
# Methodological references:
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.), Ch. 5 & 6.
#   - Wei, W. W. S. (2019). Multivariate Time Series Analysis and Applications, Ch. 2.
#   - Waggoner, D. F. & Zha, T. (1999). Review of Economics and Statistics, 81(4), 639-651.
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
# 1. Load FRED-MD Data & Prepare Stationary + Levels Estimation Systems
#    Note: Since the 2017-09 FRED-MD vintage, HWI directly reports BLS JOLTS
#    Total Nonfarm Job Openings (JTSJOL) over 2000M12-present (not a proxy).
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

d1     <- function(x) c(NA, diff(x))
dlog1  <- function(x) c(NA, diff(log(x))) * 100
d2log1 <- function(x) c(NA, NA, diff(diff(log(x)))) * 100

# Stationary FRED-MD transforms
df$t_FEDFUNDS <- d1(df$FEDFUNDS)        # tcode=2: Δ policy rate (pp)
df$t_M2SL     <- d2log1(df$M2SL)        # tcode=6: Δ² log(M2SL) * 100
df$t_BAA10YM  <- df$BAA10YM             # tcode=1: corporate credit spread (pp)
df$t_INDPRO   <- dlog1(df$INDPRO)       # tcode=5: Δ log(INDPRO) * 100 (%)
df$t_HWI      <- d1(df$HWI)             # tcode=2: Δ HWI (thousands of openings)
df$t_PAYEMS   <- dlog1(df$PAYEMS)       # tcode=5: Δ log(PAYEMS) * 100 (%)

# Levels / 100*log-levels system for formal VECM(r=2)
df$L_FEDFUNDS <- df$FEDFUNDS
df$L_M2SL     <- 100 * log(df$M2SL)
df$L_BAA10YM  <- df$BAA10YM
df$L_INDPRO   <- 100 * log(df$INDPRO)
df$L_HWI      <- 100 * log(df$HWI)
df$L_PAYEMS   <- 100 * log(df$PAYEMS)

full_sample <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]
train_df    <- full_sample[full_sample$date <= as.Date("2022-10-01"), ] # T_train = 154
test_df     <- full_sample[full_sample$date >= as.Date("2022-11-01"), ] # T_test  = 46

var_names <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")
var_cols  <- c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")
lvl_cols  <- c("L_FEDFUNDS", "L_M2SL", "L_BAA10YM", "L_INDPRO", "L_HWI", "L_PAYEMS")

Y_train     <- as.matrix(train_df[, var_cols]); colnames(Y_train)     <- var_names
Y_test      <- as.matrix(test_df[, var_cols]);  colnames(Y_test)      <- var_names
Y_train_lvl <- as.matrix(train_df[, lvl_cols]); colnames(Y_train_lvl) <- var_names
Y_test_lvl  <- as.matrix(test_df[, lvl_cols]);  colnames(Y_test_lvl)  <- var_names

k    <- ncol(Y_train)
T_tr <- nrow(Y_train)
H    <- nrow(Y_test)

# Full 2020:03-2020:12 COVID-19 Pulse Dummies (10 months, per Prof. Ceylan feedback)
covid_months <- seq(as.Date("2020-03-01"), as.Date("2020-12-01"), by = "month")
covid_dum <- sapply(covid_months, function(d) as.numeric(train_df$date == d))
colnames(covid_dum) <- format(covid_months, "D%YM%m")
n_dum <- ncol(covid_dum)

cat("======================================================================\n")
cat("1. ESTIMATION SYSTEM SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Variables (k = %d)      : %s\n", k, paste(var_names, collapse = ", ")))
cat(sprintf("Pre-AI Training (T = %d): %s to %s\n", T_tr, min(train_df$date), max(train_df$date)))
cat(sprintf("Post-AI Test (H = %d)   : %s to %s\n", H, min(test_df$date), max(test_df$date)))
cat("Frequentist COVID Spec  : 2020:03-2020:12 pulse dummies (10 months)\n\n")

# ------------------------------------------------------------------------------
# 2. Lag Length Selection on Exact Common Effective Sample (Enders 2015, p. 304)
# ------------------------------------------------------------------------------
max_p <- 8
T_eff_lr <- T_tr - max_p

log_det_sigma <- numeric(max_p)
aic_vec <- numeric(max_p)
bic_vec <- numeric(max_p)
hq_vec  <- numeric(max_p)

for (p in 1:max_p) {
  idx_slice <- (max_p - p + 1):T_tr
  v_fit <- VAR(Y_train[idx_slice, ], p = p, type = "const", exogen = covid_dum[idx_slice, ])
  res_p <- residuals(v_fit)
  sig_p <- crossprod(res_p) / T_eff_lr
  log_det_sigma[p] <- log(det(sig_p))

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

p_opt <- max(2, best_bic)
cat(sprintf("Selected baseline specification: COVID-Controlled VAR(%d) [2020M3-2020M12 dummies]\n\n", p_opt))

# ------------------------------------------------------------------------------
# 3. Baseline COVID-Controlled VAR(2) + Formal Cointegrated VECM(r=2)
#    (Wei 2019 Secs. 2.5.3 & 2.7; Enders 2015 Ch. 5 & Ch. 6)
# ------------------------------------------------------------------------------
var_baseline <- VAR(Y_train, p = p_opt, type = "const", exogen = covid_dum)

roots_obj <- roots(var_baseline, modulus = TRUE)
max_eigen <- max(roots_obj)
is_stable <- all(roots_obj < 1.0)

cat("======================================================================\n")
cat("3A. STABILITY VERIFICATION (COMPANION MATRIX ROOTS)\n")
cat("======================================================================\n")
cat(sprintf("VAR(2) Max eigenvalue modulus: %.4f (< 1.0: %s)\n", max_eigen, ifelse(is_stable, "STABLE", "UNSTABLE")))
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

# Formal Cointegrated VECM(r=2) (Enders 2015, Ch. 6; Case 3 unrestricted constant, K=3, 2020M3-12 dummies)
jo_vecm <- ca.jo(Y_train_lvl, type = "trace", ecdet = "none", K = 3, spec = "transitory", dumvar = covid_dum)
vecm_r2 <- vec2var(jo_vecm, r = 2)
cat("\n======================================================================\n")
cat("3C. FORMAL COINTEGRATED VECM(r = 2) (ENDERS 2015, CH. 6)\n")
cat("======================================================================\n")
cat("Normalized Cointegrating Vectors (beta):\n")
beta_norm <- round(jo_vecm@V[, 1:2], 4)
rownames(beta_norm) <- var_names; colnames(beta_norm) <- c("ECT1", "ECT2")
print(beta_norm)
write.csv(beta_norm, "outputs/econometrics/table_vecm_cointegrating_vectors.csv")

# ------------------------------------------------------------------------------
# 4. Residual Correlation Matrix & Cholesky Ordering Audit (Enders 2015, p. 298)
# ------------------------------------------------------------------------------
E_opt <- residuals(var_baseline)
T_e   <- nrow(E_opt)
df_eq <- T_e - (1 + n_dum + k * p_opt)
Sigma_res <- crossprod(E_opt) / df_eq
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
# ------------------------------------------------------------------------------
irf_horizon <- 24
n_boot <- 500

cat("\nComputing COVID-controlled Cholesky IRFs (500 bootstrap reps, horizon = 24)...\n")
set.seed(563)
irf_primary <- irf(var_baseline, impulse = c("FEDFUNDS", "BAA10YM"),
                   response = c("HWI", "PAYEMS", "INDPRO"),
                   n.ahead = irf_horizon, boot = TRUE, runs = n_boot, ci = 0.95)

alt1_names <- c("PAYEMS", "HWI", "INDPRO", "FEDFUNDS", "M2SL", "BAA10YM")
var_alt1 <- VAR(Y_train[, alt1_names], p = p_opt, type = "const", exogen = covid_dum)
set.seed(563)
irf_alt1 <- irf(var_alt1, impulse = c("FEDFUNDS", "BAA10YM"),
                response = c("HWI", "PAYEMS", "INDPRO"),
                n.ahead = irf_horizon, boot = TRUE, runs = n_boot, ci = 0.95)

tidy_irf <- function(irf_obj, imp, resp, label = "Baseline") {
  h <- 0:(nrow(irf_obj$irf[[imp]]) - 1)
  data.frame(
    Horizon = h, Impulse = imp, Response = resp,
    Point = as.numeric(irf_obj$irf[[imp]][, resp]),
    Lower = as.numeric(irf_obj$Lower[[imp]][, resp]),
    Upper = as.numeric(irf_obj$Upper[[imp]][, resp]),
    Model = label, stringsAsFactors = FALSE
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

build_fig2 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  ggplot(irf_plot_df, aes(x = Horizon, y = Point)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = PAL$band68, alpha = 0.28) +
    geom_line(color = PAL$navy, linewidth = if (is_web) 0.85 else 0.65) +
    geom_hline(yintercept = 0, linetype = "dashed", color = PAL$gray, linewidth = 0.35) +
    facet_wrap(~PanelLabel, scales = "free_y", ncol = 2) +
    labs(x = "Forecast Horizon (Months)", y = "Orthogonalized Response") +
    th +
    theme(strip.text = element_text(face = "bold", size = if (is_web) 9.0 else 6.2,
                                    margin = margin(t = 2, b = 2, l = 1, r = 1)))
}

p_irf_print <- build_fig2(is_web = FALSE)
ggsave("outputs/figures/fig2_structural_irfs.pdf", p_irf_print, width = 3.5, height = 2.8)
ggsave("outputs/figures/fig2_structural_irfs.png", p_irf_print, width = 3.5, height = 2.8, dpi = 300)

p_irf_web <- build_fig2(is_web = TRUE)
ggsave("outputs/figures/fig2_web_structural_irfs.png", p_irf_web, width = 4.6, height = 3.6, dpi = 300)
file.copy("outputs/figures/fig2_web_structural_irfs.png", "docs/assets/images/fig2_structural_irfs.png", overwrite = TRUE)

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
    theme(strip.text = element_text(face = "bold", size = if (is_web) 9.0 else 6.4),
          legend.text = element_text(size = if (is_web) 8.8 else 6.2))
}

p_rob_print <- build_fig3(is_web = FALSE)
ggsave("outputs/figures/fig3_ordering_robustness_irf.pdf", p_rob_print, width = 3.5, height = 2.4)
ggsave("outputs/figures/fig3_ordering_robustness_irf.png", p_rob_print, width = 3.5, height = 2.4, dpi = 300)

p_rob_web <- build_fig3(is_web = TRUE)
ggsave("outputs/figures/fig3_web_ordering_robustness.png", p_rob_web, width = 4.6, height = 3.6, dpi = 300)
file.copy("outputs/figures/fig3_web_ordering_robustness.png", "docs/assets/images/fig3_ordering_robustness_irf.png", overwrite = TRUE)

# ------------------------------------------------------------------------------
# 6. Exact Waggoner-Zha (1999) Conditional Forecasting (VAR(2) & VECM(r=2))
# ------------------------------------------------------------------------------
cond_idx <- c(1, 2, 3)

compute_wz_generic <- function(A_list, c_vec, Sig_mat, Y_init, Y_eval, cond_idx) {
  p_v <- length(A_list); k_v <- ncol(Y_eval); H_v <- nrow(Y_eval); n_c <- length(cond_idx)
  P_c <- t(chol(Sig_mat))

  Y_unc <- matrix(0, H_v, k_v); colnames(Y_unc) <- colnames(Y_eval)
  for (h in 1:H_v) {
    pr <- c_vec
    for (j in 1:p_v) {
      lv <- if ((h - j) <= 0) Y_init[p_v + (h - j), ] else Y_unc[h - j, ]
      pr <- pr + A_list[[j]] %*% lv
    }
    Y_unc[h, ] <- pr
  }
  M_mats <- vector("list", H_v); M_mats[[1]] <- diag(k_v)
  for (m in 2:H_v) {
    Mm <- matrix(0, k_v, k_v)
    for (j in 1:min(p_v, m - 1)) Mm <- Mm + A_list[[j]] %*% M_mats[[m - j]]
    M_mats[[m]] <- Mm
  }
  R <- matrix(0, H_v * n_c, H_v * k_v); r_v <- numeric(H_v * n_c); rc <- 1
  for (h in 1:H_v) {
    for (ci in 1:n_c) {
      cv <- cond_idx[ci]; r_v[rc] <- Y_eval[h, cv] - Y_unc[h, cv]
      for (s in 1:h) R[rc, ((s - 1) * k_v + 1):(s * k_v)] <- (M_mats[[h - s + 1]] %*% P_c)[cv, ]
      rc <- rc + 1
    }
  }
  eta_all <- as.numeric(t(R) %*% solve(R %*% t(R), r_v))
  Y_cond <- Y_unc
  for (h in 1:H_v) {
    adj <- rep(0, k_v)
    for (s in 1:h) adj <- adj + as.numeric(M_mats[[h - s + 1]] %*% P_c %*% eta_all[((s - 1) * k_v + 1):(s * k_v)])
    Y_cond[h, ] <- Y_unc[h, ] + adj
  }
  colnames(Y_cond) <- colnames(Y_eval)
  list(Y_unc = Y_unc, Y_cond = Y_cond)
}

# (a) COVID-Controlled VAR(2) [Differences]
A_var <- Acoef(var_baseline)
c_var <- sapply(var_baseline$varresult, function(eq) coef(eq)["const"])
wz_var <- compute_wz_generic(A_var, c_var, Sigma_res, Y_train[(T_tr - p_opt + 1):T_tr, ], Y_test, cond_idx)

# (b) Formal Cointegrated VECM(r=2) [Levels]
A_vecm   <- vecm_r2$A
p_vecm   <- length(A_vecm)
c_vecm   <- as.numeric(vecm_r2$deterministic[, 1])
E_vecm   <- residuals(vecm_r2)
Sig_vecm <- crossprod(E_vecm) / (nrow(E_vecm) - k * p_vecm - 1 - n_dum)
wz_vecm  <- compute_wz_generic(A_vecm, c_vecm, Sig_vecm, Y_train_lvl[(T_tr - p_vecm + 1):T_tr, ], Y_test_lvl, cond_idx)

cat("======================================================================\n")
cat("6. WAGGONER-ZHA CONDITIONAL ACCURACY CHECK (VAR(2) & VECM(r=2))\n")
cat("======================================================================\n")
cat(sprintf("Max abs conditioning error VAR(2) : %.2e\n", max(abs(wz_var$Y_cond[, cond_idx] - Y_test[, cond_idx]))))
cat(sprintf("Max abs conditioning error VECM(2): %.2e\n", max(abs(wz_vecm$Y_cond[, cond_idx] - Y_test_lvl[, cond_idx]))))

# ------------------------------------------------------------------------------
# 7. Back-Transformation to Levels & Newey-West HAC Break Tests
#    Explicit level mappings (Prof. Ceylan Point 2):
#      HWI_{T+h}    = HWI_T + sum_{s=1}^h dHWI_{T+s}
#      PAYEMS_{T+h} = PAYEMS_T * exp((1/100) * sum_{s=1}^h dlogPAYEMS_{T+s})
# ------------------------------------------------------------------------------
dates_test <- test_df$date
base_hwi   <- train_df$HWI[T_tr]
base_pay   <- train_df$PAYEMS[T_tr]

lvl_real_hwi <- test_df$HWI
lvl_unc_hwi  <- base_hwi + cumsum(wz_var$Y_unc[, "HWI"])
lvl_cnd_hwi  <- base_hwi + cumsum(wz_var$Y_cond[, "HWI"])
lvl_vecm_hwi <- exp(wz_vecm$Y_cond[, "HWI"] / 100)
lvl_vecm_hwi_u <- exp(wz_vecm$Y_unc[, "HWI"] / 100)

lvl_real_pay <- test_df$PAYEMS
lvl_unc_pay  <- base_pay * exp(cumsum(wz_var$Y_unc[, "PAYEMS"] / 100))
lvl_cnd_pay  <- base_pay * exp(cumsum(wz_var$Y_cond[, "PAYEMS"] / 100))
lvl_vecm_pay <- exp(wz_vecm$Y_cond[, "PAYEMS"] / 100)
lvl_vecm_pay_u <- exp(wz_vecm$Y_unc[, "PAYEMS"] / 100)

err_unc_hwi  <- Y_test[, "HWI"]    - wz_var$Y_unc[, "HWI"]
err_cnd_hwi  <- Y_test[, "HWI"]    - wz_var$Y_cond[, "HWI"]
err_vecm_hwi <- diff(c(base_hwi, lvl_real_hwi)) - diff(c(base_hwi, lvl_vecm_hwi))

err_unc_pay  <- Y_test[, "PAYEMS"] - wz_var$Y_unc[, "PAYEMS"]
err_cnd_pay  <- Y_test[, "PAYEMS"] - wz_var$Y_cond[, "PAYEMS"]
err_vecm_pay <- diff(c(100 * log(base_pay), Y_test_lvl[, "PAYEMS"])) - diff(c(100 * log(base_pay), wz_vecm$Y_cond[, "PAYEMS"]))

hac_test <- function(err_vec, ref_unc, name, end_gap_pct) {
  mod <- lm(err_vec ~ 1)
  bw <- floor(4 * (length(err_vec) / 100)^(2/9)) # b = 3
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
    End_PctGap = round(end_gap_pct, 1),
    stringsAsFactors = FALSE
  )
}

pct_g <- function(r, c) 100 * (r / c - 1)

hac_table <- rbind(
  hac_test(err_unc_hwi,  err_unc_hwi, "Job Openings (HWI) [Unconditional VAR(2)]",     pct_g(lvl_real_hwi[H], lvl_unc_hwi[H])),
  hac_test(err_cnd_hwi,  err_unc_hwi, "Job Openings (HWI) [Policy-Cond. VAR(2)]",      pct_g(lvl_real_hwi[H], lvl_cnd_hwi[H])),
  hac_test(err_vecm_hwi, err_unc_hwi, "Job Openings (HWI) [Policy-Cond. VECM(r=2)]",   pct_g(lvl_real_hwi[H], lvl_vecm_hwi[H])),
  hac_test(err_unc_pay,  err_unc_pay, "Nonfarm Payrolls (PAYEMS) [Unconditional VAR(2)]",   pct_g(lvl_real_pay[H], lvl_unc_pay[H])),
  hac_test(err_cnd_pay,  err_unc_pay, "Nonfarm Payrolls (PAYEMS) [Policy-Cond. VAR(2)]",    pct_g(lvl_real_pay[H], lvl_cnd_pay[H])),
  hac_test(err_vecm_pay, err_unc_pay, "Nonfarm Payrolls (PAYEMS) [Policy-Cond. VECM(r=2)]", pct_g(lvl_real_pay[H], lvl_vecm_pay[H]))
)

cat("\n======================================================================\n")
cat("7. NEWEY-WEST HAC STRUCTURAL BREAK TESTS (VAR(2) & VECM(r=2))\n")
cat("======================================================================\n")
print(hac_table, row.names = FALSE)
write.csv(hac_table, "outputs/econometrics/table_hac_break_tests.csv", row.names = FALSE)

fmt_p <- function(p) {
  if (p < 0.001) { ex <- floor(log10(p)); sprintf("%.2f \\times 10^{%d}", p / 10^ex, ex) }
  else sprintf("%.3f", p)
}

tex_tab4 <- sprintf(
"\\begin{table}[t]
\\centering
\\caption{Newey--West HAC Tests: COVID-Controlled VAR(2) \\& VECM($r=2$)}
\\label{tab:hac_tests}
\\begin{minipage}{\\columnwidth}
\\centering
\\fontsize{6.6pt}{8.0pt}\\selectfont
\\setlength{\\tabcolsep}{1.5pt}
\\begin{tabular}{lccccc}
\\toprule
\\textbf{Target \\& Model Regime} & $\\bar{e}$ & \\textbf{HAC SE} & $t$-\\textbf{ratio} & $p$-\\textbf{value} & \\textbf{2026 Gap} \\\\
\\midrule
\\multicolumn{6}{l}{\\textit{JOLTS Job Openings (\\texttt{HWI}, Thousands of Openings/Month)}} \\\\
Unconditional VAR(2) & $%.2f$ & $%.2f$ & $%.2f$ & $%s$ & $%+.1f\\%%$ \\\\
Policy-Cond. VAR(2)  & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%s}$ & $\\mathbf{%+.1f\\%%}$ \\\\
Policy-Cond. VECM($r=2$) & $%.2f$ & $%.2f$ & $%.2f$ & $%s$ & $%+.1f\\%%$ \\\\
\\midrule
\\multicolumn{6}{l}{\\textit{Nonfarm Payrolls (\\texttt{PAYEMS}, Monthly \\%% Growth)}} \\\\
Unconditional VAR(2) & $%+.3f$ & $%.3f$ & $%.2f$ & $%s$ & $%+.1f\\%%$ \\\\
Policy-Cond. VAR(2)  & $%+.3f$ & $%.3f$ & $%.2f$ & $%s$ & $%+.1f\\%%$ \\\\
Policy-Cond. VECM($r=2$) & $\\mathbf{%+.3f}$ & $\\mathbf{%.3f}$ & $\\mathbf{%+.2f}$ & $\\mathbf{%s}$ & $\\mathbf{%+.1f\\%%}$ \\\\
\\bottomrule
\\end{tabular}
\\vskip 2pt
{\\scriptsize \\textit{Notes}: Frequentist VAR(2) and Case~3 VECM($r=2$) include 2020:03--2020:12 COVID pulse dummies, evaluated over 2022M11--2026M8 ($H=46$) conditioned on $\\mathcal{C}=\\{\\text{FEDFUNDS},\\text{M2SL},\\text{BAA10YM}\\}$ via Waggoner--Zha (1999) projection. HAC SEs use Newey--West (1987) Bartlett kernel ($b=3$).\\par}
\\end{minipage}
\\end{table}
",
hac_table$Mean_Error[1], hac_table$HAC_SE[1], hac_table$t_ratio[1], fmt_p(hac_table$p_value[1]), hac_table$End_PctGap[1],
hac_table$Mean_Error[2], hac_table$HAC_SE[2], hac_table$t_ratio[2], fmt_p(hac_table$p_value[2]), hac_table$End_PctGap[2],
hac_table$Mean_Error[3], hac_table$HAC_SE[3], hac_table$t_ratio[3], fmt_p(hac_table$p_value[3]), hac_table$End_PctGap[3],
hac_table$Mean_Error[4], hac_table$HAC_SE[4], hac_table$t_ratio[4], fmt_p(hac_table$p_value[4]), hac_table$End_PctGap[4],
hac_table$Mean_Error[5], hac_table$HAC_SE[5], hac_table$t_ratio[5], fmt_p(hac_table$p_value[5]), hac_table$End_PctGap[5],
hac_table$Mean_Error[6], hac_table$HAC_SE[6], hac_table$t_ratio[6], fmt_p(hac_table$p_value[6]), hac_table$End_PctGap[6]
)
cat(tex_tab4, file = "outputs/econometrics/table4_hac_break_tests.tex")

# ------------------------------------------------------------------------------
# 8. Publication Figure 4: Frequentist Counterfactual Levels (VAR(2) & VECM(r=2))
# ------------------------------------------------------------------------------
plot_traj_df <- rbind(
  data.frame(date = dates_test, Value = lvl_real_hwi, Type = "Realized Post-2022",        SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_cnd_hwi,  Type = "Policy-Cond. VAR(2)",       SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_vecm_hwi, Type = "Policy-Cond. VECM(r=2)",    SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_unc_hwi,  Type = "Unconditional VAR(2)",      SeriesLabel = "(a) Job Openings (HWI, thousands)"),
  data.frame(date = dates_test, Value = lvl_real_pay, Type = "Realized Post-2022",        SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_cnd_pay,  Type = "Policy-Cond. VAR(2)",       SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_vecm_pay, Type = "Policy-Cond. VECM(r=2)",    SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(date = dates_test, Value = lvl_unc_pay,  Type = "Unconditional VAR(2)",      SeriesLabel = "(b) Payroll Employment (PAYEMS, thousands)")
)
plot_traj_df$Type <- factor(plot_traj_df$Type, levels = c(
  "Realized Post-2022", "Policy-Cond. VECM(r=2)", "Policy-Cond. VAR(2)", "Unconditional VAR(2)"
))

build_fig4 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  ggplot(plot_traj_df, aes(x = date, y = Value, color = Type, linetype = Type)) +
    geom_line(linewidth = if (is_web) 0.90 else 0.68) +
    facet_wrap(~SeriesLabel, scales = "free_y", ncol = if (is_web) 2 else 1) +
    scale_color_manual(values = c(
      "Realized Post-2022"     = PAL$red,
      "Policy-Cond. VECM(r=2)" = PAL$teal,
      "Policy-Cond. VAR(2)"    = PAL$navy,
      "Unconditional VAR(2)"   = PAL$gray
    )) +
    scale_linetype_manual(values = c(
      "Realized Post-2022"     = "solid",
      "Policy-Cond. VECM(r=2)" = "solid",
      "Policy-Cond. VAR(2)"    = "dashed",
      "Unconditional VAR(2)"   = "dotted"
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
cat("\nStep 2 completed successfully (2020M3-2020M12 COVID dummies + VECM(r=2)).\n")
