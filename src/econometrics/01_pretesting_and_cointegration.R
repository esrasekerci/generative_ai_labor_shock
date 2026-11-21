#!/usr/bin/env Rscript
# ==============================================================================
# 01_pretesting_and_cointegration.R
# Step 1: Data Audit, Confirmatory Unit-Root Battery (ADF, PP, ERS DF-GLS,
#         KPSS, Zivot-Andrews), Multivariate Outlier Diagnostics (Wei 2019,
#         Sec. 2.7), and Johansen Cointegration Rank Tests (Enders 2015, Ch. 6)
#
# Methodological references:
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.), Ch. 4 & 6.
#   - Wei, W. W. S. (2019). Multivariate Time Series Analysis and Applications,
#     Ch. 2 (Secs. 2.5 & 2.7).
#   - Cryer, J. D. & Chan, K.-S. (2008). Time Series Analysis: With Applications
#     in R (2nd ed.), Ch. 8 & 11.
#   - Shumway, R. H. & Stoffer, D. S. (2011). Time Series Analysis and Its
#     Applications: With R Examples (3rd ed.), Sec. 5.3.
# ==============================================================================

suppressPackageStartupMessages({
  library(urca)
  library(vars)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("src/econometrics/theme_publication.R")

dir.create("outputs/econometrics", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load FRED-MD Data & Construct Exact Transformations (Zero Start-Loss)
# ------------------------------------------------------------------------------
raw <- read.csv("fred_data.csv", stringsAsFactors = FALSE)
tcodes_row <- raw[1, ]
df <- raw[-1, ]
df$date <- as.Date(df$sasdate, format = "%m/%d/%Y")
df <- df[!is.na(df$date), ]
df <- df[order(df$date), ]

num_cols <- c("FEDFUNDS", "M2SL", "BAA", "GS10", "INDPRO", "HWI", "PAYEMS")
for (col in num_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}
df$BAA10YM <- df$BAA - df$GS10

# Compute transformations on the FULL history (from 1959) before slicing >= 2010
# so zero observations are lost at 2010-01-01.
d1 <- function(x) c(NA, diff(x))
dlog1 <- function(x) c(NA, diff(log(x))) * 100       # monthly % log-change
d2log1 <- function(x) c(NA, NA, diff(diff(log(x)))) * 100  # monthly % acceleration

df$log_M2SL   <- log(df$M2SL)
df$log_INDPRO <- log(df$INDPRO)
df$log_HWI    <- log(df$HWI)
df$log_PAYEMS <- log(df$PAYEMS)

# Official FRED-MD stationary transforms (scaled cleanly for interpretability)
df$t_FEDFUNDS <- d1(df$FEDFUNDS)        # tcode 2: Δ FEDFUNDS (percentage pts)
df$t_M2SL     <- d2log1(df$M2SL)        # tcode 6: Δ² log(M2SL) * 100
df$t_M2SL_d1  <- dlog1(df$M2SL)         # diagnostic I(1) check: Δ log(M2SL) * 100
df$t_BAA10YM  <- df$BAA10YM             # tcode 1: BAA - GS10 level (percentage pts)
df$t_INDPRO   <- dlog1(df$INDPRO)       # tcode 5: Δ log(INDPRO) * 100
df$t_HWI      <- d1(df$HWI)             # tcode 2: Δ HWI (thousands of openings)
df$t_HWI_dlog <- dlog1(df$HWI)          # variance-stabilized check: Δ log(HWI) * 100
df$t_PAYEMS   <- dlog1(df$PAYEMS)       # tcode 5: Δ log(PAYEMS) * 100

# Sample splits
full_all  <- df[df$date >= as.Date("2010-01-01"), ]
full_comp <- full_all[!is.na(full_all$HWI), ] # 2010-01-01 to 2026-08-01 (T = 200)
train     <- full_comp[full_comp$date <= as.Date("2022-10-01"), ] # T = 154
test      <- full_comp[full_comp$date >= as.Date("2022-11-01"), ] # T = 46
precovid  <- train[train$date <= as.Date("2020-02-01"), ]         # T = 122

cat("======================================================================\n")
cat("1. DATA AUDIT SUMMARY\n")
cat("======================================================================\n")
cat(sprintf("Raw >= 2010-01-01 rows    : %d (%s to %s)\n",
            nrow(full_all), min(full_all$date), max(full_all$date)))
cat(sprintf("Complete multivariate rows: %d (%s to %s; HWI NA at 2026-09-01)\n",
            nrow(full_comp), min(full_comp$date), max(full_comp$date)))
cat(sprintf("Pre-AI training sample    : %d (%s to %s)\n",
            nrow(train), min(train$date), max(train$date)))
cat(sprintf("Post-AI test sample       : %d (%s to %s)\n",
            nrow(test), min(test$date), max(test$date)))
cat(sprintf("Pre-COVID subsample       : %d (%s to %s)\n",
            nrow(precovid), min(precovid$date), max(precovid$date)))

# Verify zero NAs in train and full_comp across all transformed variables
check_vars <- c("t_FEDFUNDS", "t_M2SL", "t_M2SL_d1", "t_BAA10YM",
                "t_INDPRO", "t_HWI", "t_HWI_dlog", "t_PAYEMS")
cat("Missing values in train    :", sum(is.na(train[, check_vars])), "\n")
cat("Missing values in full_comp:", sum(is.na(full_comp[, check_vars])), "\n\n")

# ------------------------------------------------------------------------------
# 2. Confirmatory Unit-Root & Stationarity Battery (Enders Ch. 4, Secs. 4.7-4.10)
# ------------------------------------------------------------------------------
# Tests:
#   (1) ADF (ur.df, BIC-selected lags up to Schwert max floor(12*(T/100)^0.25))
#   (2) Phillips-Perron (ur.pp, Z-tau, short truncation)
#   (3) ERS DF-GLS (ur.ers, higher power per Enders Sec. 4.10)
#   (4) KPSS (ur.kpss, reversed H0: stationarity)
#   (5) Zivot-Andrews (ur.za, endogenous break per Enders Sec. 4.8)

sig_stars <- function(stat, cv_vec, lower_tail = TRUE) {
  # cv_vec named or ordered 1pct, 5pct, 10pct
  if (is.na(stat)) return("")
  if (lower_tail) {
    if (stat <= cv_vec[1]) return("***")
    if (stat <= cv_vec[2]) return("**")
    if (stat <= cv_vec[3]) return("*")
    return("")
  } else {
    if (stat >= cv_vec[1]) return("***")
    if (stat >= cv_vec[2]) return("**")
    if (stat >= cv_vec[3]) return("*")
    return("")
  }
}

run_ur_battery <- function(x, dates, name, spec_type = c("trend", "drift")) {
  spec_type <- match.arg(spec_type)
  T_n <- length(x)
  max_lag <- floor(12 * (T_n / 100)^0.25)

  # (1) ADF (BIC)
  adf_obj <- ur.df(x, type = spec_type, lags = max_lag, selectlags = "BIC")
  adf_stat <- as.numeric(adf_obj@teststat[1, 1])
  adf_cv   <- as.numeric(adf_obj@cval[1, c("1pct", "5pct", "10pct")])
  adf_lag  <- as.integer(adf_obj@lags)
  # Joint deterministic statistic (phi3 if trend, phi1 if drift)
  phi_stat <- if (spec_type == "trend") as.numeric(adf_obj@teststat[1, 3]) else as.numeric(adf_obj@teststat[1, 2])
  phi_cv5  <- if (spec_type == "trend") as.numeric(adf_obj@cval[3, "5pct"]) else as.numeric(adf_obj@cval[2, "5pct"])

  # (2) Phillips-Perron
  pp_model <- ifelse(spec_type == "trend", "trend", "constant")
  pp_obj <- ur.pp(x, type = "Z-tau", model = pp_model, lags = "short")
  pp_stat <- as.numeric(pp_obj@teststat)
  pp_cv   <- as.numeric(pp_obj@cval[1, c("1pct", "5pct", "10pct")])

  # (3) ERS DF-GLS
  ers_model <- ifelse(spec_type == "trend", "trend", "constant")
  ers_obj <- ur.ers(x, type = "DF-GLS", model = ers_model, lag.max = min(6, max_lag))
  ers_stat <- as.numeric(ers_obj@teststat)
  ers_cv   <- as.numeric(ers_obj@cval[1, c("1pct", "5pct", "10pct")])

  # (4) KPSS (reversed null: H0 = stationary; upper-tail rejection)
  kpss_type <- ifelse(spec_type == "trend", "tau", "mu")
  kpss_obj <- ur.kpss(x, type = kpss_type, lags = "short")
  kpss_stat <- as.numeric(kpss_obj@teststat)
  kpss_cv   <- as.numeric(kpss_obj@cval[1, c("1pct", "5pct", "10pct")])

  # (5) Zivot-Andrews endogenous break test (model = "both" for trend, "intercept" for drift)
  za_model <- ifelse(spec_type == "trend", "both", "intercept")
  za_obj <- ur.za(x, model = za_model, lag = min(4, max_lag))
  za_stat <- as.numeric(za_obj@teststat)
  za_cv   <- as.numeric(za_obj@cval[c(1, 2, 3)]) # 1%, 5%, 10%
  za_bp   <- as.character(dates[za_obj@bpoint])

  # Decision summary
  ur_rej_count <- sum(c(adf_stat <= adf_cv[2], pp_stat <= pp_cv[2], ers_stat <= ers_cv[2]))
  kpss_rej <- kpss_stat >= kpss_cv[2]
  verdict <- if (ur_rej_count >= 2 && !kpss_rej) {
    "I(0) Stationary"
  } else if (ur_rej_count == 0 && kpss_rej) {
    "I(1) Nonstationary"
  } else if (ur_rej_count >= 2 && kpss_rej) {
    "Borderline / Break-Affected"
  } else if (ur_rej_count == 1 && !kpss_rej) {
    "Weakly Stationary / Near-UR"
  } else {
    "Nonstationary (Weak Evidence)"
  }

  data.frame(
    Variable = name,
    DetSpec = spec_type,
    ADF_lag = adf_lag,
    ADF_stat = round(adf_stat, 3),
    ADF_sig  = sig_stars(adf_stat, adf_cv, TRUE),
    Phi_stat = round(phi_stat, 2),
    Phi_sig5 = ifelse(phi_stat >= phi_cv5, "Y", "N"),
    PP_stat  = round(pp_stat, 3),
    PP_sig   = sig_stars(pp_stat, pp_cv, TRUE),
    ERS_stat = round(ers_stat, 3),
    ERS_sig  = sig_stars(ers_stat, ers_cv, TRUE),
    KPSS_stat = round(kpss_stat, 3),
    KPSS_sig  = sig_stars(kpss_stat, kpss_cv, FALSE),
    ZA_stat   = round(za_stat, 3),
    ZA_sig    = sig_stars(za_stat, za_cv, TRUE),
    ZA_break  = za_bp,
    Verdict   = verdict,
    stringsAsFactors = FALSE
  )
}

# Pre-AI training sample (2010-01 to 2022-10, T = 154)
level_specs <- list(
  list(col = "FEDFUNDS",   name = "FEDFUNDS (level)",   spec = "drift"),
  list(col = "log_M2SL",   name = "log(M2SL)",          spec = "trend"),
  list(col = "BAA10YM",    name = "BAA10YM (level)",    spec = "drift"),
  list(col = "log_INDPRO", name = "log(INDPRO)",        spec = "trend"),
  list(col = "HWI",        name = "HWI (level)",        spec = "trend"),
  list(col = "log_HWI",    name = "log(HWI)",           spec = "trend"),
  list(col = "log_PAYEMS", name = "log(PAYEMS)",        spec = "trend")
)

trans_specs <- list(
  list(col = "t_FEDFUNDS", name = "d(FEDFUNDS) [tc=2]",     spec = "drift"),
  list(col = "t_M2SL_d1",  name = "dlog(M2SL) [I(1) chk]",  spec = "drift"),
  list(col = "t_M2SL",     name = "d2log(M2SL) [tc=6]",     spec = "drift"),
  list(col = "t_BAA10YM",  name = "BAA10YM [tc=1]",         spec = "drift"),
  list(col = "t_INDPRO",   name = "dlog(INDPRO) [tc=5]",    spec = "drift"),
  list(col = "t_HWI",      name = "d(HWI) [tc=2]",          spec = "drift"),
  list(col = "t_HWI_dlog", name = "dlog(HWI) [alt]",        spec = "drift"),
  list(col = "t_PAYEMS",   name = "dlog(PAYEMS) [tc=5]",    spec = "drift")
)

ur_levels_train <- do.call(rbind, lapply(level_specs, function(s) {
  run_ur_battery(train[[s$col]], train$date, s$name, s$spec)
}))
ur_trans_train <- do.call(rbind, lapply(trans_specs, function(s) {
  run_ur_battery(train[[s$col]], train$date, s$name, s$spec)
}))

# Also check full complete sample (2010-01 to 2026-08, T = 200) for transparency
ur_trans_full <- do.call(rbind, lapply(trans_specs, function(s) {
  run_ur_battery(full_comp[[s$col]], full_comp$date, s$name, s$spec)
}))

cat("======================================================================\n")
cat("2A. UNIT-ROOT BATTERY: LEVELS (PRE-AI TRAINING SAMPLE, 2010M1-2022M10)\n")
cat("======================================================================\n")
print(ur_levels_train, row.names = FALSE)

cat("\n======================================================================\n")
cat("2B. UNIT-ROOT BATTERY: STATIONARY TRANSFORMS (TRAINING, 2010M1-2022M10)\n")
cat("======================================================================\n")
print(ur_trans_train, row.names = FALSE)

cat("\n======================================================================\n")
cat("2C. UNIT-ROOT BATTERY: STATIONARY TRANSFORMS (FULL, 2010M1-2026M8)\n")
cat("======================================================================\n")
print(ur_trans_full, row.names = FALSE)

write.csv(ur_levels_train, "outputs/econometrics/table_ur_levels_train.csv", row.names = FALSE)
write.csv(ur_trans_train,  "outputs/econometrics/table_ur_trans_train.csv",  row.names = FALSE)
write.csv(ur_trans_full,   "outputs/econometrics/table_ur_trans_full.csv",   row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Multivariate Outlier Audit (Wei 2019, Sec. 2.7, Eqs. 2.84-2.87)
# ------------------------------------------------------------------------------
# Using stationary FRED-MD system on pre-AI sample (2010-01 to 2022-10)
# Note: also check dlog(M2SL) vs d2log(M2SL) and d(HWI) vs dlog(HWI)
Y_train_std <- train[, c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")]
colnames(Y_train_std) <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")

# Fit preliminary VAR(2) to obtain innovations e_t = Pi(B) Z_t - c
p_prelim <- 2
var_prelim <- VAR(Y_train_std, p = p_prelim, type = "const")
E_mat <- residuals(var_prelim)             # (T - p) x k
Sigma_hat <- crossprod(E_mat) / nrow(E_mat)
Sigma_inv <- solve(Sigma_hat)
dates_e <- train$date[(p_prelim + 1):nrow(train)]
T_e <- nrow(E_mat)
k_dim <- ncol(E_mat)

# (a) Multivariate Innovational Outlier (MIO) test statistic J_{I,h} (Wei Eq. 2.87)
#     For MIO, alpha(B) = Psi(B) = Pi(B)^{-1}, so Pi(B) alpha(B) = I, Gamma_0 = -I, Gamma_j = 0 for j>=1
#     => omega_{I,h} = e_h, Sigma_{I,h} = Sigma, J_{I,h} = e_h' Sigma^{-1} e_h ~ chi^2_k!
J_IO <- apply(E_mat, 1, function(e) as.numeric(t(e) %*% Sigma_inv %*% e))
p_IO <- 1 - pchisq(J_IO, df = k_dim)

# (b) Multivariate Additive Outlier (MAO) test statistic J_{A,h} (Wei Eq. 2.86-2.87)
#     For VAR(p), Pi(B) = I - A_1 B - ... - A_p B^p, so alpha(B) = I => Gamma_0 = -I, Gamma_j = A_j (j=1..p)
A_list <- Acoef(var_prelim) # list of p matrices (k x k)
Gamma_list <- c(list(-diag(k_dim)), A_list) # index 1 -> j=0, index j+1 -> j

J_AO <- numeric(T_e)
for (h in 1:T_e) {
  max_j <- min(p_prelim, T_e - h)
  S_inv_sum <- matrix(0, k_dim, k_dim)
  rhs_sum   <- rep(0, k_dim)
  for (j in 0:max_j) {
    Gj <- Gamma_list[[j + 1]]
    S_inv_sum <- S_inv_sum + t(Gj) %*% Sigma_inv %*% Gj
    rhs_sum   <- rhs_sum + as.numeric(t(Gj) %*% Sigma_inv %*% E_mat[h + j, ])
  }
  Sigma_Ah <- solve(S_inv_sum)
  omega_Ah <- -as.numeric(Sigma_Ah %*% rhs_sum)
  J_AO[h]  <- as.numeric(t(omega_Ah) %*% S_inv_sum %*% omega_Ah)
}
p_AO <- 1 - pchisq(J_AO, df = k_dim)

# Also compute per-variable standardized residual z-scores
sd_e <- sqrt(diag(Sigma_hat))
Z_e <- sweep(E_mat, 2, sd_e, "/")
max_z_var <- apply(abs(Z_e), 1, function(z) colnames(Z_e)[which.max(z)])
max_z_val <- apply(abs(Z_e), 1, max)

outlier_df <- data.frame(
  date = dates_e,
  J_MIO = round(J_IO, 2),
  p_MIO = signif(p_IO, 3),
  J_MAO = round(J_AO, 2),
  p_MAO = signif(p_AO, 3),
  max_abs_z = round(max_z_val, 2),
  dominant_var = max_z_var,
  z_FEDFUNDS = round(Z_e[, "FEDFUNDS"], 2),
  z_M2SL     = round(Z_e[, "M2SL"], 2),
  z_BAA10YM  = round(Z_e[, "BAA10YM"], 2),
  z_INDPRO   = round(Z_e[, "INDPRO"], 2),
  z_HWI      = round(Z_e[, "HWI"], 2),
  z_PAYEMS   = round(Z_e[, "PAYEMS"], 2),
  stringsAsFactors = FALSE
)

# Bonferroni-adjusted critical value for T_e tests at alpha = 0.01 (chi^2_6)
bonf_cv_01 <- qchisq(1 - 0.01 / T_e, df = k_dim)
bonf_cv_05 <- qchisq(1 - 0.05 / T_e, df = k_dim)

cat("\n======================================================================\n")
cat("3. MULTIVARIATE OUTLIER AUDIT (WEI 2019, SEC. 2.7: MIO & MAO TESTS)\n")
cat("======================================================================\n")
cat(sprintf("Bonferroni chi^2(6) critical values: 5%% = %.2f, 1%% = %.2f (T = %d)\n",
            bonf_cv_05, bonf_cv_01, T_e))
top_outliers <- outlier_df[order(-outlier_df$J_MIO), ][1:10, ]
print(top_outliers, row.names = FALSE)
write.csv(outlier_df, "outputs/econometrics/table_multivariate_outliers.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Johansen Cointegration Analysis (Enders Ch. 6, Secs. 6.7-6.9, pp. 384-387)
# ------------------------------------------------------------------------------
# Following Enders (p. 385):
#   1) Select lag length p in levels starting from p_max = 12 (T^(1/3) -> monthly 12)
#   2) Test both Trace and Max-Eigenvalue on the pure I(1) variables:
#      (i) Core 5 I(1) system: FEDFUNDS, log(M2SL), log(INDPRO), log(HWI), log(PAYEMS)
#      (ii) Same 5 I(1) system with COVID-19 pulse dummies (2020M3-2020M5)
#      (iii) Pre-COVID subsample (2010M1-2020M2) to verify COVID distortion
#      (iv) Full 6-variable system (adding BAA10YM) as sensitivity check

Y_lvl_5 <- train[, c("FEDFUNDS", "log_M2SL", "log_INDPRO", "log_HWI", "log_PAYEMS")]
Y_lvl_6 <- train[, c("FEDFUNDS", "log_M2SL", "BAA10YM", "log_INDPRO", "log_HWI", "log_PAYEMS")]
Y_lvl_5_pre <- precovid[, c("FEDFUNDS", "log_M2SL", "log_INDPRO", "log_HWI", "log_PAYEMS")]

# Full 2020:03-2020:12 COVID pulse dummies (10 months, per Prof. Ceylan feedback)
covid_months <- seq(as.Date("2020-03-01"), as.Date("2020-12-01"), by = "month")
covid_dum <- sapply(covid_months, function(d) as.numeric(train$date == d))
colnames(covid_dum) <- format(covid_months, "D%YM%m")

cat("\n======================================================================\n")
cat("4A. LEVEL-VAR LAG SELECTION FOR JOHANSEN TEST (p_max = 12)\n")
cat("======================================================================\n")
sel_lvl5     <- VARselect(Y_lvl_5, lag.max = 12, type = "const")
sel_lvl5_cov <- VARselect(Y_lvl_5, lag.max = 12, type = "const", exogen = covid_dum)
sel_lvl5_pre <- VARselect(Y_lvl_5_pre, lag.max = 12, type = "const")
cat("Core 5 I(1) Raw (2010M1-2022M10):\n"); print(sel_lvl5$selection)
cat("Core 5 I(1) COVID-Controlled:\n");     print(sel_lvl5_cov$selection)
cat("Core 5 I(1) Pre-COVID (2010M1-2020M2):\n"); print(sel_lvl5_pre$selection)

summarize_johansen <- function(Y_mat, K_val, label, dum = NULL, ecdet = "const") {
  jt <- ca.jo(Y_mat, type = "trace", ecdet = ecdet, K = K_val, spec = "transitory", dumvar = dum)
  je <- ca.jo(Y_mat, type = "eigen", ecdet = ecdet, K = K_val, spec = "transitory", dumvar = dum)
  k <- ncol(Y_mat)
  # urca orders rows from r <= k-1 down to r = 0; reverse to r = 0 .. k-1
  idx <- rev(seq_len(k))
  res <- data.frame(
    Specification = label,
    K_lags = K_val,
    Null_H0 = paste0("r <= ", 0:(k - 1)),
    Eigenvalue = round(jt@lambda[1:k], 4),
    Trace_Stat = round(jt@teststat[idx], 2),
    Trace_5pct = round(jt@cval[idx, "5pct"], 2),
    Trace_1pct = round(jt@cval[idx, "1pct"], 2),
    Trace_Rej5 = ifelse(jt@teststat[idx] > jt@cval[idx, "5pct"], "Reject", "Fail"),
    MaxEig_Stat = round(je@teststat[idx], 2),
    MaxEig_5pct = round(je@cval[idx, "5pct"], 2),
    MaxEig_1pct = round(je@cval[idx, "1pct"], 2),
    MaxEig_Rej5 = ifelse(je@teststat[idx] > je@cval[idx, "5pct"], "Reject", "Fail"),
    stringsAsFactors = FALSE
  )
  res
}

# Use K = max(2, BIC-selected) and also test K = 3 / 4 sensitivity
K_raw <- max(2, as.integer(sel_lvl5$selection["SC(n)"]))
K_cov <- max(2, as.integer(sel_lvl5_cov$selection["SC(n)"]))
K_pre <- max(2, as.integer(sel_lvl5_pre$selection["SC(n)"]))

joh_results <- rbind(
  summarize_johansen(Y_lvl_5,     K_raw, "Core 5 I(1) Raw (2010-2022)", dum = NULL, ecdet = "const"),
  summarize_johansen(Y_lvl_5,     K_cov, "Core 5 I(1) + COVID Dummies 2020M3-12 (2010-2022)", dum = covid_dum, ecdet = "const"),
  summarize_johansen(Y_lvl_5_pre, K_pre, "Core 5 I(1) Pre-COVID (2010-2020M2)", dum = NULL, ecdet = "const"),
  summarize_johansen(Y_lvl_6,     3,     "All 6 Vars + COVID Dummies 2020M3-12 [Case 3 Drift]", dum = covid_dum, ecdet = "none")
)

cat("\n======================================================================\n")
cat("4B. JOHANSEN COINTEGRATION RANK TESTS (TRACE & MAX-EIGENVALUE)\n")
cat("======================================================================\n")
print(joh_results, row.names = FALSE)
write.csv(joh_results, "outputs/econometrics/table_johansen_cointegration.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Publication-Grade Diagnostic Figure: Multivariate Outlier Audit (COVID-19)
#    Unified punto for IEEE print (3.5 in) and GitHub.io web widescreen (8.4 in)
# ------------------------------------------------------------------------------
covid_rows <- outlier_df[outlier_df$date %in% as.Date(c("2020-03-01", "2020-04-01", "2020-05-01")), ]
z_long <- do.call(rbind, lapply(c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS"), function(v) {
  data.frame(
    Month = format(covid_rows$date, "%Y-%m"),
    Variable = v,
    Z_score = covid_rows[[paste0("z_", v)]],
    stringsAsFactors = FALSE
  )
}))
z_long$Variable <- factor(z_long$Variable, levels = c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS"))

build_fig1 <- function(is_web = FALSE) {
  txt_sz <- if (is_web) 3.0 else 2.2
  pt_sz  <- if (is_web) 1.7 else 1.0
  lw_seg <- if (is_web) 0.45 else 0.30
  th     <- if (is_web) theme_web_wide() else theme_ieee_print()
  ttl_a  <- "(a) MIO Statistic (Tsay-Pena-Pankratz)"
  ttl_b  <- "(b) Standardized Innovations (Mar-May 2020)"

  p_mio <- ggplot(outlier_df, aes(x = date, y = J_MIO)) +
    geom_segment(aes(xend = date, yend = 0), color = PAL$slate, linewidth = lw_seg, alpha = 0.85) +
    geom_point(aes(color = J_MIO > bonf_cv_01), size = pt_sz) +
    scale_color_manual(values = c("FALSE" = PAL$slate, "TRUE" = PAL$red), guide = "none") +
    geom_hline(yintercept = bonf_cv_01, linetype = "dashed", color = PAL$red, linewidth = 0.45) +
    annotate("text", x = as.Date("2010-06-01"), y = bonf_cv_01 + 10,
             label = sprintf("Bonferroni 1%% CV (%.2f)", bonf_cv_01),
             hjust = 0, size = txt_sz, color = PAL$red, fontface = "italic") +
    annotate("text", x = as.Date("2020-01-01"), y = max(outlier_df$J_MIO) * 0.92,
             label = "COVID-19\n(Mar-May 2020)", hjust = 1.08, size = txt_sz,
             color = PAL$red, fontface = "bold", lineheight = 0.92) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(title = ttl_a, x = NULL, y = expression(paste("MIO ", J[I*","*h]))) +
    th +
    theme(plot.title = element_text(face = "bold", size = if (is_web) 10.0 else 7.2, margin = margin(b = 2)))

  p_z <- ggplot(z_long, aes(x = Variable, y = Z_score, fill = Month)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    geom_hline(yintercept = c(-3, 3), linetype = "dotted", color = PAL$gray, linewidth = 0.4) +
    geom_hline(yintercept = 0, color = PAL$slate, linewidth = 0.3) +
    scale_fill_manual(values = c("2020-03" = PAL$navy, "2020-04" = PAL$red, "2020-05" = PAL$amber)) +
    labs(title = ttl_b, x = NULL, y = "Residual z-score", fill = "Month") +
    th +
    theme(
      plot.title = element_text(face = "bold", size = if (is_web) 10.0 else 7.2, margin = margin(b = 2)),
      legend.position = "bottom",
      axis.text.x = element_text(size = if (is_web) 8.8 else 6.0)
    )

  if (is_web) {
    (p_mio | p_z) + plot_layout(widths = c(1.1, 1.0))
  } else {
    (p_mio / p_z) + plot_layout(heights = c(1.05, 1.0))
  }
}

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/assets/images", recursive = TRUE, showWarnings = FALSE)

fig_out_print <- build_fig1(is_web = FALSE)
ggsave("outputs/figures/fig1_covid_outlier_audit.pdf", fig_out_print, width = 3.5, height = 3.3)
ggsave("outputs/figures/fig1_covid_outlier_audit.png", fig_out_print, width = 3.5, height = 3.3, dpi = 300)

fig_out_web <- build_fig1(is_web = TRUE)
ggsave("outputs/figures/fig1_web_outlier_audit.png", fig_out_web, width = 8.4, height = 3.9, dpi = 300)
file.copy("outputs/figures/fig1_web_outlier_audit.png", "docs/assets/images/fig1_covid_outlier_audit.png", overwrite = TRUE)

cat("\nSaved all Step 1 tables and unified outlier diagnostic figures to outputs/\n")
