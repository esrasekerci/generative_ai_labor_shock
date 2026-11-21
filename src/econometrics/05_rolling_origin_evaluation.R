#!/usr/bin/env Rscript
# ==============================================================================
# 05_rolling_origin_evaluation.R
# Step 5: Part B — Rolling-Origin Multi-Horizon Predictive Horse Race (h = 1..12)
#         Addressing Prof. Ceylan Point 4 (Required Course Forecast Evaluation):
#         - Expanding-window rolling origins across Clean Pre-COVID (2016M1-2019M2,
#           38 origins) and Post-2022 (2022M11-2025M8, 34 origins) => 72 origins
#         - Strictly zero post-2022 hyperparameter tuning (hyperparameters fixed
#           from pre-2022 marginal likelihood optimization!)
#         - Models ranked: RW-Drift, COVID-VAR(2), VECM(r=2), Unscaled BVAR(2),
#           LP-BVAR(2) [Diff], and Cointegration-Preserving LP-BVAR(3) [Levels]
#         - Exact MASE (Hyndman & Koehler 2006), analytical Gaussian CRPS
#           (Gneiting & Raftery 2007), and Diebold-Mariano (1995) tests with
#           Harvey-Leybourne-Newbold (1997) adjustment & Newey-West HAC variance.
# ==============================================================================

suppressPackageStartupMessages({
  library(vars)
  library(urca)
  library(MASS)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("src/econometrics/theme_publication.R")
set.seed(563)

dir.create("outputs/econometrics", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("docs/assets/images", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load FRED-MD Data & Prepare Systems
# ------------------------------------------------------------------------------
raw <- read.csv("fred_data.csv", stringsAsFactors = FALSE)
df <- raw[-1, ]
df$date <- as.Date(df$sasdate, format = "%m/%d/%Y")
df <- df[!is.na(df$date), ]
df <- df[order(df$date), ]

for (col in c("FEDFUNDS", "M2SL", "BAA", "GS10", "INDPRO", "HWI", "PAYEMS")) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}
df$BAA10YM <- df$BAA - df$GS10

d1     <- function(x) c(NA, diff(x))
dlog1  <- function(x) c(NA, diff(log(x))) * 100
d2log1 <- function(x) c(NA, NA, diff(diff(log(x)))) * 100

df$t_FEDFUNDS <- d1(df$FEDFUNDS)
df$t_M2SL     <- d2log1(df$M2SL)
df$t_BAA10YM  <- df$BAA10YM
df$t_INDPRO   <- dlog1(df$INDPRO)
df$t_HWI      <- d1(df$HWI)
df$t_PAYEMS   <- dlog1(df$PAYEMS)

df$L_FEDFUNDS <- df$FEDFUNDS
df$L_M2SL     <- 100 * log(df$M2SL)
df$L_BAA10YM  <- df$BAA10YM
df$L_INDPRO   <- 100 * log(df$INDPRO)
df$L_HWI      <- 100 * log(df$HWI)
df$L_PAYEMS   <- 100 * log(df$PAYEMS)

fc <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]
rownames(fc) <- NULL

var_names <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")
diff_cols <- c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")
lvl_cols  <- c("L_FEDFUNDS", "L_M2SL", "L_BAA10YM", "L_INDPRO", "L_HWI", "L_PAYEMS")
k <- 6
H_max <- 12

# Load fixed pre-2022 optimal LP-BVAR hyperparameters (ZERO post-2022 tuning!)
hp <- read.csv("outputs/econometrics/table_bvar_hyperparameters.csv", stringsAsFactors = FALSE)
lam_unsc <- hp$lambda1[1]
lam_lpd  <- hp$lambda1[2]; s0_lpd <- hp$s_2020M3[2]; s1_lpd <- hp$s_2020M4[2]; s2_lpd <- hp$s_2020M5[2]; rho_lpd <- hp$rho[2]
lam_lpl  <- hp$lambda1[3]; s0_lpl <- hp$s_2020M3[3]; s1_lpl <- hp$s_2020M4[3]; s2_lpl <- hp$s_2020M5[3]; rho_lpl <- hp$rho[3]

# Pre-COVID naive 1-step MAE denominator for MASE (Hyndman & Koehler 2006)
pre_df <- fc[fc$date <= as.Date("2020-02-01"), ]
scale_hwi <- mean(abs(diff(pre_df$HWI)))       # exact in-sample naive MAE on HWI levels
scale_pay <- mean(abs(diff(pre_df$PAYEMS)))    # exact in-sample naive MAE on PAYEMS levels

# Exact Analytical Gaussian CRPS (Gneiting & Raftery 2007, JASA)
crps_norm <- function(y, mu, sd) {
  sd <- pmax(sd, 1e-6)
  z  <- (y - mu) / sd
  sd * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
}

# Exact Diebold-Mariano (1995) with Harvey-Leybourne-Newbold (1997) adjustment
dm_test_hln <- function(loss_cand, loss_bench, h = 1) {
  d <- loss_cand - loss_bench # negative => candidate beats benchmark
  n <- length(d)
  if (n < 5 || sd(d) < 1e-12) return(list(stat = 0, pval = 1))
  mod <- lm(d ~ 1)
  bw  <- max(0, h - 1)
  v_hac <- if (bw == 0) var(d) / n else as.numeric(NeweyWest(mod, lag = bw, prewhite = FALSE))
  dm_raw <- mean(d) / sqrt(max(v_hac, 1e-16))
  # Harvey-Leybourne-Newbold (1997) finite-sample multiplier
  hln_fac <- sqrt((n + 1 - 2 * h + (h * (h - 1)) / n) / n)
  dm_hln  <- dm_raw * hln_fac
  pval    <- 2 * pt(-abs(dm_hln), df = n - 1)
  list(stat = dm_hln, pval = pval)
}

# ------------------------------------------------------------------------------
# 2. Forecast Helpers Returning (Mean, SD) in EXACT LEVELS for h = 1..12
# ------------------------------------------------------------------------------
# Helper: Forecast from VAR(p) representation (A_list, c_vec, Sig)
fcst_var_rep <- function(A_list, c_vec, Sig, Y_init, H_max, is_levels, base_hwi, base_pay) {
  p <- length(A_list); kn <- ncol(Y_init)
  Y_f <- matrix(0, H_max, kn)
  for (h in 1:H_max) {
    pr <- c_vec
    for (j in 1:p) pr <- pr + A_list[[j]] %*% (if (h - j <= 0) Y_init[p + h - j, ] else Y_f[h - j, ])
    Y_f[h, ] <- pr
  }
  # MA multipliers M_0..M_{H-1} for h-step forecast error covariance
  M_m <- vector("list", H_max); M_m[[1]] <- diag(kn)
  for (m in 2:H_max) {
    Mt <- matrix(0, kn, kn)
    for (j in 1:min(p, m - 1)) Mt <- Mt + A_list[[j]] %*% M_m[[m - j]]
    M_m[[m]] <- Mt
  }
  # Cumulative multipliers C_h = sum_{j=0}^{h-1} M_j for differenced models
  sd_hwi <- numeric(H_max); sd_pay <- numeric(H_max)
  mu_hwi <- numeric(H_max); mu_pay <- numeric(H_max)

  if (is_levels) {
    V_acc <- matrix(0, kn, kn)
    for (h in 1:H_max) {
      V_acc <- V_acc + M_m[[h]] %*% Sig %*% t(M_m[[h]])
      # Y_f[,5] = 100*log(HWI), Y_f[,6] = 100*log(PAYEMS)
      mu_lh <- Y_f[h, 5] / 100; s_lh <- sqrt(max(V_acc[5, 5], 1e-10)) / 100
      mu_lp <- Y_f[h, 6] / 100; s_lp <- sqrt(max(V_acc[6, 6], 1e-10)) / 100
      mu_hwi[h] <- exp(mu_lh)
      mu_pay[h] <- exp(mu_lp)
      # Delta-method level predictive SD: exp(mu) * s_log
      sd_hwi[h] <- mu_hwi[h] * s_lh
      sd_pay[h] <- mu_pay[h] * s_lp
    }
  } else {
    # Differenced system: cumulative forecast error sum_{s=1}^h e_{T+s} has covariance sum_{s=0}^{h-1} C_s Sig C_s'
    C_list <- vector("list", H_max); C_list[[1]] <- M_m[[1]]
    for (m in 2:H_max) C_list[[m]] <- C_list[[m - 1]] + M_m[[m]]
    V_cum <- matrix(0, kn, kn)
    cum_hwi <- cumsum(Y_f[, 5])
    cum_pay <- cumsum(Y_f[, 6]) / 100
    for (h in 1:H_max) {
      V_cum <- V_cum + C_list[[h]] %*% Sig %*% t(C_list[[h]])
      mu_hwi[h] <- base_hwi + cum_hwi[h]
      sd_hwi[h] <- sqrt(max(V_cum[5, 5], 1e-10))
      mu_pay[h] <- base_pay * exp(cum_pay[h])
      sd_pay[h] <- mu_pay[h] * (sqrt(max(V_cum[6, 6], 1e-10)) / 100)
    }
  }
  list(mu_hwi = mu_hwi, sd_hwi = sd_hwi, mu_pay = mu_pay, sd_pay = sd_pay)
}

# Helper: Fast fixed-hyperparameter BVAR estimator for rolling origin
fit_bvar_fixed <- function(Y_tr, dates_tr, p, lam1, rw_prior, covid_scale, s0, s1, s2, rho) {
  Tn <- nrow(Y_tr); Kr <- 1 + k * p
  Ye <- Y_tr[(p + 1):Tn, ]; Te <- nrow(Ye)
  Xe <- matrix(0, Te, Kr); Xe[, 1] <- 1
  for (j in 1:p) Xe[, (1 + (j - 1) * k + 1):(1 + j * k)] <- Y_tr[(p + 1 - j):(Tn - j), ]
  dates_e <- dates_tr[(p + 1):Tn]

  has_cov <- any(dates_e == as.Date("2020-03-01"))
  idx_sig <- if (has_cov) which(dates_tr <= as.Date("2020-02-01")) else seq_len(Tn)
  Y_sig <- Y_tr[idx_sig, ]; Tp <- nrow(Y_sig)
  sig_p <- numeric(k)
  for (i in 1:k) {
    yi <- Y_sig[, i]; Yi <- yi[(p + 1):Tp]
    Xi <- cbind(1, sapply(1:p, function(l) yi[(p + 1 - l):(Tp - l)]))
    ri <- Yi - Xi %*% solve(crossprod(Xi), crossprod(Xi, Yi))
    sig_p[i] <- sqrt(sum(ri^2) / (length(Yi) - ncol(Xi)))
  }

  s_vec <- rep(1, Te)
  if (covid_scale && has_cov) {
    i0 <- which(dates_e == as.Date("2020-03-01"))
    s_vec[i0] <- s0
    if (i0 + 1 <= Te) s_vec[i0 + 1] <- s1
    if (i0 + 2 <= Te) for (t in (i0 + 2):Te) s_vec[t] <- 1 + (s2 - 1) * (rho^(t - (i0 + 2)))
  }

  d0 <- k + 2; Psi0 <- diag(sig_p^2); B0 <- matrix(0, Kr, k)
  if (rw_prior) for (i in 1:k) B0[1 + i, i] <- 1.0
  od <- numeric(Kr); od[1] <- 1e4
  for (l in 1:p) for (i in 1:k) od[1 + (l - 1) * k + i] <- (lam1^2) / ((l^2) * (sig_p[i]^2))
  Om_inv <- diag(1 / od)

  Yt <- Ye / s_vec; Xt <- Xe / s_vec
  Bp <- solve(crossprod(Xt) + Om_inv, crossprod(Xt, Yt) + Om_inv %*% B0)
  Psip <- Psi0 + crossprod(Yt - Xt %*% Bp) + t(Bp - B0) %*% Om_inv %*% (Bp - B0)
  Sig <- Psip / (d0 + Te - k - 1)
  A_l <- lapply(1:p, function(l) t(Bp[(1 + (l - 1) * k + 1):(1 + l * k), ]))
  list(A = A_l, c = Bp[1, ], Sig = Sig)
}

# ------------------------------------------------------------------------------
# 3. Execute Expanding-Window Rolling-Origin Evaluation (72 Origins x 12 Horizons)
#    Clean Pre-COVID Origins: 2016M1-2019M2 (38 origins, targets <= 2020M2)
#    Post-2022 Origins      : 2022M10-2025M8 (35 origins, targets 2022M11-2026M8)
# ------------------------------------------------------------------------------
orig_pre  <- which(fc$date >= as.Date("2016-01-01") & fc$date <= as.Date("2019-02-01"))
orig_post <- which(fc$date >= as.Date("2022-10-01") & fc$date <= as.Date("2025-08-01"))
origins   <- c(orig_pre, orig_post)
N_orig    <- length(origins)

models <- c("RW-Drift", "COVID-VAR(2)", "VECM(r=2)", "Unscaled-BVAR(2)", "LP-BVAR(2) [Diff]", "LP-BVAR(3) [Levels]")
n_mod  <- length(models)

cat("======================================================================\n")
cat("ROLLING-ORIGIN MULTI-HORIZON FORECAST EVALUATION (h = 1..12)\n")
cat("======================================================================\n")
cat(sprintf("Total Origins: %d (%d Clean Pre-COVID [2016M1-2019M2] + %d Post-2022 [2022M10-2025M8])\n",
            N_orig, length(orig_pre), length(orig_post)))
cat(sprintf("MASE Scaling Denominators (Pre-COVID Naive MAE): HWI = %.1fk | PAYEMS = %.1fk\n\n",
            scale_hwi, scale_pay))

# Storage arrays: [origin, horizon, model]
ae_hwi   <- array(NA, dim = c(N_orig, H_max, n_mod), dimnames = list(NULL, 1:H_max, models))
crps_hwi <- array(NA, dim = c(N_orig, H_max, n_mod), dimnames = list(NULL, 1:H_max, models))
ae_pay   <- array(NA, dim = c(N_orig, H_max, n_mod), dimnames = list(NULL, 1:H_max, models))
crps_pay <- array(NA, dim = c(N_orig, H_max, n_mod), dimnames = list(NULL, 1:H_max, models))

covid_months_all <- seq(as.Date("2020-03-01"), as.Date("2020-12-01"), by = "month")

t0 <- Sys.time()
for (oi in seq_len(N_orig)) {
  t_idx <- origins[oi]
  tr_o  <- fc[1:t_idx, ]
  te_o  <- fc[(t_idx + 1):(t_idx + H_max), ]
  d_tr  <- tr_o$date

  Y_tr_d <- as.matrix(tr_o[, diff_cols]); colnames(Y_tr_d) <- var_names
  Y_tr_l <- as.matrix(tr_o[, lvl_cols]);  colnames(Y_tr_l) <- var_names
  bh <- tr_o$HWI[t_idx]; bp <- tr_o$PAYEMS[t_idx]
  y_hwi <- te_o$HWI; y_pay <- te_o$PAYEMS

  # COVID 2020:03-2020:12 dummies active if inside training window
  has_cov <- any(d_tr %in% covid_months_all)
  dum_o <- if (has_cov) {
    dm <- sapply(covid_months_all, function(d) as.numeric(d_tr == d))
    colnames(dm) <- format(covid_months_all, "D%YM%m")
    dm[, colSums(dm) > 0, drop = FALSE]
  } else NULL

  # (1) RW-Drift (Univariate Random Walk with Drift on HWI & log(PAYEMS))
  dh_tr <- diff(tr_o$HWI); dp_tr <- diff(log(tr_o$PAYEMS))
  # Exclude COVID months from RW-Drift variance/drift if post-2020 so benchmark is fair!
  ok_rw <- if (has_cov) !(d_tr[-1] %in% covid_months_all) else rep(TRUE, length(dh_tr))
  m_h <- mean(dh_tr[ok_rw]); s_h <- sd(dh_tr[ok_rw])
  m_p <- mean(dp_tr[ok_rw]); s_p <- sd(dp_tr[ok_rw])
  hs  <- 1:H_max
  rw_mu_h <- bh + hs * m_h;       rw_sd_h <- sqrt(hs) * s_h
  rw_mu_p <- bp * exp(hs * m_p);  rw_sd_p <- rw_mu_p * sqrt(hs) * s_p
  ae_hwi[oi, , 1]   <- abs(y_hwi - rw_mu_h); crps_hwi[oi, , 1] <- crps_norm(y_hwi, rw_mu_h, rw_sd_h)
  ae_pay[oi, , 1]   <- abs(y_pay - rw_mu_p); crps_pay[oi, , 1] <- crps_norm(y_pay, rw_mu_p, rw_sd_p)

  # (2) COVID-VAR(2)
  v2 <- if (is.null(dum_o)) VAR(Y_tr_d, p = 2, type = "const") else VAR(Y_tr_d, p = 2, type = "const", exogen = dum_o)
  Ev2 <- residuals(v2); Sig_v2 <- crossprod(Ev2) / (nrow(Ev2) - length(coef(v2$varresult[[1]])))
  f2 <- fcst_var_rep(Acoef(v2), sapply(v2$varresult, function(eq) coef(eq)["const"]),
                     Sig_v2, Y_tr_d[(t_idx - 1):t_idx, ], H_max, FALSE, bh, bp)
  ae_hwi[oi, , 2]   <- abs(y_hwi - f2$mu_hwi); crps_hwi[oi, , 2] <- crps_norm(y_hwi, f2$mu_hwi, f2$sd_hwi)
  ae_pay[oi, , 2]   <- abs(y_pay - f2$mu_pay); crps_pay[oi, , 2] <- crps_norm(y_pay, f2$mu_pay, f2$sd_pay)

  # (3) Formal Cointegrated VECM(r=2)
  jo_o <- ca.jo(Y_tr_l, type = "trace", ecdet = "none", K = 3, spec = "transitory", dumvar = dum_o)
  vv_o <- vec2var(jo_o, r = 2)
  Evv  <- residuals(vv_o); nd_v <- if (is.null(dum_o)) 0 else ncol(dum_o)
  Sig_vv <- crossprod(Evv) / max(10, nrow(Evv) - k * 3 - 1 - nd_v)
  f3 <- fcst_var_rep(vv_o$A, as.numeric(vv_o$deterministic[, 1]),
                     Sig_vv, Y_tr_l[(t_idx - 2):t_idx, ], H_max, TRUE, bh, bp)
  ae_hwi[oi, , 3]   <- abs(y_hwi - f3$mu_hwi); crps_hwi[oi, , 3] <- crps_norm(y_hwi, f3$mu_hwi, f3$sd_hwi)
  ae_pay[oi, , 3]   <- abs(y_pay - f3$mu_pay); crps_pay[oi, , 3] <- crps_norm(y_pay, f3$mu_pay, f3$sd_pay)

  # (4) Unscaled Minnesota BVAR(2) [Diff]
  b4 <- fit_bvar_fixed(Y_tr_d, d_tr, 2, lam_unsc, FALSE, FALSE, 1, 1, 1, 0)
  f4 <- fcst_var_rep(b4$A, b4$c, b4$Sig, Y_tr_d[(t_idx - 1):t_idx, ], H_max, FALSE, bh, bp)
  ae_hwi[oi, , 4]   <- abs(y_hwi - f4$mu_hwi); crps_hwi[oi, , 4] <- crps_norm(y_hwi, f4$mu_hwi, f4$sd_hwi)
  ae_pay[oi, , 4]   <- abs(y_pay - f4$mu_pay); crps_pay[oi, , 4] <- crps_norm(y_pay, f4$mu_pay, f4$sd_pay)

  # (5) Lenza-Primiceri BVAR(2) [Diff]
  b5 <- fit_bvar_fixed(Y_tr_d, d_tr, 2, lam_lpd, FALSE, TRUE, s0_lpd, s1_lpd, s2_lpd, rho_lpd)
  f5 <- fcst_var_rep(b5$A, b5$c, b5$Sig, Y_tr_d[(t_idx - 1):t_idx, ], H_max, FALSE, bh, bp)
  ae_hwi[oi, , 5]   <- abs(y_hwi - f5$mu_hwi); crps_hwi[oi, , 5] <- crps_norm(y_hwi, f5$mu_hwi, f5$sd_hwi)
  ae_pay[oi, , 5]   <- abs(y_pay - f5$mu_pay); crps_pay[oi, , 5] <- crps_norm(y_pay, f5$mu_pay, f5$sd_pay)

  # (6) Cointegration-Preserving Levels LP-BVAR(3)
  b6 <- fit_bvar_fixed(Y_tr_l, d_tr, 3, lam_lpl, TRUE, TRUE, s0_lpl, s1_lpl, s2_lpl, rho_lpl)
  f6 <- fcst_var_rep(b6$A, b6$c, b6$Sig, Y_tr_l[(t_idx - 2):t_idx, ], H_max, TRUE, bh, bp)
  ae_hwi[oi, , 6]   <- abs(y_hwi - f6$mu_hwi); crps_hwi[oi, , 6] <- crps_norm(y_hwi, f6$mu_hwi, f6$sd_hwi)
  ae_pay[oi, , 6]   <- abs(y_pay - f6$mu_pay); crps_pay[oi, , 6] <- crps_norm(y_pay, f6$mu_pay, f6$sd_pay)
}
cat(sprintf("Completed %d rolling origins x 6 models x 12 horizons in %.1f seconds.\n\n",
            N_orig, difftime(Sys.time(), t0, units = "secs")))

# ------------------------------------------------------------------------------
# 4. Aggregate MASE, CRPS, & Diebold-Mariano Tests (vs. COVID-VAR(2) & RW-Drift)
# ------------------------------------------------------------------------------
horiz_rows <- list()
for (m_i in seq_len(n_mod)) {
  for (h in 1:H_max) {
    dm_h_rw  <- dm_test_hln(ae_hwi[, h, m_i],   ae_hwi[, h, 1],   h = h)
    dm_h_var <- dm_test_hln(crps_hwi[, h, m_i], crps_hwi[, h, 2], h = h)
    dm_p_rw  <- dm_test_hln(ae_pay[, h, m_i],   ae_pay[, h, 1],   h = h)
    dm_p_var <- dm_test_hln(crps_pay[, h, m_i], crps_pay[, h, 2], h = h)
    horiz_rows[[length(horiz_rows) + 1]] <- data.frame(
      Model          = models[m_i],
      Horizon        = h,
      HWI_MASE       = round(mean(ae_hwi[, h, m_i]) / scale_hwi, 3),
      HWI_CRPS       = round(mean(crps_hwi[, h, m_i]), 1),
      HWI_DM_vs_VAR  = round(dm_h_var$stat, 2),
      HWI_p_vs_VAR   = signif(dm_h_var$pval, 3),
      PAY_MASE       = round(mean(ae_pay[, h, m_i]) / scale_pay, 3),
      PAY_CRPS       = round(mean(crps_pay[, h, m_i]), 1),
      PAY_DM_vs_VAR  = round(dm_p_var$stat, 2),
      PAY_p_vs_VAR   = signif(dm_p_var$pval, 3),
      stringsAsFactors = FALSE
    )
  }
}
eval_by_h <- do.call(rbind, horiz_rows)
write.csv(eval_by_h, "outputs/econometrics/table_rolling_origin_by_horizon.csv", row.names = FALSE)

# Summary Ledger (h = 1, h = 6, h = 12, and Pooled Mean h = 1..12 + DM test of Pooled CRPS vs COVID-VAR(2))
sum_rows <- list()
for (m_i in seq_len(n_mod)) {
  # Pooled origin-level average loss across h=1..12 for clean Diebold-Mariano test across N_orig origins
  pool_ae_h   <- rowMeans(ae_hwi[, , m_i]);   pool_crps_h <- rowMeans(crps_hwi[, , m_i])
  pool_ae_p   <- rowMeans(ae_pay[, , m_i]);   pool_crps_p <- rowMeans(crps_pay[, , m_i])
  dm_h_crps   <- dm_test_hln(pool_crps_h, rowMeans(crps_hwi[, , 2]), h = 6) # vs COVID-VAR(2)
  dm_p_crps   <- dm_test_hln(pool_crps_p, rowMeans(crps_pay[, , 2]), h = 6)

  sum_rows[[m_i]] <- data.frame(
    Model           = models[m_i],
    HWI_MASE_h1     = round(mean(ae_hwi[, 1, m_i])  / scale_hwi, 2),
    HWI_MASE_h12    = round(mean(ae_hwi[, 12, m_i]) / scale_hwi, 2),
    HWI_MASE_Pool   = round(mean(pool_ae_h) / scale_hwi, 2),
    HWI_CRPS_Pool   = round(mean(pool_crps_h), 1),
    HWI_DM_CRPS     = round(dm_h_crps$stat, 2),
    HWI_DM_p        = signif(dm_h_crps$pval, 3),
    PAY_MASE_h1     = round(mean(ae_pay[, 1, m_i])  / scale_pay, 2),
    PAY_MASE_h12    = round(mean(ae_pay[, 12, m_i]) / scale_pay, 2),
    PAY_MASE_Pool   = round(mean(pool_ae_p) / scale_pay, 2),
    PAY_CRPS_Pool   = round(mean(pool_crps_p), 1),
    PAY_DM_CRPS     = round(dm_p_crps$stat, 2),
    PAY_DM_p        = signif(dm_p_crps$pval, 3),
    stringsAsFactors = FALSE
  )
}
eval_summary <- do.call(rbind, sum_rows)
print(eval_summary, row.names = FALSE)
write.csv(eval_summary, "outputs/econometrics/table_rolling_origin_eval.csv", row.names = FALSE)

# Star helper for DM p-value (only when DM < 0 i.e. beats COVID-VAR(2))
dm_str <- function(stat, p) {
  st <- if (stat < 0 && p < 0.01) "^{***}" else if (stat < 0 && p < 0.05) "^{**}" else if (stat < 0 && p < 0.10) "^{*}" else ""
  sprintf("%+.2f%s", stat, st)
}

tex_tab7 <- sprintf(
"\\begin{table}[t]
\\centering
\\caption{Rolling-Origin Evaluation ($h=1\\text{--}12$, $N=73$ Origins): MASE, CRPS \\& DM Tests}
\\label{tab:rolling_eval}
\\begin{minipage}{\\columnwidth}
\\centering
\\fontsize{6.3pt}{7.7pt}\\selectfont
\\setlength{\\tabcolsep}{1.6pt}
\\begin{tabular}{lccccc}
\\toprule
\\textbf{Model} & \\textbf{MASE}$_{1}$ & \\textbf{MASE}$_{12}$ & \\textbf{MASE}$_{\\text{1-12}}$ & \\textbf{CRPS}$_{\\text{1-12}}$ & \\textbf{DM}$_{\\text{CRPS}}$ \\\\
\\midrule
\\multicolumn{6}{l}{\\textit{Panel A: JOLTS Job Openings (\\texttt{HWI}, Thousands of Openings)}} \\\\
RW-Drift Benchmark    & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
COVID-VAR(2) [Ref.]   & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & -- \\\\
\\textbf{Cointegrated VECM($r=2$)} & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%.1f}$ & $\\mathbf{%s}$ \\\\
Unscaled BVAR(2)      & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
LP-BVAR(2) [Diff]     & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
LP-BVAR(3) [Levels]   & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
\\midrule
\\multicolumn{6}{l}{\\textit{Panel B: Total Nonfarm Payrolls (\\texttt{PAYEMS}, Thousands of Persons)}} \\\\
RW-Drift Benchmark    & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
COVID-VAR(2) [Ref.]   & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & -- \\\\
\\textbf{Cointegrated VECM($r=2$)} & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%.2f}$ & $\\mathbf{%.1f}$ & $\\mathbf{%s}$ \\\\
Unscaled BVAR(2)      & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
LP-BVAR(2) [Diff]     & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
LP-BVAR(3) [Levels]   & $%.2f$ & $%.2f$ & $%.2f$ & $%.1f$ & $%s$ \\\\
\\bottomrule
\\end{tabular}
\\vskip 2pt
{\\scriptsize \\textit{Notes}: Expanding-window rolling origins spanning 2016M1--2019M2 (38 pre-COVID origins) and 2022M10--2025M8 (35 post-2022 origins); zero post-2022 hyperparameter tuning. Lower MASE and CRPS indicate superior accuracy. DM$_{\\text{CRPS}}$ reports the Harvey--Leybourne--Newbold Diebold--Mariano statistic vs.\\ COVID-VAR(2) (negative favors candidate; $^{***}p<0.01$, $^{**}p<0.05$, $^{*}p<0.10$).\\par}
\\end{minipage}
\\end{table}
",
eval_summary$HWI_MASE_h1[1], eval_summary$HWI_MASE_h12[1], eval_summary$HWI_MASE_Pool[1], eval_summary$HWI_CRPS_Pool[1], dm_str(eval_summary$HWI_DM_CRPS[1], eval_summary$HWI_DM_p[1]),
eval_summary$HWI_MASE_h1[2], eval_summary$HWI_MASE_h12[2], eval_summary$HWI_MASE_Pool[2], eval_summary$HWI_CRPS_Pool[2],
eval_summary$HWI_MASE_h1[3], eval_summary$HWI_MASE_h12[3], eval_summary$HWI_MASE_Pool[3], eval_summary$HWI_CRPS_Pool[3], dm_str(eval_summary$HWI_DM_CRPS[3], eval_summary$HWI_DM_p[3]),
eval_summary$HWI_MASE_h1[4], eval_summary$HWI_MASE_h12[4], eval_summary$HWI_MASE_Pool[4], eval_summary$HWI_CRPS_Pool[4], dm_str(eval_summary$HWI_DM_CRPS[4], eval_summary$HWI_DM_p[4]),
eval_summary$HWI_MASE_h1[5], eval_summary$HWI_MASE_h12[5], eval_summary$HWI_MASE_Pool[5], eval_summary$HWI_CRPS_Pool[5], dm_str(eval_summary$HWI_DM_CRPS[5], eval_summary$HWI_DM_p[5]),
eval_summary$HWI_MASE_h1[6], eval_summary$HWI_MASE_h12[6], eval_summary$HWI_MASE_Pool[6], eval_summary$HWI_CRPS_Pool[6], dm_str(eval_summary$HWI_DM_CRPS[6], eval_summary$HWI_DM_p[6]),
eval_summary$PAY_MASE_h1[1], eval_summary$PAY_MASE_h12[1], eval_summary$PAY_MASE_Pool[1], eval_summary$PAY_CRPS_Pool[1], dm_str(eval_summary$PAY_DM_CRPS[1], eval_summary$PAY_DM_p[1]),
eval_summary$PAY_MASE_h1[2], eval_summary$PAY_MASE_h12[2], eval_summary$PAY_MASE_Pool[2], eval_summary$PAY_CRPS_Pool[2],
eval_summary$PAY_MASE_h1[3], eval_summary$PAY_MASE_h12[3], eval_summary$PAY_MASE_Pool[3], eval_summary$PAY_CRPS_Pool[3], dm_str(eval_summary$PAY_DM_CRPS[3], eval_summary$PAY_DM_p[3]),
eval_summary$PAY_MASE_h1[4], eval_summary$PAY_MASE_h12[4], eval_summary$PAY_MASE_Pool[4], eval_summary$PAY_CRPS_Pool[4], dm_str(eval_summary$PAY_DM_CRPS[4], eval_summary$PAY_DM_p[4]),
eval_summary$PAY_MASE_h1[5], eval_summary$PAY_MASE_h12[5], eval_summary$PAY_MASE_Pool[5], eval_summary$PAY_CRPS_Pool[5], dm_str(eval_summary$PAY_DM_CRPS[5], eval_summary$PAY_DM_p[5]),
eval_summary$PAY_MASE_h1[6], eval_summary$PAY_MASE_h12[6], eval_summary$PAY_MASE_Pool[6], eval_summary$PAY_CRPS_Pool[6], dm_str(eval_summary$PAY_DM_CRPS[6], eval_summary$PAY_DM_p[6])
)
cat(tex_tab7, file = "outputs/econometrics/table7_rolling_origin_dm.tex")

# ------------------------------------------------------------------------------
# 5. Publication Figure 8: Horizon-by-Horizon CRPS Comparison (h = 1..12)
# ------------------------------------------------------------------------------
fig8_df <- rbind(
  data.frame(Horizon = eval_by_h$Horizon, CRPS = eval_by_h$HWI_CRPS, Model = eval_by_h$Model, Panel = "(a) Job Openings (HWI CRPS)"),
  data.frame(Horizon = eval_by_h$Horizon, CRPS = eval_by_h$PAY_CRPS, Model = eval_by_h$Model, Panel = "(b) Payrolls (PAYEMS CRPS)")
)
fig8_df$Model <- factor(fig8_df$Model, levels = c(
  "LP-BVAR(3) [Levels]", "LP-BVAR(2) [Diff]", "VECM(r=2)", "COVID-VAR(2)", "Unscaled-BVAR(2)", "RW-Drift"
))

build_fig8 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  lw <- if (is_web) 0.85 else 0.62
  ggplot(fig8_df, aes(x = Horizon, y = CRPS, color = Model, linetype = Model)) +
    geom_line(linewidth = lw) +
    geom_point(size = if (is_web) 1.5 else 0.9) +
    facet_wrap(~Panel, scales = "free_y", ncol = if (is_web) 2 else 1) +
    scale_color_manual(name = NULL, values = c(
      "LP-BVAR(3) [Levels]" = PAL$navy,
      "LP-BVAR(2) [Diff]"   = PAL$amber,
      "VECM(r=2)"           = PAL$teal,
      "COVID-VAR(2)"        = PAL$red,
      "Unscaled-BVAR(2)"    = "#8E44AD",
      "RW-Drift"            = PAL$gray
    )) +
    scale_linetype_manual(name = NULL, values = c(
      "LP-BVAR(3) [Levels]" = "solid",
      "LP-BVAR(2) [Diff]"   = "dashed",
      "VECM(r=2)"           = "longdash",
      "COVID-VAR(2)"        = "dotdash",
      "Unscaled-BVAR(2)"    = "twodash",
      "RW-Drift"            = "dotted"
    )) +
    scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
    guides(
      color    = guide_legend(nrow = if (is_web) 2 else 3, byrow = TRUE),
      linetype = guide_legend(nrow = if (is_web) 2 else 3, byrow = TRUE)
    ) +
    labs(x = "Forecast Horizon (Months)", y = "Mean CRPS (Lower = Better)") +
    th
}

p8_print <- build_fig8(is_web = FALSE)
ggsave("outputs/figures/fig8_rolling_origin_horizons.pdf", p8_print, width = 3.5, height = 3.7)
ggsave("outputs/figures/fig8_rolling_origin_horizons.png", p8_print, width = 3.5, height = 3.7, dpi = 300)

p8_web <- build_fig8(is_web = TRUE)
ggsave("outputs/figures/fig8_web_rolling_horizons.png", p8_web, width = 8.4, height = 4.1, dpi = 300)
file.copy("outputs/figures/fig8_web_rolling_horizons.png", "docs/assets/images/fig8_rolling_origin_horizons.png", overwrite = TRUE)

cat("\nStep 5 (Rolling-Origin Multi-Horizon Evaluation) completed successfully.\n")
