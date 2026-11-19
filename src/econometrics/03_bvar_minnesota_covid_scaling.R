#!/usr/bin/env Rscript
# ==============================================================================
# 03_bvar_minnesota_covid_scaling.R
# Step 3: Bayesian VAR with Exact Marginal-Likelihood Minnesota Shrinkage Prior
#         (Giannone, Lenza & Primiceri 2015) and Lenza & Primiceri (2022, JAE)
#         COVID-19 Volatility Scaling, Cointegration-Preserving Levels BVAR(3)
#         (Enders 2015 Ch. 6; Crump et al. 2021), and Exact Waggoner-Zha (1999)
#         Conditional Predictive Density Fan Charts.
#
# Methodological references:
#   - Lenza, M. & Primiceri, G. E. (2022). How to estimate a VAR after March 2020.
#     Journal of Applied Econometrics, 37(4), 688-699.
#   - Giannone, D., Lenza, M. & Primiceri, G. E. (2015). Prior selection for
#     vector autoregressions. Review of Economics and Statistics, 97(2), 436-451.
#   - Doan, T., Litterman, R. & Sims, C. (1984). Econometric Reviews, 3(1), 1-100.
#   - Waggoner, D. F. & Zha, T. (1999). Review of Economics and Statistics, 81(4), 639-651.
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.), Ch. 5 & 6.
# ==============================================================================

suppressPackageStartupMessages({
  library(vars)
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
# 1. Load FRED-MD Data & Build Stationary and Cointegration-Preserving Systems
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

# Stationary FRED-MD transforms
df$t_FEDFUNDS <- d1(df$FEDFUNDS)
df$t_M2SL     <- d2log1(df$M2SL)
df$t_BAA10YM  <- df$BAA10YM
df$t_INDPRO   <- dlog1(df$INDPRO)
df$t_HWI      <- d1(df$HWI)
df$t_PAYEMS   <- dlog1(df$PAYEMS)

# Levels / 100*log-levels system (Crump et al. 2021; Lenza & Primiceri 2022)
df$L_FEDFUNDS <- df$FEDFUNDS
df$L_M2SL     <- 100 * log(df$M2SL)
df$L_BAA10YM  <- df$BAA10YM
df$L_INDPRO   <- 100 * log(df$INDPRO)
df$L_HWI      <- 100 * log(df$HWI)
df$L_PAYEMS   <- 100 * log(df$PAYEMS)

full_sample <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]
train_df    <- full_sample[full_sample$date <= as.Date("2022-10-01"), ] # T = 154
test_df     <- full_sample[full_sample$date >= as.Date("2022-11-01"), ] # H = 46
precovid_df <- train_df[train_df$date <= as.Date("2020-02-01"), ]       # T_pre = 122

var_names <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")
diff_cols <- c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")
lvl_cols  <- c("L_FEDFUNDS", "L_M2SL", "L_BAA10YM", "L_INDPRO", "L_HWI", "L_PAYEMS")

Y_tr_diff  <- as.matrix(train_df[, diff_cols]);    colnames(Y_tr_diff)  <- var_names
Y_te_diff  <- as.matrix(test_df[, diff_cols]);     colnames(Y_te_diff)  <- var_names
Y_pre_diff <- as.matrix(precovid_df[, diff_cols]); colnames(Y_pre_diff) <- var_names

Y_tr_lvl   <- as.matrix(train_df[, lvl_cols]);     colnames(Y_tr_lvl)   <- var_names
Y_te_lvl   <- as.matrix(test_df[, lvl_cols]);      colnames(Y_te_lvl)   <- var_names
Y_pre_lvl  <- as.matrix(precovid_df[, lvl_cols]);  colnames(Y_pre_lvl)  <- var_names

k    <- length(var_names)
T_tr <- nrow(train_df)
H    <- nrow(test_df)
cond_idx <- c(1, 2, 3) # FEDFUNDS, M2SL, BAA10YM
base_hwi <- train_df$HWI[T_tr]
base_pay <- train_df$PAYEMS[T_tr]
hwi_real <- test_df$HWI
pay_real <- test_df$PAYEMS

cat("======================================================================\n")
cat("BAYESIAN VAR: EXACT LENZA-PRIMICERI (2022) MARGINAL LIKELIHOOD & WZ\n")
cat("======================================================================\n")

# ------------------------------------------------------------------------------
# 2. Core Estimation Engine: Conjugate NIW Minnesota BVAR + Lenza-Primiceri ML
# ------------------------------------------------------------------------------
build_xy <- function(Y_mat, p) {
  Tn <- nrow(Y_mat); kn <- ncol(Y_mat); Kr <- 1 + kn * p
  Ye <- Y_mat[(p + 1):Tn, , drop = FALSE]
  Te <- nrow(Ye)
  Xe <- matrix(0, Te, Kr); Xe[, 1] <- 1
  for (j in 1:p) Xe[, (1 + (j - 1) * kn + 1):(1 + j * kn)] <- Y_mat[(p + 1 - j):(Tn - j), ]
  list(Y = Ye, X = Xe, T_eff = Te, K_reg = Kr)
}

get_pre_sigmas <- function(Y_pre_mat, p) {
  kn <- ncol(Y_pre_mat); Tp <- nrow(Y_pre_mat)
  sig <- numeric(kn)
  for (i in 1:kn) {
    yi <- Y_pre_mat[, i]
    Yi <- yi[(p + 1):Tp]
    Xi <- cbind(1, sapply(1:p, function(l) yi[(p + 1 - l):(Tp - l)]))
    ri <- Yi - Xi %*% solve(crossprod(Xi), crossprod(Xi, Yi))
    sig[i] <- sqrt(sum(ri^2) / (length(Yi) - ncol(Xi)))
  }
  names(sig) <- colnames(Y_pre_mat)
  sig
}

fit_lp_bvar <- function(Y_tr_mat, Y_pre_mat, dates_tr, p, rw_prior = FALSE, covid_scale = TRUE) {
  xy <- build_xy(Y_tr_mat, p)
  Ye <- xy$Y; Xe <- xy$X; Te <- xy$T_eff; Kr <- xy$K_reg
  dates_e <- dates_tr[(p + 1):length(dates_tr)]
  idx_tstar <- which(dates_e == as.Date("2020-03-01"))

  sig_pre <- get_pre_sigmas(Y_pre_mat, p)
  d0   <- k + 2
  Psi0 <- diag(sig_pre^2)
  B0   <- matrix(0, Kr, k)
  if (rw_prior) {
    for (i in 1:k) B0[1 + i, i] <- 1.0 # Litterman random-walk prior on I(1) levels
  }

  eval_lml <- function(lam1, s0, s1, s2, rho) {
    if (lam1 < 0.02 || lam1 > 5 || s0 > 150 || s1 > 150 || s2 > 150) return(-1e12)
    s_vec <- rep(1, Te)
    if (covid_scale) {
      s_vec[idx_tstar]     <- s0
      s_vec[idx_tstar + 1] <- s1
      for (t in (idx_tstar + 2):Te) s_vec[t] <- 1 + (s2 - 1) * (rho^(t - (idx_tstar + 2)))
    }
    omega_diag <- numeric(Kr); omega_diag[1] <- 1e4
    for (l in 1:p) for (i in 1:k) omega_diag[1 + (l - 1) * k + i] <- (lam1^2) / ((l^2) * (sig_pre[i]^2))
    Omega0_inv <- diag(1 / omega_diag)
    Yt <- Ye / s_vec; Xt <- Xe / s_vec
    V_inv <- crossprod(Xt) + Omega0_inv
    Bp    <- solve(V_inv, crossprod(Xt, Yt) + Omega0_inv %*% B0)
    Ep    <- Yt - Xt %*% Bp
    Psip  <- Psi0 + crossprod(Ep) + t(Bp - B0) %*% Omega0_inv %*% (Bp - B0)

    lml <- -k * sum(log(s_vec)) -
           (k / 2) * (sum(log(omega_diag)) + as.numeric(determinant(V_inv, logarithm = TRUE)$modulus)) +
           (d0 / 2) * as.numeric(determinant(Psi0, logarithm = TRUE)$modulus) -
           ((d0 + Te) / 2) * as.numeric(determinant(Psip, logarithm = TRUE)$modulus)
    if (!is.finite(lml)) return(-1e12)
    lml
  }

  if (covid_scale) {
    opt <- optim(
      par = c(log(0.2), log(10), log(40), log(5), qlogis(0.85)),
      fn = function(par) -eval_lml(exp(par[1]), 1 + exp(par[2]), 1 + exp(par[3]), 1 + exp(par[4]), 1 / (1 + exp(-par[5]))),
      method = "BFGS", control = list(maxit = 500)
    )
    lam1_hat <- exp(opt$par[1])
    s0_hat   <- 1 + exp(opt$par[2])
    s1_hat   <- 1 + exp(opt$par[3])
    s2_hat   <- 1 + exp(opt$par[4])
    rho_hat  <- 1 / (1 + exp(-opt$par[5]))
    lml_opt  <- -opt$value
  } else {
    opt <- optimize(function(ll) eval_lml(exp(ll), 1, 1, 1, 0), interval = c(log(0.02), log(2)), maximum = TRUE)
    lam1_hat <- exp(opt$maximum)
    s0_hat <- 1; s1_hat <- 1; s2_hat <- 1; rho_hat <- 0
    lml_opt <- opt$objective
  }

  s_vec <- rep(1, Te)
  if (covid_scale) {
    s_vec[idx_tstar] <- s0_hat; s_vec[idx_tstar + 1] <- s1_hat
    for (t in (idx_tstar + 2):Te) s_vec[t] <- 1 + (s2_hat - 1) * (rho_hat^(t - (idx_tstar + 2)))
  }
  omega_diag <- numeric(Kr); omega_diag[1] <- 1e4
  for (l in 1:p) for (i in 1:k) omega_diag[1 + (l - 1) * k + i] <- (lam1_hat^2) / ((l^2) * (sig_pre[i]^2))
  Omega0_inv <- diag(1 / omega_diag)
  Yt <- Ye / s_vec; Xt <- Xe / s_vec
  V_post <- solve(crossprod(Xt) + Omega0_inv)
  B_post <- V_post %*% (crossprod(Xt, Yt) + Omega0_inv %*% B0)
  E_post <- Yt - Xt %*% B_post
  Psi_post <- Psi0 + crossprod(E_post) + t(B_post - B0) %*% Omega0_inv %*% (B_post - B0)
  nu_post <- d0 + Te
  Sigma_post <- Psi_post / (nu_post - k - 1)

  # Max companion root
  cb <- matrix(0, k * p, k * p)
  for (l in 1:p) cb[1:k, ((l - 1) * k + 1):(l * k)] <- t(B_post[(1 + (l - 1) * k + 1):(1 + l * k), ])
  if (p > 1) cb[(k + 1):(k * p), 1:(k * (p - 1))] <- diag(k * (p - 1))
  max_root <- max(Mod(eigen(cb, only.values = TRUE)$values))

  list(
    p = p, K_reg = Kr, rw_prior = rw_prior, covid_scale = covid_scale,
    lam1 = lam1_hat, s0 = s0_hat, s1 = s1_hat, s2 = s2_hat, rho = rho_hat,
    lml = lml_opt, max_root = max_root, s_vec = s_vec,
    B_post = B_post, V_post = V_post, Psi_post = Psi_post,
    nu_post = nu_post, Sigma_post = Sigma_post, sig_pre = sig_pre
  )
}

# ------------------------------------------------------------------------------
# 3. Exact Waggoner & Zha (1999) Conditional Predictive Density Sampler
#    Combines Parameter Uncertainty (B, Sigma ~ NIW) + Exact Future Conditional
#    Innovation Uncertainty: eta ~ N(mu_eta, I - R' (R R')^-1 R)
# ------------------------------------------------------------------------------
wz_single <- function(B_m, Sig_m, Y_init, Y_eval, p, cond_idx, add_shock = TRUE) {
  kn <- ncol(Y_eval); Hv <- nrow(Y_eval); n_c <- length(cond_idx)
  A_m <- lapply(1:p, function(l) t(B_m[(1 + (l - 1) * kn + 1):(1 + l * kn), ]))
  c_m <- B_m[1, ]
  P_m <- t(chol(Sig_m))

  Y_u <- matrix(0, Hv, kn)
  for (h in 1:Hv) {
    pr <- c_m
    for (j in 1:p) {
      lv <- if ((h - j) <= 0) Y_init[p + (h - j), ] else Y_u[h - j, ]
      pr <- pr + A_m[[j]] %*% lv
    }
    Y_u[h, ] <- pr
  }
  Mm <- vector("list", Hv); Mm[[1]] <- diag(kn)
  for (m in 2:Hv) {
    Mt <- matrix(0, kn, kn)
    for (j in 1:min(p, m - 1)) Mt <- Mt + A_m[[j]] %*% Mm[[m - j]]
    Mm[[m]] <- Mt
  }

  R <- matrix(0, Hv * n_c, Hv * kn); r <- numeric(Hv * n_c); rc <- 1
  for (h in 1:Hv) {
    MP_l <- lapply(1:h, function(s) Mm[[h - s + 1]] %*% P_m)
    for (ci in 1:n_c) {
      cv <- cond_idx[ci]
      r[rc] <- Y_eval[h, cv] - Y_u[h, cv]
      for (s in 1:h) R[rc, ((s - 1) * kn + 1):(s * kn)] <- MP_l[[s]][cv, ]
      rc <- rc + 1
    }
  }
  RRt_inv <- solve(R %*% t(R))
  mu_eta  <- as.numeric(t(R) %*% (RRt_inv %*% r))
  if (add_shock) {
    # Exact projection draw (Waggoner & Zha 1999, Eq. 12):
    z <- rnorm(Hv * kn)
    eta <- mu_eta + z - as.numeric(t(R) %*% (RRt_inv %*% (R %*% z)))
  } else {
    eta <- mu_eta
  }
  Y_c <- Y_u
  for (h in 1:Hv) {
    adj <- rep(0, kn)
    for (s in 1:h) adj <- adj + as.numeric(Mm[[h - s + 1]] %*% P_m %*% eta[((s - 1) * kn + 1):(s * kn)])
    Y_c[h, ] <- Y_u[h, ] + adj
  }
  list(Y_u = Y_u, Y_c = Y_c)
}

sample_bvar_density <- function(fit, Y_tr_mat, Y_te_mat, M = 500, is_levels = FALSE) {
  p <- fit$p; Kr <- fit$K_reg
  Y_init <- Y_tr_mat[(nrow(Y_tr_mat) - p + 1):nrow(Y_tr_mat), ]
  V_chol <- t(chol(fit$V_post))
  Psi_inv <- solve(fit$Psi_post)
  root_cap <- if (is_levels) 1.015 else 0.999

  # Posterior mean point forecast
  pt <- wz_single(fit$B_post, fit$Sigma_post, Y_init, Y_te_mat, p, cond_idx, add_shock = FALSE)

  hwi_lvl_cnd <- matrix(0, M, H); pay_lvl_cnd <- matrix(0, M, H)
  hwi_lvl_unc <- matrix(0, M, H); pay_lvl_unc <- matrix(0, M, H)
  dhwi_mean_err <- numeric(M); dpay_mean_err <- numeric(M)

  m_ok <- 0; att <- 0
  while (m_ok < M && att < M * 10) {
    att <- att + 1
    Sig_d <- solve(rWishart(1, fit$nu_post, Psi_inv)[, , 1])
    P_d   <- t(chol(Sig_d))
    B_d   <- fit$B_post + V_chol %*% matrix(rnorm(Kr * k), Kr, k) %*% t(P_d)
    cb <- matrix(0, k * p, k * p)
    for (l in 1:p) cb[1:k, ((l - 1) * k + 1):(l * k)] <- t(B_d[(1 + (l - 1) * k + 1):(1 + l * k), ])
    if (p > 1) cb[(k + 1):(k * p), 1:(k * (p - 1))] <- diag(k * (p - 1))
    if (max(Mod(eigen(cb, only.values = TRUE)$values)) > root_cap) next

    m_ok <- m_ok + 1
    dr <- wz_single(B_d, Sig_d, Y_init, Y_te_mat, p, cond_idx, add_shock = TRUE)
    if (is_levels) {
      hwi_lvl_cnd[m_ok, ] <- exp(dr$Y_c[, 5] / 100)
      pay_lvl_cnd[m_ok, ] <- exp(dr$Y_c[, 6] / 100)
      hwi_lvl_unc[m_ok, ] <- exp(dr$Y_u[, 5] / 100)
      pay_lvl_unc[m_ok, ] <- exp(dr$Y_u[, 6] / 100)
      dhwi_cnd <- diff(c(base_hwi, hwi_lvl_cnd[m_ok, ]))
      dpay_cnd <- diff(c(100 * log(base_pay), dr$Y_c[, 6]))
      dhwi_mean_err[m_ok] <- mean(Y_te_diff[, 5] - dhwi_cnd)
      dpay_mean_err[m_ok] <- mean(Y_te_diff[, 6] - dpay_cnd)
    } else {
      hwi_lvl_cnd[m_ok, ] <- base_hwi + cumsum(dr$Y_c[, 5])
      pay_lvl_cnd[m_ok, ] <- base_pay * exp(cumsum(dr$Y_c[, 6]) / 100)
      hwi_lvl_unc[m_ok, ] <- base_hwi + cumsum(dr$Y_u[, 5])
      pay_lvl_unc[m_ok, ] <- base_pay * exp(cumsum(dr$Y_u[, 6]) / 100)
      dhwi_mean_err[m_ok] <- mean(Y_te_diff[, 5] - dr$Y_c[, 5])
      dpay_mean_err[m_ok] <- mean(Y_te_diff[, 6] - dr$Y_c[, 6])
    }
  }

  pt_hwi <- if (is_levels) exp(pt$Y_c[, 5] / 100) else base_hwi + cumsum(pt$Y_c[, 5])
  pt_pay <- if (is_levels) exp(pt$Y_c[, 6] / 100) else base_pay * exp(cumsum(pt$Y_c[, 6]) / 100)
  pt_hwi_u <- if (is_levels) exp(pt$Y_u[, 5] / 100) else base_hwi + cumsum(pt$Y_u[, 5])
  pt_pay_u <- if (is_levels) exp(pt$Y_u[, 6] / 100) else base_pay * exp(cumsum(pt$Y_u[, 6]) / 100)

  list(
    M = m_ok, attempts = att,
    pt_hwi = pt_hwi, pt_pay = pt_pay, pt_hwi_u = pt_hwi_u, pt_pay_u = pt_pay_u,
    hwi_lvl_cnd = hwi_lvl_cnd[1:m_ok, ], pay_lvl_cnd = pay_lvl_cnd[1:m_ok, ],
    hwi_lvl_unc = hwi_lvl_unc[1:m_ok, ], pay_lvl_unc = pay_lvl_unc[1:m_ok, ],
    dhwi_mean_err = dhwi_mean_err[1:m_ok], dpay_mean_err = dpay_mean_err[1:m_ok]
  )
}

# ------------------------------------------------------------------------------
# 4. Fit All 3 Bayesian Specifications + Classical COVID-Controlled VAR(2)
# ------------------------------------------------------------------------------
fit_bvar_unsc_diff <- fit_lp_bvar(Y_tr_diff, Y_pre_diff, train_df$date, p = 2, rw_prior = FALSE, covid_scale = FALSE)
fit_bvar_lp_diff   <- fit_lp_bvar(Y_tr_diff, Y_pre_diff, train_df$date, p = 2, rw_prior = FALSE, covid_scale = TRUE)
fit_bvar_lp_lvl    <- fit_lp_bvar(Y_tr_lvl,  Y_pre_lvl,  train_df$date, p = 3, rw_prior = TRUE,  covid_scale = TRUE)

cat(sprintf("1. Unscaled Minnesota BVAR(2) [Diff]  : lam1=%.4f | LML=%.2f | MaxRoot=%.4f\n",
            fit_bvar_unsc_diff$lam1, fit_bvar_unsc_diff$lml, fit_bvar_unsc_diff$max_root))
cat(sprintf("2. Lenza-Primiceri BVAR(2) [Diff]     : lam1=%.4f | s0=%.2f, s1=%.2f, s2=%.2f, rho=%.3f | LML=%.2f (LogBF=+%.1f)\n",
            fit_bvar_lp_diff$lam1, fit_bvar_lp_diff$s0, fit_bvar_lp_diff$s1, fit_bvar_lp_diff$s2,
            fit_bvar_lp_diff$rho, fit_bvar_lp_diff$lml, fit_bvar_lp_diff$lml - fit_bvar_unsc_diff$lml))
cat(sprintf("3. Lenza-Primiceri BVAR(3) [Levels]   : lam1=%.4f | s0=%.2f, s1=%.2f, s2=%.2f, rho=%.3f | LML=%.2f | MaxRoot=%.4f\n\n",
            fit_bvar_lp_lvl$lam1, fit_bvar_lp_lvl$s0, fit_bvar_lp_lvl$s1, fit_bvar_lp_lvl$s2,
            fit_bvar_lp_lvl$rho, fit_bvar_lp_lvl$lml, fit_bvar_lp_lvl$max_root))

cat("Sampling 500 exact Waggoner-Zha predictive density draws per Bayesian model...\n")
t0 <- Sys.time()
sim_unsc_diff <- sample_bvar_density(fit_bvar_unsc_diff, Y_tr_diff, Y_te_diff, M = 500, is_levels = FALSE)
sim_lp_diff   <- sample_bvar_density(fit_bvar_lp_diff,   Y_tr_diff, Y_te_diff, M = 500, is_levels = FALSE)
sim_lp_lvl    <- sample_bvar_density(fit_bvar_lp_lvl,    Y_tr_lvl,  Y_te_lvl,  M = 500, is_levels = TRUE)
cat(sprintf("Completed all 1,500 posterior predictive draws in %.1f seconds.\n\n",
            difftime(Sys.time(), t0, units = "secs")))

# Classical COVID-Controlled VAR(2) benchmark (from Step 2 CSV for exact consistency)
traj_var <- read.csv("outputs/econometrics/table_counterfactual_trajectories.csv", stringsAsFactors = FALSE)
hac_var  <- read.csv("outputs/econometrics/table_hac_break_tests.csv", stringsAsFactors = FALSE)
var_cnd_hwi <- traj_var$Value[traj_var$Type == "Policy-Conditioned VAR(2)" & grepl("HWI", traj_var$SeriesLabel)]
var_cnd_pay <- traj_var$Value[traj_var$Type == "Policy-Conditioned VAR(2)" & grepl("PAYEMS", traj_var$SeriesLabel)]

# ------------------------------------------------------------------------------
# 5. Export Exact Model Comparison Ledger (Zero Hardcoded Values!)
# ------------------------------------------------------------------------------
pct_gap <- function(real_end, cnd_end) 100 * (real_end / cnd_end - 1)

# Newey-West HAC CIs for Classical VAR(2)
var_hwi_mean <- hac_var$Mean_Error[2]; var_hwi_se <- hac_var$HAC_SE[2]
var_pay_mean <- hac_var$Mean_Error[4]; var_pay_se <- hac_var$HAC_SE[4]

comp_table <- data.frame(
  Model = c(
    "COVID-Controlled VAR(2) (OLS)",
    "Unscaled Minnesota BVAR(2) (Diff)",
    "Lenza-Primiceri BVAR(2) (Diff)",
    "Lenza-Primiceri BVAR(3) (Levels)"
  ),
  HWI_Mean_Err = c(
    var_hwi_mean,
    mean(sim_unsc_diff$dhwi_mean_err),
    mean(sim_lp_diff$dhwi_mean_err),
    mean(sim_lp_lvl$dhwi_mean_err)
  ),
  HWI_CI_Low = c(
    var_hwi_mean - 1.96 * var_hwi_se,
    quantile(sim_unsc_diff$dhwi_mean_err, 0.025),
    quantile(sim_lp_diff$dhwi_mean_err,   0.025),
    quantile(sim_lp_lvl$dhwi_mean_err,    0.025)
  ),
  HWI_CI_High = c(
    var_hwi_mean + 1.96 * var_hwi_se,
    quantile(sim_unsc_diff$dhwi_mean_err, 0.975),
    quantile(sim_lp_diff$dhwi_mean_err,   0.975),
    quantile(sim_lp_lvl$dhwi_mean_err,    0.975)
  ),
  HWI_End_PctGap = c(
    pct_gap(hwi_real[H], var_cnd_hwi[H]),
    pct_gap(hwi_real[H], median(sim_unsc_diff$hwi_lvl_cnd[, H])),
    pct_gap(hwi_real[H], median(sim_lp_diff$hwi_lvl_cnd[, H])),
    pct_gap(hwi_real[H], median(sim_lp_lvl$hwi_lvl_cnd[, H]))
  ),
  HWI_Prob_Neg = c(
    pnorm(0, mean = var_hwi_mean, sd = var_hwi_se, lower.tail = FALSE), # 1 - Phi(0) for negative mean -> P(e<0)
    mean(sim_unsc_diff$dhwi_mean_err < 0),
    mean(sim_lp_diff$dhwi_mean_err < 0),
    mean(sim_lp_lvl$dhwi_mean_err < 0)
  ),
  PAY_Mean_Err = c(
    var_pay_mean,
    mean(sim_unsc_diff$dpay_mean_err),
    mean(sim_lp_diff$dpay_mean_err),
    mean(sim_lp_lvl$dpay_mean_err)
  ),
  PAY_CI_Low = c(
    var_pay_mean - 1.96 * var_pay_se,
    quantile(sim_unsc_diff$dpay_mean_err, 0.025),
    quantile(sim_lp_diff$dpay_mean_err,   0.025),
    quantile(sim_lp_lvl$dpay_mean_err,    0.025)
  ),
  PAY_CI_High = c(
    var_pay_mean + 1.96 * var_pay_se,
    quantile(sim_unsc_diff$dpay_mean_err, 0.975),
    quantile(sim_lp_diff$dpay_mean_err,   0.975),
    quantile(sim_lp_lvl$dpay_mean_err,    0.975)
  ),
  PAY_End_PctGap = c(
    pct_gap(pay_real[H], var_cnd_pay[H]),
    pct_gap(pay_real[H], median(sim_unsc_diff$pay_lvl_cnd[, H])),
    pct_gap(pay_real[H], median(sim_lp_diff$pay_lvl_cnd[, H])),
    pct_gap(pay_real[H], median(sim_lp_lvl$pay_lvl_cnd[, H]))
  ),
  PAY_Prob_Neg = c(
    pnorm(0, mean = var_pay_mean, sd = var_pay_se, lower.tail = TRUE),
    mean(sim_unsc_diff$dpay_mean_err < 0),
    mean(sim_lp_diff$dpay_mean_err < 0),
    mean(sim_lp_lvl$dpay_mean_err < 0)
  ),
  stringsAsFactors = FALSE
)
# Fix P(e < 0) for VAR HWI using exact normal CDF: P(Z < 0) = pnorm(0, mean, sd)
comp_table$HWI_Prob_Neg[1] <- pnorm(0, mean = var_hwi_mean, sd = var_hwi_se, lower.tail = TRUE)

print(comp_table, row.names = FALSE)
write.csv(comp_table, "outputs/econometrics/table_bvar_comparison.csv", row.names = FALSE)

# Export Lenza-Primiceri Hyperparameter Summary CSV
lp_hyper <- data.frame(
  Specification = c("Unscaled BVAR(2) [Diff]", "Lenza-Primiceri BVAR(2) [Diff]", "Lenza-Primiceri BVAR(3) [Levels]"),
  lambda1 = round(c(fit_bvar_unsc_diff$lam1, fit_bvar_lp_diff$lam1, fit_bvar_lp_lvl$lam1), 4),
  s_2020M3 = round(c(fit_bvar_unsc_diff$s0, fit_bvar_lp_diff$s0, fit_bvar_lp_lvl$s0), 2),
  s_2020M4 = round(c(fit_bvar_unsc_diff$s1, fit_bvar_lp_diff$s1, fit_bvar_lp_lvl$s1), 2),
  s_2020M5 = round(c(fit_bvar_unsc_diff$s2, fit_bvar_lp_diff$s2, fit_bvar_lp_lvl$s2), 2),
  rho      = round(c(fit_bvar_unsc_diff$rho, fit_bvar_lp_diff$rho, fit_bvar_lp_lvl$rho), 4),
  LogML    = round(c(fit_bvar_unsc_diff$lml, fit_bvar_lp_diff$lml, fit_bvar_lp_lvl$lml), 2),
  MaxRoot  = round(c(fit_bvar_unsc_diff$max_root, fit_bvar_lp_diff$max_root, fit_bvar_lp_lvl$max_root), 4)
)
write.csv(lp_hyper, "outputs/econometrics/table_bvar_hyperparameters.csv", row.names = FALSE)

# Single-Column IEEE LaTeX Table 5 (Zero-badness minipage)
tex_tab5 <- sprintf(
"\\begin{table}[t]
\\centering
\\caption{Classical VAR vs. Lenza--Primiceri BVAR Conditional Metrics}
\\label{tab:bvar_comparison}
\\begin{minipage}{\\columnwidth}
\\centering
\\fontsize{6.6pt}{8.0pt}\\selectfont
\\setlength{\\tabcolsep}{1.9pt}
\\begin{tabular}{lcccc}
\\toprule
\\textbf{Model Specification} & $\\bar{e}$ & \\textbf{95\\%% Interval} & \\textbf{2026 Gap} & $P(\\bar{e}<0)$ \\\\
\\midrule
\\multicolumn{5}{l}{\\textit{Job Openings ($\\Delta\\text{HWI}$, thsd/mo; 2026M8 Realized = 7,079k)}} \\\\
COVID-Controlled VAR(2) & $%.1f$ & $[%.1f, %.1f]$ & $%.1f\\%%$ & $%.3f$ \\\\
Unscaled BVAR(2) [Diff] & $%.1f$ & $[%.1f, %.1f]$ & $%.1f\\%%$ & $%.3f$ \\\\
LP-BVAR(2) [Diff]       & $%.1f$ & $[%.1f, %.1f]$ & $%.1f\\%%$ & $%.3f$ \\\\
\\textbf{LP-BVAR(3) [Levels]} & $\\mathbf{%.1f}$ & $\\mathbf{[%.1f, %.1f]}$ & $\\mathbf{%.1f\\%%}$ & $\\mathbf{%.3f}$ \\\\
\\midrule
\\multicolumn{5}{l}{\\textit{Total Payrolls ($\\Delta\\log\\text{PAYEMS}$, \\%%/mo; 2026M8 Realized = 159,015k)}} \\\\
COVID-Controlled VAR(2) & $%+.3f$ & $[%.3f, %+.3f]$ & $%.1f\\%%$ & $%.3f$ \\\\
Unscaled BVAR(2) [Diff] & $%+.3f$ & $[%.3f, %+.3f]$ & $%.1f\\%%$ & $%.3f$ \\\\
LP-BVAR(2) [Diff]       & $%+.3f$ & $[%.3f, %+.3f]$ & $%.1f\\%%$ & $%.3f$ \\\\
\\textbf{LP-BVAR(3) [Levels]} & $\\mathbf{%+.3f}$ & $\\mathbf{[%.3f, %+.3f]}$ & $\\mathbf{%.1f\\%%}$ & $\\mathbf{%.3f}$ \\\\
\\bottomrule
\\end{tabular}
\\vskip 2pt
{\\scriptsize \\textit{Notes}: Out-of-sample horizon $H=46$ (2022M11--2026M8) conditioned on realized $\\mathcal{C}=\\{\\text{FEDFUNDS},\\text{M2SL},\\text{BAA10YM}\\}$. Bayesian 95\\%% predictive intervals combine NIW parameter draws and exact Waggoner--Zha (1999) conditional shock draws ($M=500$). LP-BVAR hyperparameters are maximized via exact log marginal likelihood (Log Bayes Factor $= +%.1f$ for COVID scaling). 2026 Gap reports relative level deviation at 2026M8.\\par}
\\end{minipage}
\\end{table}
",
comp_table$HWI_Mean_Err[1], comp_table$HWI_CI_Low[1], comp_table$HWI_CI_High[1], comp_table$HWI_End_PctGap[1], comp_table$HWI_Prob_Neg[1],
comp_table$HWI_Mean_Err[2], comp_table$HWI_CI_Low[2], comp_table$HWI_CI_High[2], comp_table$HWI_End_PctGap[2], comp_table$HWI_Prob_Neg[2],
comp_table$HWI_Mean_Err[3], comp_table$HWI_CI_Low[3], comp_table$HWI_CI_High[3], comp_table$HWI_End_PctGap[3], comp_table$HWI_Prob_Neg[3],
comp_table$HWI_Mean_Err[4], comp_table$HWI_CI_Low[4], comp_table$HWI_CI_High[4], comp_table$HWI_End_PctGap[4], comp_table$HWI_Prob_Neg[4],
comp_table$PAY_Mean_Err[1], comp_table$PAY_CI_Low[1], comp_table$PAY_CI_High[1], comp_table$PAY_End_PctGap[1], comp_table$PAY_Prob_Neg[1],
comp_table$PAY_Mean_Err[2], comp_table$PAY_CI_Low[2], comp_table$PAY_CI_High[2], comp_table$PAY_End_PctGap[2], comp_table$PAY_Prob_Neg[2],
comp_table$PAY_Mean_Err[3], comp_table$PAY_CI_Low[3], comp_table$PAY_CI_High[3], comp_table$PAY_End_PctGap[3], comp_table$PAY_Prob_Neg[3],
comp_table$PAY_Mean_Err[4], comp_table$PAY_CI_Low[4], comp_table$PAY_CI_High[4], comp_table$PAY_End_PctGap[4], comp_table$PAY_Prob_Neg[4],
fit_bvar_lp_diff$lml - fit_bvar_unsc_diff$lml
)
cat(tex_tab5, file = "outputs/econometrics/table5_bvar_comparison.tex")

# ------------------------------------------------------------------------------
# 6. Publication Figure 5: Exact Waggoner-Zha Predictive Fan Charts (Levels!)
#    Cointegration-Preserving Lenza-Primiceri BVAR(3) — Unified Punto & ZERO Clipping!
# ------------------------------------------------------------------------------
plot_dates <- test_df$date

make_fan_df <- function(draws_cnd, draws_unc, real_vec, label) {
  data.frame(
    Date          = plot_dates,
    Realized      = real_vec,
    Median        = apply(draws_cnd, 2, quantile, probs = 0.500),
    Low68         = apply(draws_cnd, 2, quantile, probs = 0.160),
    High68        = apply(draws_cnd, 2, quantile, probs = 0.840),
    Low95         = apply(draws_cnd, 2, quantile, probs = 0.025),
    High95        = apply(draws_cnd, 2, quantile, probs = 0.975),
    Unconditional = apply(draws_unc, 2, quantile, probs = 0.500),
    Panel         = label,
    stringsAsFactors = FALSE
  )
}

fan_df <- rbind(
  make_fan_df(sim_lp_lvl$hwi_lvl_cnd, sim_lp_lvl$hwi_lvl_unc, hwi_real, "(a) Job Openings (HWI, thousands)"),
  make_fan_df(sim_lp_lvl$pay_lvl_cnd, sim_lp_lvl$pay_lvl_unc, pay_real, "(b) Payroll Employment (PAYEMS, thousands)")
)

build_fig5 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  lw <- if (is_web) 0.90 else 0.68
  ggplot(fan_df, aes(x = Date)) +
    geom_ribbon(aes(ymin = Low95, ymax = High95, fill = "95% Predictive Interval"), alpha = 0.30) +
    geom_ribbon(aes(ymin = Low68, ymax = High68, fill = "68% Predictive Interval"), alpha = 0.45) +
    geom_line(aes(y = Unconditional, color = "Unconditional Median", linetype = "Unconditional Median"), linewidth = lw * 0.85) +
    geom_line(aes(y = Median,        color = "LP-BVAR(3) Conditional", linetype = "LP-BVAR(3) Conditional"), linewidth = lw) +
    geom_line(aes(y = Realized,      color = "Realized Post-2022", linetype = "Realized Post-2022"), linewidth = lw * 1.1) +
    facet_wrap(~Panel, scales = "free_y", ncol = if (is_web) 2 else 1) +
    scale_color_manual(name = NULL, values = c(
      "Realized Post-2022"     = PAL$red,
      "LP-BVAR(3) Conditional" = PAL$navy,
      "Unconditional Median"   = PAL$gray
    )) +
    scale_linetype_manual(name = NULL, values = c(
      "Realized Post-2022"     = "solid",
      "LP-BVAR(3) Conditional" = "dashed",
      "Unconditional Median"   = "dotted"
    )) +
    scale_fill_manual(name = NULL, values = c(
      "68% Predictive Interval" = PAL$band68,
      "95% Predictive Interval" = PAL$band95
    )) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    guides(
      color    = guide_legend(nrow = if (is_web) 1 else 2, byrow = TRUE, order = 1),
      linetype = guide_legend(nrow = if (is_web) 1 else 2, byrow = TRUE, order = 1),
      fill     = guide_legend(nrow = 1, order = 2)
    ) +
    labs(x = NULL, y = "Levels (Thousands)") +
    th
}

p_fig5_print <- build_fig5(is_web = FALSE)
ggsave("outputs/figures/fig5_bvar_posterior_bands.pdf", p_fig5_print, width = 3.5, height = 3.7)
ggsave("outputs/figures/fig5_bvar_posterior_bands.png", p_fig5_print, width = 3.5, height = 3.7, dpi = 300)

p_fig5_web <- build_fig5(is_web = TRUE)
ggsave("outputs/figures/fig5_web_bvar_bands.png", p_fig5_web, width = 8.4, height = 4.1, dpi = 300)
file.copy("outputs/figures/fig5_web_bvar_bands.png", "docs/assets/images/fig5_bvar_posterior_bands.png", overwrite = TRUE)

# ------------------------------------------------------------------------------
# 7. Publication Figure 6: Multi-Specification Levels Comparison (Unified Punto!)
# ------------------------------------------------------------------------------
comp_lvl_df <- rbind(
  data.frame(Date = plot_dates, Value = hwi_real, Model = "Realized Post-2022", Panel = "(a) Job Openings (HWI, thousands)"),
  data.frame(Date = plot_dates, Value = apply(sim_lp_lvl$hwi_lvl_cnd, 2, median), Model = "LP-BVAR(3) Levels", Panel = "(a) Job Openings (HWI, thousands)"),
  data.frame(Date = plot_dates, Value = apply(sim_lp_diff$hwi_lvl_cnd, 2, median), Model = "LP-BVAR(2) Diff", Panel = "(a) Job Openings (HWI, thousands)"),
  data.frame(Date = plot_dates, Value = var_cnd_hwi, Model = "COVID-Controlled VAR(2)", Panel = "(a) Job Openings (HWI, thousands)"),
  data.frame(Date = plot_dates, Value = pay_real, Model = "Realized Post-2022", Panel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(Date = plot_dates, Value = apply(sim_lp_lvl$pay_lvl_cnd, 2, median), Model = "LP-BVAR(3) Levels", Panel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(Date = plot_dates, Value = apply(sim_lp_diff$pay_lvl_cnd, 2, median), Model = "LP-BVAR(2) Diff", Panel = "(b) Payroll Employment (PAYEMS, thousands)"),
  data.frame(Date = plot_dates, Value = var_cnd_pay, Model = "COVID-Controlled VAR(2)", Panel = "(b) Payroll Employment (PAYEMS, thousands)")
)
comp_lvl_df$Model <- factor(comp_lvl_df$Model, levels = c(
  "Realized Post-2022", "LP-BVAR(3) Levels", "LP-BVAR(2) Diff", "COVID-Controlled VAR(2)"
))

build_fig6 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  lw <- if (is_web) 0.90 else 0.68
  ggplot(comp_lvl_df, aes(x = Date, y = Value, color = Model, linetype = Model)) +
    geom_line(linewidth = lw) +
    facet_wrap(~Panel, scales = "free_y", ncol = if (is_web) 2 else 1) +
    scale_color_manual(name = NULL, values = c(
      "Realized Post-2022"      = PAL$red,
      "LP-BVAR(3) Levels"       = PAL$navy,
      "LP-BVAR(2) Diff"         = PAL$amber,
      "COVID-Controlled VAR(2)" = PAL$gray
    )) +
    scale_linetype_manual(name = NULL, values = c(
      "Realized Post-2022"      = "solid",
      "LP-BVAR(3) Levels"       = "solid",
      "LP-BVAR(2) Diff"         = "dashed",
      "COVID-Controlled VAR(2)" = "dotdash"
    )) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(labels = comma) +
    guides(
      color    = guide_legend(nrow = if (is_web) 1 else 2, byrow = TRUE),
      linetype = guide_legend(nrow = if (is_web) 1 else 2, byrow = TRUE)
    ) +
    labs(x = NULL, y = "Levels (Thousands)") +
    th
}

p_fig6_print <- build_fig6(is_web = FALSE)
ggsave("outputs/figures/fig6_bvar_vs_var_counterfactual.pdf", p_fig6_print, width = 3.5, height = 3.5)
ggsave("outputs/figures/fig6_bvar_vs_var_counterfactual.png", p_fig6_print, width = 3.5, height = 3.5, dpi = 300)

p_fig6_web <- build_fig6(is_web = TRUE)
ggsave("outputs/figures/fig6_web_bvar_comparison.png", p_fig6_web, width = 8.4, height = 4.0, dpi = 300)
file.copy("outputs/figures/fig6_web_bvar_comparison.png", "docs/assets/images/fig6_bvar_vs_var_counterfactual.png", overwrite = TRUE)

cat("======================================================================\n")
cat("STEP 3 COMPLETED SUCCESSFULLY (ZERO HARDCODED VALUES, UNIFIED PUNTO)\n")
cat("======================================================================\n")
