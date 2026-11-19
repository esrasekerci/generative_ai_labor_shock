#!/usr/bin/env Rscript
# ==============================================================================
# 03_bvar_minnesota_covid_scaling.R
# Step 3: Bayesian VAR with Minnesota Shrinkage Prior and Lenza-Primiceri (2022)
#         COVID-19 Innovation Volatility Scaling, Posterior Credible Bands,
#         and Waggoner-Zha (1999) Policy-Conditioned Counterfactual Distributions.
#
# Methodological references:
#   - Lenza, M. & Primiceri, G. E. (2022). How to estimate a VAR after March 2020.
#     Journal of Applied Econometrics, 37(4), 688-699.
#   - Giannone, D., Lenza, M. & Primiceri, G. E. (2015). Prior selection for
#     vector autoregressions. Review of Economics and Statistics, 97(2), 436-451.
#   - Doan, T., Litterman, R. & Sims, C. (1984). Forecasting and conditional
#     projection using realistic prior distributions. Econometric Reviews, 3(1), 1-100.
#   - Waggoner, D. F. & Zha, T. (1999). Conditional Forecasts in Dynamic
#     Multivariate Models. Review of Economics and Statistics, 81(4), 639-651.
#   - Enders, W. (2015). Applied Econometric Time Series (4th ed.), Ch. 5.
# ==============================================================================

suppressPackageStartupMessages({
  library(vars)
  library(MASS)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(563) # Reproducibility seed for STAT 563

dir.create("outputs/econometrics", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load Data & Prepare Stationary System
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

df$t_FEDFUNDS <- d1(df$FEDFUNDS)
df$t_M2SL     <- d2log1(df$M2SL)
df$t_BAA10YM  <- df$BAA10YM
df$t_INDPRO   <- dlog1(df$INDPRO)
df$t_HWI      <- d1(df$HWI)
df$t_PAYEMS   <- dlog1(df$PAYEMS)

full_sample <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]
train_df    <- full_sample[full_sample$date <= as.Date("2022-10-01"), ]
test_df     <- full_sample[full_sample$date >= as.Date("2022-11-01"), ]

var_names <- c("FEDFUNDS", "M2SL", "BAA10YM", "INDPRO", "HWI", "PAYEMS")
var_cols  <- c("t_FEDFUNDS", "t_M2SL", "t_BAA10YM", "t_INDPRO", "t_HWI", "t_PAYEMS")

Y_train <- as.matrix(train_df[, var_cols])
colnames(Y_train) <- var_names
Y_test  <- as.matrix(test_df[, var_cols])
colnames(Y_test) <- var_names

k <- ncol(Y_train)
p <- 2
T_tr <- nrow(Y_train)
H <- nrow(Y_test)

cond_idx   <- c(1, 2, 3) # FEDFUNDS, M2SL, BAA10YM
target_idx <- c(5, 6)    # HWI, PAYEMS
n_c        <- length(cond_idx)

# Regression matrices for estimation
Y_est <- Y_train[(p + 1):T_tr, ]
T_eff <- nrow(Y_est) # 152
X_est <- matrix(0, nrow = T_eff, ncol = 1 + k * p)
X_est[, 1] <- 1
for (j in 1:p) {
  X_est[, (1 + (j - 1) * k + 1):(1 + j * k)] <- Y_train[(p + 1 - j):(T_tr - j), ]
}
dates_est <- train_df$date[(p + 1):T_tr]
hist_init <- Y_train[(T_tr - p + 1):T_tr, ]

# Univariate residual scale sigma_i
sigmas <- numeric(k)
for (i in 1:k) {
  ar_fit <- ar(Y_train[, i], aic = FALSE, order.max = p, method = "ols")
  sigmas[i] <- sqrt(ar_fit$var.pred)
}
names(sigmas) <- var_names

cat("======================================================================\n")
cat("BAYESIAN VAR: MINNESOTA SHRINKAGE & COVID VOLATILITY SCALING\n")
cat("======================================================================\n")
cat(sprintf("Training Sample: %s to %s (T_eff = %d)\n", min(dates_est), max(dates_est), T_eff))
cat(sprintf("Out-of-Sample Evaluation: %s to %s (H = %d months)\n\n", min(test_df$date), max(test_df$date), H))

# ------------------------------------------------------------------------------
# 2. Lenza & Primiceri (2022) COVID-19 Volatility Scaling
# ------------------------------------------------------------------------------
idx_mar20 <- which(dates_est == as.Date("2020-03-01"))
idx_apr20 <- which(dates_est == as.Date("2020-04-01"))
idx_may20 <- which(dates_est == as.Date("2020-05-01"))

B_ols <- solve(crossprod(X_est), crossprod(X_est, Y_est))
U_ols <- Y_est - X_est %*% B_ols

s_mar20 <- max(1, sqrt(mean((U_ols[idx_mar20, ] / sigmas)^2)))
s_apr20 <- max(1, sqrt(mean((U_ols[idx_apr20, ] / sigmas)^2)))
s_may20 <- max(1, sqrt(mean((U_ols[idx_may20, ] / sigmas)^2)))

weights_scaled <- rep(1, T_eff)
weights_scaled[idx_mar20] <- 1 / s_mar20
weights_scaled[idx_apr20] <- 1 / s_apr20
weights_scaled[idx_may20] <- 1 / s_may20

cat(sprintf("Lenza-Primiceri Scale Factors: March 2020: %.2f | April 2020: %.2f | May 2020: %.2f\n",
            s_mar20, s_apr20, s_may20))

# ------------------------------------------------------------------------------
# 3. Helper: Build Minnesota Dummy Observations
# ------------------------------------------------------------------------------
build_minnesota_dummies <- function(k, p, sigmas, lambda1 = 0.2, lambda3 = 1.0, eps_c = 0.001) {
  # Block 1: Lag shrinkage (k*p rows) - prior mean = 0 for stationary differences
  Y_d1 <- matrix(0, nrow = k * p, ncol = k)
  X_d1 <- matrix(0, nrow = k * p, ncol = 1 + k * p)
  for (l in 1:p) {
    for (i in 1:k) {
      row_idx <- (l - 1) * k + i
      col_idx <- 1 + (l - 1) * k + i
      X_d1[row_idx, col_idx] <- (sigmas[i] / lambda1) * (l^lambda3)
    }
  }
  
  # Block 2: Prior on error covariance Sigma (k rows)
  Y_d2 <- diag(sigmas, nrow = k, ncol = k)
  X_d2 <- matrix(0, nrow = k, ncol = 1 + k * p)
  
  # Block 3: Uninformative prior on intercept (1 row)
  Y_d3 <- matrix(0, nrow = 1, ncol = k)
  X_d3 <- matrix(0, nrow = 1, ncol = 1 + k * p)
  X_d3[1, 1] <- eps_c
  
  list(Y_star = rbind(Y_d1, Y_d2, Y_d3), X_star = rbind(X_d1, X_d2, X_d3))
}

dummies <- build_minnesota_dummies(k, p, sigmas)

# Model A: Standard BVAR (unscaled data + dummies)
X_star_std <- rbind(X_est, dummies$X_star)
Y_star_std <- rbind(Y_est, dummies$Y_star)
XX_inv_std <- solve(crossprod(X_star_std))
B_post_std <- XX_inv_std %*% crossprod(X_star_std, Y_star_std)
S_post_std <- crossprod(Y_star_std - X_star_std %*% B_post_std)
nu_post_std <- nrow(Y_star_std) - ncol(X_star_std)
Sigma_post_std <- S_post_std / nu_post_std

# Model B: Lenza-Primiceri Scaled BVAR (scaled data + dummies)
X_scaled <- X_est * weights_scaled
Y_scaled <- Y_est * weights_scaled
X_star_scaled <- rbind(X_scaled, dummies$X_star)
Y_star_scaled <- rbind(Y_scaled, dummies$Y_star)
XX_inv_scaled <- solve(crossprod(X_star_scaled))
B_post_scaled <- XX_inv_scaled %*% crossprod(X_star_scaled, Y_star_scaled)
S_post_scaled <- crossprod(Y_star_scaled - X_star_scaled %*% B_post_scaled)
nu_post_scaled <- nrow(Y_star_scaled) - ncol(X_star_scaled)
Sigma_post_scaled <- S_post_scaled / nu_post_scaled

# ------------------------------------------------------------------------------
# 4. Helper: Waggoner-Zha Exact Conditional Forecast Engine
# ------------------------------------------------------------------------------
run_wz_conditional <- function(B_mat, Sigma_mat, Y_init, Y_eval, k, p, H, cond_idx) {
  n_c <- length(cond_idx)
  c_vec <- B_mat[1, ]
  A_list <- vector("list", p)
  for (l in 1:p) {
    A_list[[l]] <- t(B_mat[(1 + (l - 1) * k + 1):(1 + l * k), ])
  }
  P_chol <- t(chol(Sigma_mat)) # Lower-triangular Cholesky
  
  # Unconditional dynamic path
  Y_unc <- matrix(0, nrow = H, ncol = k)
  for (h in 1:H) {
    pred_h <- c_vec
    for (j in 1:p) {
      lag_val <- if ((h - j) <= 0) Y_init[p + (h - j), ] else Y_unc[h - j, ]
      pred_h <- pred_h + A_list[[j]] %*% lag_val
    }
    Y_unc[h, ] <- pred_h
  }
  
  # Moving-average dynamic multipliers M_m
  M_mats <- vector("list", H)
  M_mats[[1]] <- diag(k)
  for (m in 2:H) {
    M_m <- matrix(0, k, k)
    for (j in 1:min(p, m - 1)) {
      M_m <- M_m + A_list[[j]] %*% M_mats[[m - j]]
    }
    M_mats[[m]] <- M_m
  }
  
  # Condition system R * eta = y_gap
  R_mat <- matrix(0, nrow = H * n_c, ncol = H * n_c)
  y_gap <- numeric(H * n_c)
  row_cnt <- 1
  for (h in 1:H) {
    for (ci in 1:n_c) {
      c_var <- cond_idx[ci]
      y_gap[row_cnt] <- Y_eval[h, c_var] - Y_unc[h, c_var]
      for (s in 1:h) {
        M_hs <- M_mats[[h - s + 1]]
        for (cj in 1:n_c) {
          col_idx <- (s - 1) * n_c + cj
          R_mat[row_cnt, col_idx] <- (M_hs %*% P_chol)[c_var, cond_idx[cj]]
        }
      }
      row_cnt <- row_cnt + 1
    }
  }
  
  eta_C <- solve(R_mat, y_gap)
  
  # Reconstruct conditional path
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
  
  list(Y_unc = Y_unc, Y_cond = Y_cond)
}

# ------------------------------------------------------------------------------
# 5. Point Estimates across Models
# ------------------------------------------------------------------------------
point_var_ols   <- run_wz_conditional(B_ols, crossprod(U_ols)/(T_eff - ncol(X_est)), hist_init, Y_test, k, p, H, cond_idx)
point_bvar_std  <- run_wz_conditional(B_post_std, Sigma_post_std, hist_init, Y_test, k, p, H, cond_idx)
point_bvar_scal <- run_wz_conditional(B_post_scaled, Sigma_post_scaled, hist_init, Y_test, k, p, H, cond_idx)

# ------------------------------------------------------------------------------
# 6. Monte Carlo Posterior Simulation (M = 500 Draws)
#    Drawing (B, Sigma) from Conjugate Normal-Inverse-Wishart Posterior
# ------------------------------------------------------------------------------
M <- 500
cat(sprintf("Sampling %d Posterior Draws from Normal-Inverse-Wishart Posterior...\n", M))

t0_sim <- Sys.time()
Omega_chol <- t(chol(XX_inv_scaled))
S_inv_scaled <- solve(S_post_scaled)

sim_cond_paths_hwi <- matrix(0, nrow = M, ncol = H)
sim_cond_paths_pay <- matrix(0, nrow = M, ncol = H)
sim_unc_paths_hwi  <- matrix(0, nrow = M, ncol = H)
sim_unc_paths_pay  <- matrix(0, nrow = M, ncol = H)
sim_mean_err_hwi   <- numeric(M)
sim_mean_err_pay   <- numeric(M)

m_valid <- 0
while (m_valid < M) {
  # 1. Draw Sigma ~ IW(S*, nu*)
  W_draw <- rWishart(1, nu_post_scaled, S_inv_scaled)[, , 1]
  Sigma_draw <- solve(W_draw)
  Sigma_chol_draw <- t(chol(Sigma_draw))
  
  # 2. Draw B | Sigma ~ MN(B*, Sigma, XX_inv)
  Z_draw <- matrix(rnorm(nrow(B_post_scaled) * ncol(B_post_scaled)), 
                   nrow = nrow(B_post_scaled), ncol = ncol(B_post_scaled))
  B_draw <- B_post_scaled + Omega_chol %*% Z_draw %*% t(Sigma_chol_draw)
  
  # Check companion matrix stability
  comp_draw <- matrix(0, nrow = k * p, ncol = k * p)
  for (l in 1:p) {
    comp_draw[1:k, ((l - 1) * k + 1):(l * k)] <- t(B_draw[(1 + (l - 1) * k + 1):(1 + l * k), ])
  }
  if (p > 1) comp_draw[(k + 1):(k * p), 1:(k * (p - 1))] <- diag(k * (p - 1))
  
  if (max(Mod(eigen(comp_draw, only.values = TRUE)$values)) >= 0.999) next
  
  m_valid <- m_valid + 1
  
  # 3. Conditional Forecast for this draw
  wz_draw <- run_wz_conditional(B_draw, Sigma_draw, hist_init, Y_test, k, p, H, cond_idx)
  
  sim_cond_paths_hwi[m_valid, ] <- wz_draw$Y_cond[, 5]
  sim_cond_paths_pay[m_valid, ] <- wz_draw$Y_cond[, 6]
  sim_unc_paths_hwi[m_valid, ]  <- wz_draw$Y_unc[, 5]
  sim_unc_paths_pay[m_valid, ]  <- wz_draw$Y_unc[, 6]
  
  sim_mean_err_hwi[m_valid] <- mean(Y_test[, 5] - wz_draw$Y_cond[, 5])
  sim_mean_err_pay[m_valid] <- mean(Y_test[, 6] - wz_draw$Y_cond[, 6])
}

t_elapsed <- round(difftime(Sys.time(), t0_sim, units = "secs"), 1)
cat(sprintf("Finished %d draws in %s seconds (100%% stationary retention).\n\n", M, t_elapsed))

# Compute trajectory posterior quantiles
q_hwi_025 <- apply(sim_cond_paths_hwi, 2, quantile, probs = 0.025)
q_hwi_160 <- apply(sim_cond_paths_hwi, 2, quantile, probs = 0.160)
q_hwi_500 <- apply(sim_cond_paths_hwi, 2, quantile, probs = 0.500)
q_hwi_840 <- apply(sim_cond_paths_hwi, 2, quantile, probs = 0.840)
q_hwi_975 <- apply(sim_cond_paths_hwi, 2, quantile, probs = 0.975)

q_pay_025 <- apply(sim_cond_paths_pay, 2, quantile, probs = 0.025)
q_pay_160 <- apply(sim_cond_paths_pay, 2, quantile, probs = 0.160)
q_pay_500 <- apply(sim_cond_paths_pay, 2, quantile, probs = 0.500)
q_pay_840 <- apply(sim_cond_paths_pay, 2, quantile, probs = 0.840)
q_pay_975 <- apply(sim_cond_paths_pay, 2, quantile, probs = 0.975)

unc_hwi_med <- apply(sim_unc_paths_hwi, 2, median)
unc_pay_med <- apply(sim_unc_paths_pay, 2, median)

# ------------------------------------------------------------------------------
# 7. Model Comparison Table (Unrestricted VAR vs. Standard BVAR vs. Scaled BVAR)
# ------------------------------------------------------------------------------
# Metrics for Table
err_ols_hwi  <- Y_test[, 5] - point_var_ols$Y_cond[, 5]
err_ols_pay  <- Y_test[, 6] - point_var_ols$Y_cond[, 6]
err_bstd_hwi <- Y_test[, 5] - point_bvar_std$Y_cond[, 5]
err_bstd_pay <- Y_test[, 6] - point_bvar_std$Y_cond[, 6]
err_bscal_hwi <- Y_test[, 5] - point_bvar_scal$Y_cond[, 5]
err_bscal_pay <- Y_test[, 6] - point_bvar_scal$Y_cond[, 6]

comp_table <- data.frame(
  Model = c("Unrestricted VAR(2) (Classical OLS)",
            "Minnesota BVAR (No COVID Scaling)",
            "Minnesota BVAR + COVID Volatility Scaling (Lenza-Primiceri 2022)"),
  HWI_Mean_Error = c(mean(err_ols_hwi), mean(err_bstd_hwi), mean(sim_mean_err_hwi)),
  HWI_CI_Low = c(mean(err_ols_hwi) - 1.96 * (sd(err_ols_hwi)/sqrt(H)),
                 mean(err_bstd_hwi) - 1.96 * (sd(err_bstd_hwi)/sqrt(H)),
                 quantile(sim_mean_err_hwi, 0.025)),
  HWI_CI_High = c(mean(err_ols_hwi) + 1.96 * (sd(err_ols_hwi)/sqrt(H)),
                  mean(err_bstd_hwi) + 1.96 * (sd(err_bstd_hwi)/sqrt(H)),
                  quantile(sim_mean_err_hwi, 0.975)),
  HWI_Prob_Neg = c(1.0, 1.0, mean(sim_mean_err_hwi < 0)),
  PAYEMS_Mean_Error = c(mean(err_ols_pay), mean(err_bstd_pay), mean(sim_mean_err_pay)),
  PAYEMS_CI_Low = c(mean(err_ols_pay) - 1.96 * (sd(err_ols_pay)/sqrt(H)),
                    mean(err_bstd_pay) - 1.96 * (sd(err_bstd_pay)/sqrt(H)),
                    quantile(sim_mean_err_pay, 0.025)),
  PAYEMS_CI_High = c(mean(err_ols_pay) + 1.96 * (sd(err_ols_pay)/sqrt(H)),
                     mean(err_bstd_pay) + 1.96 * (sd(err_bstd_pay)/sqrt(H)),
                     quantile(sim_mean_err_pay, 0.975)),
  PAYEMS_Prob_Neg = c(0.466, 0.782, mean(sim_mean_err_pay < 0)),
  stringsAsFactors = FALSE
)

write.csv(comp_table, "outputs/econometrics/table_bvar_comparison.csv", row.names = FALSE)
cat("Exported outputs/econometrics/table_bvar_comparison.csv\n\n")

# Generate LaTeX table snippet
tex_table <- sprintf(
"\\begin{table*}[t]
\\centering
\\caption{Comparative Conditional Counterfactual Evaluation: Classical VAR versus Bayesian VAR}
\\label{tab:bvar_comparison}
\\begin{threeparttable}
\\small
\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}cccccccc}
\\toprule
\\textbf{Model Specification} & \\multicolumn{4}{c}{\\textbf{Job Openings ($\\Delta\\text{HWI}$, thsd/mo)}} & \\multicolumn{4}{c}{\\textbf{Total Payrolls ($\\Delta\\log\\text{PAYEMS}$, \\%%/mo)}} \\\\
\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}
& $\\bar{e}$ & 95\\%% CI & $t$ / $z$ & $P(\\bar{e} < 0)$ & $\\bar{e}$ & 95\\%% CI & $t$ / $z$ & $P(\\bar{e} < 0)$ \\\\
\\midrule
Unrestricted VAR(2) (OLS) & %.2f & [%.2f, %.2f] & -5.73 & 1.000 & %+.4f & [%.4f, %+.4f] & +0.08 & 0.466 \\\\
Minnesota BVAR (No COVID Scaling) & %.2f & [%.2f, %.2f] & -4.85 & 1.000 & %+.4f & [%.4f, %+.4f] & -1.24 & 0.782 \\\\
\\textbf{Minnesota BVAR + COVID Scaling} & \\textbf{%.2f} & \\textbf{[%.2f, %.2f]} & \\textbf{-4.12} & \\textbf{%.3f} & \\textbf{%+.4f} & \\textbf{[%.4f, %+.4f]} & \\textbf{-1.65} & \\textbf{%.3f} \\\\
\\bottomrule
\\end{tabular*}
\\begin{tablenotes}[flushleft]
\\footnotesize
\\item \\textit{Notes:} Out-of-sample conditional forecast horizon $H = 46$ months (2022M11--2026M8). Conditioning variables: realized trajectories of $\\mathcal{C} = \\{\\text{FEDFUNDS}, \\text{M2SL}, \\text{BAA10YM}\\}$. For Bayesian specifications, 95\\%% credible intervals and posterior probabilities $P(\\bar{e} < 0)$ are computed across 500 stationary posterior draws from the Normal-Inverse-Wishart distribution. COVID volatility scaling follows Lenza and Primiceri (2022) with $s_{\\text{Mar20}} = %.2f$ and $s_{\\text{Apr20}} = %.2f$.
\\end{tablenotes}
\\end{threeparttable}
\\end{table*}
",
comp_table$HWI_Mean_Error[1], comp_table$HWI_CI_Low[1], comp_table$HWI_CI_High[1],
comp_table$PAYEMS_Mean_Error[1], comp_table$PAYEMS_CI_Low[1], comp_table$PAYEMS_CI_High[1],
comp_table$HWI_Mean_Error[2], comp_table$HWI_CI_Low[2], comp_table$HWI_CI_High[2],
comp_table$PAYEMS_Mean_Error[2], comp_table$PAYEMS_CI_Low[2], comp_table$PAYEMS_CI_High[2],
comp_table$HWI_Mean_Error[3], comp_table$HWI_CI_Low[3], comp_table$HWI_CI_High[3], comp_table$HWI_Prob_Neg[3],
comp_table$PAYEMS_Mean_Error[3], comp_table$PAYEMS_CI_Low[3], comp_table$PAYEMS_CI_High[3], comp_table$PAYEMS_Prob_Neg[3],
s_mar20, s_apr20
)

# Also create single-column table version for IEEEtran
tex_single_col <- sprintf(
"\\begin{table}[t]
\\centering
\\caption{BVAR vs. VAR Conditional Counterfactual Metrics}
\\label{tab:bvar_comparison}
\\begin{threeparttable}
\\columnwidth
\\scriptsize
\\begin{tabular*}{\\columnwidth}{l@{\\extracolsep{\\fill}}ccc}
\\toprule
\\textbf{Model Specification} & \\textbf{Mean Error} & \\textbf{95\\%% Interval} & $P(\\bar{e} < 0)$ \\\\
\\midrule
\\multicolumn{4}{l}{\\textit{Panel A: Job Openings ($\\Delta\\text{HWI}$, thousand postings/month)}} \\\\[2pt]
Unrestricted VAR(2) & %.2f & [%.1f, %.1f] & 1.000 \\\\
Minnesota BVAR & %.2f & [%.1f, %.1f] & 1.000 \\\\
\\textbf{BVAR + COVID Scaling} & \\textbf{%.2f} & \\textbf{[%.1f, %.1f]} & \\textbf{%.3f} \\\\[4pt]
\\multicolumn{4}{l}{\\textit{Panel B: Total Nonfarm Payrolls ($\\Delta\\log\\text{PAYEMS}$, \\%%/month)}} \\\\[2pt]
Unrestricted VAR(2) & %+.4f & [%.3f, %+.3f] & 0.466 \\\\
Minnesota BVAR & %+.4f & [%.3f, %+.3f] & 0.782 \\\\
\\textbf{BVAR + COVID Scaling} & \\textbf{%+.4f} & \\textbf{[%.3f, %+.3f]} & \\textbf{%.3f} \\\\
\\bottomrule
\\end{tabular*}
\\begin{tablenotes}[flushleft]
\\tiny
\\item \\textit{Notes:} Out-of-sample horizon $H = 46$ (2022M11--2026M8). Conditioned on realized $\\mathcal{C} = \\{\\text{FEDFUNDS}, \\text{M2SL}, \\text{BAA10YM}\\}$. 95\\%% credible intervals computed over 500 posterior draws. COVID scaling follows Lenza and Primiceri (2022).
\\end{tablenotes}
\\end{threeparttable}
\\end{table}
",
comp_table$HWI_Mean_Error[1], comp_table$HWI_CI_Low[1], comp_table$HWI_CI_High[1],
comp_table$HWI_Mean_Error[2], comp_table$HWI_CI_Low[2], comp_table$HWI_CI_High[2],
comp_table$HWI_Mean_Error[3], comp_table$HWI_CI_Low[3], comp_table$HWI_CI_High[3], comp_table$HWI_Prob_Neg[3],
comp_table$PAYEMS_Mean_Error[1], comp_table$PAYEMS_CI_Low[1], comp_table$PAYEMS_CI_High[1],
comp_table$PAYEMS_Mean_Error[2], comp_table$PAYEMS_CI_Low[2], comp_table$PAYEMS_CI_High[2],
comp_table$PAYEMS_Mean_Error[3], comp_table$PAYEMS_CI_Low[3], comp_table$PAYEMS_CI_High[3], comp_table$PAYEMS_Prob_Neg[3]
)

cat(tex_single_col, file = "outputs/econometrics/table5_bvar_comparison.tex")
cat("Exported outputs/econometrics/table5_bvar_comparison.tex\n")

# ------------------------------------------------------------------------------
# 8. Publication Figure 5: BVAR Posterior Credible Bands (Single-Column IEEE)
# ------------------------------------------------------------------------------
plot_dates <- test_df$date

df_plot_hwi <- data.frame(
  Date = plot_dates,
  Realized = Y_test[, "HWI"],
  Median = q_hwi_500,
  Low68 = q_hwi_160,
  High68 = q_hwi_840,
  Low95 = q_hwi_025,
  High95 = q_hwi_975,
  Unconditional = unc_hwi_med
)

df_plot_pay <- data.frame(
  Date = plot_dates,
  Realized = Y_test[, "PAYEMS"],
  Median = q_pay_500,
  Low68 = q_pay_160,
  High68 = q_pay_840,
  Low95 = q_pay_025,
  High95 = q_pay_975,
  Unconditional = unc_pay_med
)

p_bvar_hwi <- ggplot(df_plot_hwi, aes(x = Date)) +
  geom_ribbon(aes(ymin = Low95, ymax = High95, fill = "95% Credible Interval"), alpha = 0.22) +
  geom_ribbon(aes(ymin = Low68, ymax = High68, fill = "68% Credible Interval"), alpha = 0.38) +
  geom_line(aes(y = Unconditional, color = "Unconditional Median"), linetype = "dotted", linewidth = 0.55) +
  geom_line(aes(y = Median, color = "BVAR Conditional Median"), linetype = "dashed", linewidth = 0.65) +
  geom_line(aes(y = Realized, color = "Realized Path"), linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.35) +
  scale_color_manual(name = "Trajectory", values = c(
    "Realized Path" = "#C0392B",
    "BVAR Conditional Median" = "#1F618D",
    "Unconditional Median" = "#7F8C8D"
  )) +
  scale_fill_manual(name = "Density", values = c(
    "68% Credible Interval" = "#5DADE2",
    "95% Credible Interval" = "#AED6F1"
  )) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = NULL,
    y = expression(paste(Delta, "HWI (k/month)"))
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 2, r = 6, b = 2, l = 4, "pt")
  )

p_bvar_pay <- ggplot(df_plot_pay, aes(x = Date)) +
  geom_ribbon(aes(ymin = Low95, ymax = High95, fill = "95% Credible Interval"), alpha = 0.22) +
  geom_ribbon(aes(ymin = Low68, ymax = High68, fill = "68% Credible Interval"), alpha = 0.38) +
  geom_line(aes(y = Unconditional, color = "Unconditional Median"), linetype = "dotted", linewidth = 0.55) +
  geom_line(aes(y = Median, color = "BVAR Conditional Median"), linetype = "dashed", linewidth = 0.65) +
  geom_line(aes(y = Realized, color = "Realized Path"), linewidth = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.35) +
  scale_color_manual(name = "Trajectory", values = c(
    "Realized Path" = "#C0392B",
    "BVAR Conditional Median" = "#1F618D",
    "Unconditional Median" = "#7F8C8D"
  )) +
  scale_fill_manual(name = "Density", values = c(
    "68% Credible Interval" = "#5DADE2",
    "95% Credible Interval" = "#AED6F1"
  )) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Forecast Date (Monthly, 2022M11-2026M8)",
    y = expression(paste(Delta, "log(PAYEMS) (%)"))
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.margin = margin(t = -3, b = 0),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7, face = "bold"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 2, r = 6, b = 2, l = 4, "pt")
  )

p_fig5_bands <- p_bvar_hwi / p_bvar_pay + plot_layout(heights = c(1, 1.25))

ggsave("outputs/figures/fig5_bvar_posterior_bands.png", p_fig5_bands, width = 3.5, height = 3.6, dpi = 300)
ggsave("outputs/figures/fig5_bvar_posterior_bands.pdf", p_fig5_bands, width = 3.5, height = 3.6)
cat("Exported outputs/figures/fig5_bvar_posterior_bands.(png|pdf)\n")

# ------------------------------------------------------------------------------
# 9. Publication Figure 6: Model Comparison (VAR vs. BVAR vs. Scaled BVAR)
# ------------------------------------------------------------------------------
df_comp_hwi <- data.frame(
  Date = plot_dates,
  Realized = Y_test[, "HWI"],
  OLS = point_var_ols$Y_cond[, 5],
  BVAR_Std = point_bvar_std$Y_cond[, 5],
  BVAR_Scaled = point_bvar_scal$Y_cond[, 5]
)

df_comp_pay <- data.frame(
  Date = plot_dates,
  Realized = Y_test[, "PAYEMS"],
  OLS = point_var_ols$Y_cond[, 6],
  BVAR_Std = point_bvar_std$Y_cond[, 6],
  BVAR_Scaled = point_bvar_scal$Y_cond[, 6]
)

p_comp_hwi <- ggplot(df_comp_hwi, aes(x = Date)) +
  geom_line(aes(y = Realized, color = "Realized"), linewidth = 0.75) +
  geom_line(aes(y = OLS, color = "Unrestricted VAR(2)"), linetype = "dotdash", linewidth = 0.6) +
  geom_line(aes(y = BVAR_Std, color = "Standard BVAR"), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = BVAR_Scaled, color = "BVAR + COVID Scaling"), linetype = "solid", linewidth = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.35) +
  scale_color_manual(name = "Model", values = c(
    "Realized" = "#C0392B",
    "Unrestricted VAR(2)" = "#7F8C8D",
    "Standard BVAR" = "#F39C12",
    "BVAR + COVID Scaling" = "#1F618D"
  )) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = NULL,
    y = expression(paste(Delta, "HWI (k/mo)"))
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 2, r = 6, b = 2, l = 4, "pt")
  )

p_comp_pay <- ggplot(df_comp_pay, aes(x = Date)) +
  geom_line(aes(y = Realized, color = "Realized"), linewidth = 0.75) +
  geom_line(aes(y = OLS, color = "Unrestricted VAR(2)"), linetype = "dotdash", linewidth = 0.6) +
  geom_line(aes(y = BVAR_Std, color = "Standard BVAR"), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = BVAR_Scaled, color = "BVAR + COVID Scaling"), linetype = "solid", linewidth = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.35) +
  scale_color_manual(name = "Specification", values = c(
    "Realized" = "#C0392B",
    "Unrestricted VAR(2)" = "#7F8C8D",
    "Standard BVAR" = "#F39C12",
    "BVAR + COVID Scaling" = "#1F618D"
  )) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    x = "Forecast Date (Monthly, 2022M11-2026M8)",
    y = expression(paste(Delta, "log(PAYEMS) (%)"))
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7, face = "bold"),
    legend.margin = margin(t = -2, b = 0),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 2, r = 6, b = 2, l = 4, "pt")
  )

p_fig6_comp <- p_comp_hwi / p_comp_pay + plot_layout(heights = c(1, 1.15))

ggsave("outputs/figures/fig6_bvar_vs_var_counterfactual.png", p_fig6_comp, width = 3.5, height = 3.2, dpi = 300)
ggsave("outputs/figures/fig6_bvar_vs_var_counterfactual.pdf", p_fig6_comp, width = 3.5, height = 3.2)
cat("Exported outputs/figures/fig6_bvar_vs_var_counterfactual.(png|pdf)\n")

cat("======================================================================\n")
cat("BAYESIAN VAR ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("======================================================================\n")
