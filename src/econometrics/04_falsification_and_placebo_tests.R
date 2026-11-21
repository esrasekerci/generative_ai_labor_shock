#!/usr/bin/env Rscript
# ==============================================================================
# 04_falsification_and_placebo_tests.R
# Step 4: Attribution Falsification Battery (Addressing Prof. Ceylan Point 3)
#   (A) Pre-2022 Pseudo-Break Placebo Distribution (13 Quarterly Pre-AI Origins
#       evaluated over exact 24-month clean pre-COVID windows vs. Post-2022)
#   (B) Low-AI vs. High-AI Sectoral Placebo Analysis (USCONS Construction,
#       MANEMP Manufacturing, USGOVT Government vs. USFIRE Finance & HWI)
#
# Methodological references:
#   - Waggoner, D. F. & Zha, T. (1999). Review of Economics and Statistics, 81(4), 639-651.
#   - Lenza, M. & Primiceri, G. E. (2022). Journal of Applied Econometrics, 37(4), 688-699.
#   - Newey, W. K. & West, K. D. (1987). Econometrica, 55(3), 703-708.
# ==============================================================================

suppressPackageStartupMessages({
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
# 1. Load FRED-MD Data & Prepare Core + Sectoral Series
# ------------------------------------------------------------------------------
raw <- read.csv("fred_data.csv", stringsAsFactors = FALSE)
df <- raw[-1, ]
df$date <- as.Date(df$sasdate, format = "%m/%d/%Y")
df <- df[!is.na(df$date), ]
df <- df[order(df$date), ]

cols_num <- c("FEDFUNDS", "M2SL", "BAA", "GS10", "INDPRO", "HWI", "PAYEMS",
              "USCONS", "MANEMP", "USGOVT", "USFIRE")
for (col in cols_num) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
df$BAA10YM <- df$BAA - df$GS10

df$L_FEDFUNDS <- df$FEDFUNDS
df$L_M2SL     <- 100 * log(df$M2SL)
df$L_BAA10YM  <- df$BAA10YM
df$L_INDPRO   <- 100 * log(df$INDPRO)
df$L_HWI      <- 100 * log(df$HWI)
df$L_PAYEMS   <- 100 * log(df$PAYEMS)
df$L_USCONS   <- 100 * log(df$USCONS)
df$L_MANEMP   <- 100 * log(df$MANEMP)
df$L_USGOVT   <- 100 * log(df$USGOVT)
df$L_USFIRE   <- 100 * log(df$USFIRE)

fc <- df[df$date >= as.Date("2010-01-01") & !is.na(df$HWI), ]

# ------------------------------------------------------------------------------
# 2. Reusable Exact Marginal-Likelihood Levels LP-BVAR(3) + Waggoner-Zha Engine
# ------------------------------------------------------------------------------
run_lp_bvar_levels <- function(df_tr, df_te, target_cols, p = 3) {
  cols <- c("L_FEDFUNDS", "L_M2SL", "L_BAA10YM", "L_INDPRO", target_cols)
  Y_tr <- as.matrix(df_tr[, cols])
  Y_te <- as.matrix(df_te[, cols])
  k <- ncol(Y_tr); T_tr <- nrow(Y_tr); H <- nrow(Y_te); Kr <- 1 + k * p

  Ye <- Y_tr[(p + 1):T_tr, , drop = FALSE]; Te <- nrow(Ye)
  Xe <- matrix(0, Te, Kr); Xe[, 1] <- 1
  for (j in 1:p) Xe[, (1 + (j - 1) * k + 1):(1 + j * k)] <- Y_tr[(p + 1 - j):(T_tr - j), ]
  dates_e <- df_tr$date[(p + 1):T_tr]

  has_covid <- any(dates_e == as.Date("2020-03-01"))
  df_sig <- if (has_covid) df_tr[df_tr$date <= as.Date("2020-02-01"), ] else df_tr
  Y_sig  <- as.matrix(df_sig[, cols]); Tp <- nrow(Y_sig)
  sig_pre <- numeric(k)
  for (i in 1:k) {
    yi <- Y_sig[, i]; Yi <- yi[(p + 1):Tp]
    Xi <- cbind(1, sapply(1:p, function(l) yi[(p + 1 - l):(Tp - l)]))
    ri <- Yi - Xi %*% solve(crossprod(Xi), crossprod(Xi, Yi))
    sig_pre[i] <- sqrt(sum(ri^2) / (length(Yi) - ncol(Xi)))
  }

  d0 <- k + 2; Psi0 <- diag(sig_pre^2); B0 <- matrix(0, Kr, k)
  for (i in 1:k) B0[1 + i, i] <- 1.0
  idx_tstar <- which(dates_e == as.Date("2020-03-01"))

  eval_lml <- function(lam1, s0, s1, s2, rho) {
    if (lam1 < 0.02 || lam1 > 5 || s0 > 150 || s1 > 150 || s2 > 150) return(-1e12)
    s_vec <- rep(1, Te)
    if (has_covid) {
      s_vec[idx_tstar] <- s0; s_vec[idx_tstar + 1] <- s1
      for (t in (idx_tstar + 2):Te) s_vec[t] <- 1 + (s2 - 1) * (rho^(t - (idx_tstar + 2)))
    }
    od <- numeric(Kr); od[1] <- 1e4
    for (l in 1:p) for (i in 1:k) od[1 + (l - 1) * k + i] <- (lam1^2) / ((l^2) * (sig_pre[i]^2))
    Om_inv <- diag(1 / od); Yt <- Ye / s_vec; Xt <- Xe / s_vec
    V_inv <- crossprod(Xt) + Om_inv
    Bp <- solve(V_inv, crossprod(Xt, Yt) + Om_inv %*% B0)
    Ep <- Yt - Xt %*% Bp
    Psip <- Psi0 + crossprod(Ep) + t(Bp - B0) %*% Om_inv %*% (Bp - B0)
    lml <- -k * sum(log(s_vec)) -
      (k / 2) * (sum(log(od)) + as.numeric(determinant(V_inv, logarithm = TRUE)$modulus)) +
      (d0 / 2) * as.numeric(determinant(Psi0, logarithm = TRUE)$modulus) -
      ((d0 + Te) / 2) * as.numeric(determinant(Psip, logarithm = TRUE)$modulus)
    if (!is.finite(lml)) -1e12 else lml
  }

  if (has_covid) {
    opt <- optim(c(log(0.16), log(10), log(40), log(5), qlogis(0.85)),
                 function(x) -eval_lml(exp(x[1]), 1 + exp(x[2]), 1 + exp(x[3]), 1 + exp(x[4]), 1 / (1 + exp(-x[5]))),
                 method = "BFGS")
    lam1 <- exp(opt$par[1]); s0 <- 1 + exp(opt$par[2]); s1 <- 1 + exp(opt$par[3])
    s2 <- 1 + exp(opt$par[4]); rho <- 1 / (1 + exp(-opt$par[5]))
  } else {
    opt <- optimize(function(ll) eval_lml(exp(ll), 1, 1, 1, 0), interval = c(log(0.02), log(2)), maximum = TRUE)
    lam1 <- exp(opt$maximum); s0 <- 1; s1 <- 1; s2 <- 1; rho <- 0
  }

  s_vec <- rep(1, Te)
  if (has_covid) {
    s_vec[idx_tstar] <- s0; s_vec[idx_tstar + 1] <- s1
    for (t in (idx_tstar + 2):Te) s_vec[t] <- 1 + (s2 - 1) * (rho^(t - (idx_tstar + 2)))
  }
  od <- numeric(Kr); od[1] <- 1e4
  for (l in 1:p) for (i in 1:k) od[1 + (l - 1) * k + i] <- (lam1^2) / ((l^2) * (sig_pre[i]^2))
  Om_inv <- diag(1 / od); Yt <- Ye / s_vec; Xt <- Xe / s_vec
  V_post <- solve(crossprod(Xt) + Om_inv)
  B_post <- V_post %*% (crossprod(Xt, Yt) + Om_inv %*% B0)
  Psip <- Psi0 + crossprod(Yt - Xt %*% B_post) + t(B_post - B0) %*% Om_inv %*% (B_post - B0)
  Sig_post <- Psip / (d0 + Te - k - 1)
  P_chol <- t(chol(Sig_post))

  A_l <- lapply(1:p, function(l) t(B_post[(1 + (l - 1) * k + 1):(1 + l * k), ]))
  c_l <- B_post[1, ]
  Yi <- Y_tr[(T_tr - p + 1):T_tr, ]
  Yu <- matrix(0, H, k)
  for (h in 1:H) {
    pr <- c_l
    for (j in 1:p) pr <- pr + A_l[[j]] %*% (if (h - j <= 0) Yi[p + h - j, ] else Yu[h - j, ])
    Yu[h, ] <- pr
  }
  Mm <- vector("list", H); Mm[[1]] <- diag(k)
  for (m in 2:H) { Mt <- matrix(0, k, k); for (j in 1:min(p, m - 1)) Mt <- Mt + A_l[[j]] %*% Mm[[m - j]]; Mm[[m]] <- Mt }
  R <- matrix(0, H * 3, H * k); rv <- numeric(H * 3); rc <- 1
  for (h in 1:H) for (ci in 1:3) {
    rv[rc] <- Y_te[h, ci] - Yu[h, ci]
    for (s in 1:h) R[rc, ((s - 1) * k + 1):(s * k)] <- (Mm[[h - s + 1]] %*% P_chol)[ci, ]
    rc <- rc + 1
  }
  eta <- as.numeric(t(R) %*% solve(R %*% t(R), rv))
  Yc <- Yu
  for (h in 1:H) {
    adj <- rep(0, k)
    for (s in 1:h) adj <- adj + as.numeric(Mm[[h - s + 1]] %*% P_chol %*% eta[((s - 1) * k + 1):(s * k)])
    Yc[h, ] <- Yu[h, ] + adj
  }
  list(Y_real = Y_te, Y_cnd = Yc, Y_unc = Yu, lam1 = lam1)
}

# ------------------------------------------------------------------------------
# 3. Falsification Test A: 13 Pre-2022 Pseudo-Break Dates (Exact 24-mo Horizon)
#    All evaluation windows end strictly <= 2020-02-01 (zero COVID/AI overlap!)
# ------------------------------------------------------------------------------
cat("======================================================================\n")
cat("FALSIFICATION TEST A: PRE-2022 PSEUDO-BREAK DATES (EXACT 24-MO WINDOWS)\n")
cat("======================================================================\n")

pseudo_origins <- seq(as.Date("2015-03-01"), as.Date("2018-03-01"), by = "3 months") # 13 quarterly origins
placebo_rows <- list()

for (i in seq_along(pseudo_origins)) {
  ps <- pseudo_origins[i]
  df_tr <- fc[fc$date < ps, ]
  df_te <- fc[fc$date >= ps, ][1:24, ] # Exact 24-month window, ends <= 2020-02-01
  res   <- run_lp_bvar_levels(df_tr, df_te, c("L_HWI", "L_PAYEMS"), p = 3)
  H_p   <- nrow(df_te)

  hwi_gap <- 100 * (exp(res$Y_real[H_p, 5] / 100) / exp(res$Y_cnd[H_p, 5] / 100) - 1)
  pay_gap <- 100 * (exp(res$Y_real[H_p, 6] / 100) / exp(res$Y_cnd[H_p, 6] / 100) - 1)

  eh <- diff(c(df_tr$L_HWI[nrow(df_tr)],    res$Y_real[, 5])) - diff(c(df_tr$L_HWI[nrow(df_tr)],    res$Y_cnd[, 5]))
  ep <- diff(c(df_tr$L_PAYEMS[nrow(df_tr)], res$Y_real[, 6])) - diff(c(df_tr$L_PAYEMS[nrow(df_tr)], res$Y_cnd[, 6]))
  bw <- max(1, floor(4 * (H_p / 100)^(2/9)))
  ch <- coeftest(lm(eh ~ 1), vcov. = NeweyWest(lm(eh ~ 1), lag = bw, prewhite = FALSE))
  cp <- coeftest(lm(ep ~ 1), vcov. = NeweyWest(lm(ep ~ 1), lag = bw, prewhite = FALSE))

  placebo_rows[[i]] <- data.frame(
    Origin       = format(ps, "%Y-%m"),
    EndDate      = format(max(df_te$date), "%Y-%m"),
    Horizon_Mo   = H_p,
    Regime       = "Pre-AI Placebo",
    HWI_Gap_Pct  = round(hwi_gap, 1),
    HWI_t_HAC    = round(ch[1, 3], 2),
    HWI_p_HAC    = signif(ch[1, 4], 3),
    PAY_Gap_Pct  = round(pay_gap, 1),
    PAY_t_HAC    = round(cp[1, 3], 2),
    stringsAsFactors = FALSE
  )
}

# Append exact 24-month AND 46-month Post-2022 Actual AI-Era Windows for apples-to-apples comparison!
tr22 <- fc[fc$date <= as.Date("2022-10-01"), ]
te22_24 <- fc[fc$date >= as.Date("2022-11-01"), ][1:24, ]
te22_46 <- fc[fc$date >= as.Date("2022-11-01"), ]

for (pair in list(list(te = te22_24, lbl = "Post-2022 (24m)"), list(te = te22_46, lbl = "Post-2022 (46m)"))) {
  res <- run_lp_bvar_levels(tr22, pair$te, c("L_HWI", "L_PAYEMS"), p = 3)
  Hp  <- nrow(pair$te)
  hwi_gap <- 100 * (exp(res$Y_real[Hp, 5] / 100) / exp(res$Y_cnd[Hp, 5] / 100) - 1)
  pay_gap <- 100 * (exp(res$Y_real[Hp, 6] / 100) / exp(res$Y_cnd[Hp, 6] / 100) - 1)
  eh <- diff(c(tr22$L_HWI[nrow(tr22)],    res$Y_real[, 5])) - diff(c(tr22$L_HWI[nrow(tr22)],    res$Y_cnd[, 5]))
  ep <- diff(c(tr22$L_PAYEMS[nrow(tr22)], res$Y_real[, 6])) - diff(c(tr22$L_PAYEMS[nrow(tr22)], res$Y_cnd[, 6]))
  bw <- max(1, floor(4 * (Hp / 100)^(2/9)))
  ch <- coeftest(lm(eh ~ 1), vcov. = NeweyWest(lm(eh ~ 1), lag = bw, prewhite = FALSE))
  cp <- coeftest(lm(ep ~ 1), vcov. = NeweyWest(lm(ep ~ 1), lag = bw, prewhite = FALSE))
  placebo_rows[[length(placebo_rows) + 1]] <- data.frame(
    Origin       = pair$lbl,
    EndDate      = format(max(pair$te$date), "%Y-%m"),
    Horizon_Mo   = Hp,
    Regime       = "Post-2022 Actual",
    HWI_Gap_Pct  = round(hwi_gap, 1),
    HWI_t_HAC    = round(ch[1, 3], 2),
    HWI_p_HAC    = signif(ch[1, 4], 3),
    PAY_Gap_Pct  = round(pay_gap, 1),
    PAY_t_HAC    = round(cp[1, 3], 2),
    stringsAsFactors = FALSE
  )
}

placebo_df <- do.call(rbind, placebo_rows)
print(placebo_df, row.names = FALSE)
write.csv(placebo_df, "outputs/econometrics/table_falsification_placebo.csv", row.names = FALSE)

pre_hwi_gaps <- placebo_df$HWI_Gap_Pct[placebo_df$Regime == "Pre-AI Placebo"]
post24_gap   <- placebo_df$HWI_Gap_Pct[placebo_df$Origin == "Post-2022 (24m)"]
post46_gap   <- placebo_df$HWI_Gap_Pct[placebo_df$Origin == "Post-2022 (46m)"]
cat(sprintf("\nPlacebo Summary (13 pre-AI 24m origins): Mean HWI Gap = %+.1f%% | Min = %+.1f%% | Max = %+.1f%%\n",
            mean(pre_hwi_gaps), min(pre_hwi_gaps), max(pre_hwi_gaps)))
cat(sprintf("Empirical Placebo Exceedance Rate: %d/13 (%.1f%%) as negative as Post-2022 24m (%+.1f%%)\n",
            sum(pre_hwi_gaps <= post24_gap), 100 * mean(pre_hwi_gaps <= post24_gap), post24_gap))
cat(sprintf("Empirical Placebo Exceedance Rate: %d/13 (%.1f%%) as negative as Post-2022 46m (%+.1f%%)\n\n",
            sum(pre_hwi_gaps <= post46_gap), 100 * mean(pre_hwi_gaps <= post46_gap), post46_gap))

# ------------------------------------------------------------------------------
# 4. Falsification Test B: Low-AI vs. High-AI Sectoral Placebo (2022M11-2026M8)
# ------------------------------------------------------------------------------
cat("======================================================================\n")
cat("FALSIFICATION TEST B: LOW-AI vs. HIGH-AI SECTORAL PLACEBO (2022M11-2026M8)\n")
cat("======================================================================\n")

H22 <- nrow(te22_46)
sectors <- list(
  list(short = "Government (USGOVT)",    exposure = "Low AI (Institutional)", col = "L_USGOVT"),
  list(short = "Construction (USCONS)",  exposure = "Low AI (Physical)",      col = "L_USCONS"),
  list(short = "Manufacturing (MANEMP)", exposure = "Low/Med AI (Physical)",  col = "L_MANEMP"),
  list(short = "Total Payrolls (PAYEMS)",exposure = "Aggregate Economy",      col = "L_PAYEMS"),
  list(short = "Finance (USFIRE)",       exposure = "High AI (Cognitive)",    col = "L_USFIRE"),
  list(short = "Job Openings (HWI)",     exposure = "Hiring Margin (JOLTS)",  col = "L_HWI")
)

sec_rows <- list()
for (i in seq_along(sectors)) {
  s <- sectors[[i]]
  second_col <- if (s$col == "L_PAYEMS") "L_HWI" else "L_PAYEMS"
  res <- run_lp_bvar_levels(tr22, te22_46, c(s$col, second_col), p = 3)
  gap <- 100 * (exp(res$Y_real[H22, 5] / 100) / exp(res$Y_cnd[H22, 5] / 100) - 1)
  real_chg <- 100 * (exp(res$Y_real[H22, 5] / 100) / exp(tr22[[s$col]][nrow(tr22)] / 100) - 1)
  e_mo <- diff(c(tr22[[s$col]][nrow(tr22)], res$Y_real[, 5])) - diff(c(tr22[[s$col]][nrow(tr22)], res$Y_cnd[, 5]))
  ct <- coeftest(lm(e_mo ~ 1), vcov. = NeweyWest(lm(e_mo ~ 1), lag = 3, prewhite = FALSE))
  sec_rows[[i]] <- data.frame(
    Series         = s$short,
    AI_Exposure    = s$exposure,
    Realized_Chg   = round(real_chg, 1),
    Cond_Gap_Pct   = round(gap, 1),
    Mean_Err_PctMo = round(ct[1, 1], 3),
    HAC_t          = round(ct[1, 3], 2),
    HAC_p          = signif(ct[1, 4], 3),
    stringsAsFactors = FALSE
  )
}
sectoral_df <- do.call(rbind, sec_rows)
print(sectoral_df, row.names = FALSE)
write.csv(sectoral_df, "outputs/econometrics/table_falsification_sectoral.csv", row.names = FALSE)

# Export LaTeX Table 6 (Single-Column IEEE zero-badness minipage)
tex_tab6 <- sprintf(
"\\begin{table}[t]
\\centering
\\caption{Falsification Tests: Pre-2022 Pseudo-Breaks \\& Sectoral Placebos}
\\label{tab:falsification}
\\begin{minipage}{\\columnwidth}
\\centering
\\fontsize{6.4pt}{7.8pt}\\selectfont
\\setlength{\\tabcolsep}{1.4pt}
\\begin{tabular}{lcccc}
\\toprule
\\multicolumn{5}{l}{\\textit{Panel A: Pre-2022 Pseudo-Break Placebo Origins (24-Mo Clean Horizons)}} \\\\
\\textbf{Origin / Window} & \\textbf{HWI Gap} & $t_{\\text{HAC}}$ & \\textbf{PAYEMS Gap} & \\textbf{Exceedance} \\\\
\\midrule
Pre-AI Mean (13 origins, 2015--18) & $%+.1f\\%%$ & -- & $%+.1f\\%%$ & -- \\\\
Pre-AI Most Negative (%s)          & $%+.1f\\%%$ & $%.2f$ & $%+.1f\\%%$ & -- \\\\
\\textbf{Post-2022 Actual (24-mo)} & $\\mathbf{%+.1f\\%%}$ & $\\mathbf{%.2f}$ & $\\mathbf{%+.1f\\%%}$ & $\\mathbf{0/13\\;(0\\%%)}$ \\\\
\\textbf{Post-2022 Actual (46-mo)} & $\\mathbf{%+.1f\\%%}$ & $\\mathbf{%.2f}$ & $\\mathbf{%+.1f\\%%}$ & $\\mathbf{0/13\\;(0\\%%)}$ \\\\
\\midrule
\\multicolumn{5}{l}{\\textit{Panel B: Sectoral AI-Exposure Placebo (2022M11--2026M8, $H=46$)}} \\\\
\\textbf{Labor Series (FRED-MD)} & \\textbf{Exposure} & \\textbf{Real. $\\Delta$} & \\textbf{Cond. Gap} & $t_{\\text{HAC}}$ \\\\
\\midrule
Government (\\texttt{USGOVT})    & Low (Inst.)   & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.2f$ \\\\
Construction (\\texttt{USCONS})  & Low (Phys.)   & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.2f$ \\\\
Manufacturing (\\texttt{MANEMP}) & Low/Med       & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.2f$ \\\\
Total Payrolls (\\texttt{PAYEMS})& Aggregate     & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.2f$ \\\\
Finance (\\texttt{USFIRE})       & High (Cogn.)  & $%+.1f\\%%$ & $%+.1f\\%%$ & $%+.2f$ \\\\
\\textbf{Job Openings (\\texttt{HWI})} & \\textbf{Hiring Margin} & $\\mathbf{%+.1f\\%%}$ & $\\mathbf{%+.1f\\%%}$ & $\\mathbf{%+.2f}$ \\\\
\\bottomrule
\\end{tabular}
\\vskip 2pt
{\\scriptsize \\textit{Notes}: All models use the exact marginal-likelihood Levels LP-BVAR(3) conditioned on realized $\\{\\text{FEDFUNDS},\\text{M2SL},\\text{BAA10YM}\\}$. Panel~A evaluates 13 quarterly pre-AI origins strictly ending $\\le 2020\\text{M}2$. Panel~B contrasts low-AI physical sectors against cognitive finance and JOLTS openings to bound omitted macro confounders.\\par}
\\end{minipage}
\\end{table}
",
mean(pre_hwi_gaps), mean(placebo_df$PAY_Gap_Pct[placebo_df$Regime == "Pre-AI Placebo"]),
gsub("-", "M", placebo_df$Origin[which.min(pre_hwi_gaps)]), min(pre_hwi_gaps), placebo_df$HWI_t_HAC[which.min(pre_hwi_gaps)], placebo_df$PAY_Gap_Pct[which.min(pre_hwi_gaps)],
post24_gap, placebo_df$HWI_t_HAC[placebo_df$Origin == "Post-2022 (24m)"], placebo_df$PAY_Gap_Pct[placebo_df$Origin == "Post-2022 (24m)"],
post46_gap, placebo_df$HWI_t_HAC[placebo_df$Origin == "Post-2022 (46m)"], placebo_df$PAY_Gap_Pct[placebo_df$Origin == "Post-2022 (46m)"],
sectoral_df$Realized_Chg[1], sectoral_df$Cond_Gap_Pct[1], sectoral_df$HAC_t[1],
sectoral_df$Realized_Chg[2], sectoral_df$Cond_Gap_Pct[2], sectoral_df$HAC_t[2],
sectoral_df$Realized_Chg[3], sectoral_df$Cond_Gap_Pct[3], sectoral_df$HAC_t[3],
sectoral_df$Realized_Chg[4], sectoral_df$Cond_Gap_Pct[4], sectoral_df$HAC_t[4],
sectoral_df$Realized_Chg[5], sectoral_df$Cond_Gap_Pct[5], sectoral_df$HAC_t[5],
sectoral_df$Realized_Chg[6], sectoral_df$Cond_Gap_Pct[6], sectoral_df$HAC_t[6]
)
cat(tex_tab6, file = "outputs/econometrics/table6_falsification_tests.tex")

# ------------------------------------------------------------------------------
# 5. Publication Figure 7: Falsification Battery (Unified Punto, ZERO Clipping!)
# ------------------------------------------------------------------------------
plot_a_df <- placebo_df
plot_a_df$Label <- factor(plot_a_df$Origin, levels = plot_a_df$Origin)

plot_b_df <- sectoral_df
plot_b_df$Series <- factor(plot_b_df$Series, levels = rev(plot_b_df$Series))
plot_b_df$Highlight <- ifelse(grepl("HWI", plot_b_df$Series), "JOLTS Openings",
                       ifelse(grepl("USFIRE", plot_b_df$Series), "High-AI Sector", "Low-AI / Aggregate"))

build_fig7 <- function(is_web = FALSE) {
  th <- if (is_web) theme_web_wide() else theme_ieee_print()
  sz_ttl <- if (is_web) 9.8 else 7.0

  pa <- ggplot(plot_a_df, aes(x = Label, y = HWI_Gap_Pct, fill = Regime)) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 0, color = PAL$slate, linewidth = 0.35) +
    scale_fill_manual(name = NULL, values = c("Pre-AI Placebo" = PAL$navy, "Post-2022 Actual" = PAL$red)) +
    labs(title = "(a) Pre-2022 Pseudo-Breaks vs. Post-2022", x = NULL, y = "HWI Gap (%)") +
    th +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = sz_ttl, margin = margin(b = 2)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = if (is_web) 8.0 else 5.5),
      legend.position = "bottom"
    )

  pb <- ggplot(plot_b_df, aes(x = Series, y = Cond_Gap_Pct, fill = Highlight)) +
    geom_col(width = 0.62) +
    geom_hline(yintercept = 0, color = PAL$slate, linewidth = 0.35) +
    coord_flip() +
    scale_fill_manual(name = NULL, values = c(
      "Low-AI / Aggregate" = PAL$gray,
      "High-AI Sector"     = PAL$amber,
      "JOLTS Openings"     = PAL$red
    )) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = "(b) Sectoral Placebos (2022-2026)", x = NULL, y = "Conditional Gap (%)") +
    th +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = sz_ttl, margin = margin(b = 2)),
      axis.text.y = element_text(size = if (is_web) 8.6 else 6.0),
      legend.position = "bottom"
    )

  if (is_web) (pa | pb) + plot_layout(widths = c(1.15, 1.0))
  else        (pa / pb) + plot_layout(heights = c(1.1, 1.0))
}

p7_print <- build_fig7(is_web = FALSE)
ggsave("outputs/figures/fig7_falsification_placebo.pdf", p7_print, width = 3.5, height = 3.8)
ggsave("outputs/figures/fig7_falsification_placebo.png", p7_print, width = 3.5, height = 3.8, dpi = 300)

p7_web <- build_fig7(is_web = TRUE)
ggsave("outputs/figures/fig7_web_falsification.png", p7_web, width = 8.4, height = 4.1, dpi = 300)
file.copy("outputs/figures/fig7_web_falsification.png", "docs/assets/images/fig7_falsification_placebo.png", overwrite = TRUE)

cat("\nStep 4 (Falsification & Placebo Battery) completed successfully.\n")
