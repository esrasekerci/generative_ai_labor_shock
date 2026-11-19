# Model Implementation Roadmap

This directory houses the four model families specified in Section 3 of the STAT 563 Track A project proposal (`propejct_proposal/project_proposal_edakalarica_esrasekerci.pdf`).

## Planned Modules

1. **`var_baseline.py`** — Unrestricted VAR($p$)
   - Stationary FRED-MD transformations (`tcode`: 5 for `INDPRO`/`PAYEMS`, 2 for `HWI`/`FEDFUNDS`, 6 for `M2SL`, 1 for `BAA10YM = BAA - GS10`).
   - Lag order $p$ selected via BIC on the pre-2022 training sample (`2010-01` – `2022-10`, $T_{\text{train}} = 154$).
   - Unconditional multi-step baseline forecasts over the post-2022 evaluation window.

2. **`bvar_conditional.py` (or R counterpart)** — Bayesian VAR with Minnesota Priors
   - Shrinkage via Minnesota prior hyperparameters to control over-parameterization.
   - Structural validation through Cholesky-identified Impulse Response Functions (IRFs) verifying standard monetary transmission.
   - Conditional forecasting following Waggoner & Zha (1999) and Bańbura, Giannone, & Lenza (2015), conditioning on the realized post-2022 trajectories of `FEDFUNDS`, `M2SL`, and `BAA10YM` to isolate unexplained forecast errors in `PAYEMS` and `HWI`.

3. **`deep_learning_tft_lstm.py`** — Temporal Fusion Transformer (TFT) / LSTM
   - Trained exclusively on pre-2022 multivariate data (`2010-01` – `2022-10`).
   - Treats post-2022 monetary and financial variables (`FEDFUNDS`, `M2SL`, `BAA10YM`) as known future covariates.
   - Hyperparameter selection (hidden dimension, layers, learning rate, dropout) via expanding-window time-series cross-validation.

4. **`foundation_zero_shot.py`** — Zero-Shot Foundation Models
   - Evaluates pretrained time-series foundation models (Chronos, Ansari et al., 2024; TimesFM) prompted with the pre-2022 context without task-specific fine-tuning.
   - Benchmarks zero-shot structural-break sensitivity against the conditional BVAR and trained deep learning models.
