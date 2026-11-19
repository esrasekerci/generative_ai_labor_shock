# Isolating the Generative AI Labor Shock from Monetary Tightening

**Middle East Technical University — Department of Statistics**
**STAT 563: Multivariate Time Series Analysis · Track A (Empirical Project)**
**Authors:** Eda Kal Arıca & Esra Şekerci

[![Track A](https://img.shields.io/badge/STAT%20563-Track%20A%20Empirical-B5121B?style=flat-square)](#)
[![Data: FRED-MD](https://img.shields.io/badge/Data-FRED--MD%20(k%3D6)-102A43?style=flat-square)](https://www.stlouisfed.org/research/economists/mccracken/fred-databases)
[![GitHub Pages](https://img.shields.io/badge/Live%20Site-GitHub%20Pages-15803D?style=flat-square)](https://esrasekerci.github.io/generative_ai_labor_shock/)

**Live Project Website:** https://esrasekerci.github.io/generative_ai_labor_shock/

---

## Research Question

> *Does the systematic forecast error in U.S. employment (`PAYEMS`) and job openings (`HWI`) post-2022 reveal a structural break attributable to the generative AI shock, and how do classical Bayesian VARs compare to modern deep learning and zero-shot foundation models in isolating this effect?*

The public release of generative AI in late 2022 coincided directly with the Federal Reserve's monetary tightening cycle. Because simple intervention dummies cannot separate simultaneous macroeconomic shocks, this project uses a **multivariate conditional forecasting framework** across classical econometric and modern learning-based architectures to isolate endogenous policy-induced cooling from any exogenous structural break in U.S. labor demand.

---

## Dataset (`FRED-MD`, $k = 6$)

Source: [Federal Reserve Economic Data — Monthly Database (FRED-MD)](https://www.stlouisfed.org/research/economists/mccracken/fred-databases) (McCracken & Ng, 2016).

| Series | Description | Role in System | Units | FRED-MD `tcode` |
| :--- | :--- | :--- | :--- | :--- |
| `INDPRO` | Industrial Production Index | Real economic activity | Index (`2017 = 100`) | `5` ($\Delta \ln x_t$) |
| `PAYEMS` | All Employees, Total Nonfarm | Aggregate employment (**Target**) | Thousands of persons | `5` ($\Delta \ln x_t$) |
| `HWI` | Help-Wanted Index for United States | Job openings proxy (**Target**) | Index / Thousands | `2` ($\Delta x_t$) |
| `FEDFUNDS` | Effective Federal Funds Rate | Monetary policy instrument (**Conditioning**) | `% p.a.` | `2` ($\Delta x_t$) |
| `M2SL` | M2 Money Stock | Broad liquidity aggregate (**Conditioning**) | Billions of dollars | `6` ($\Delta^2 \ln x_t$) |
| `BAA10YM` | Moody's Seasoned Baa minus 10Y Treasury | Credit risk / financial stress (**Conditioning**) | `% p.a.` (`BAA − GS10`) | `1` ($x_t$, level) |

- **Training Window:** January 2010 – October 2022 ($T_{\text{train}} = 154$ monthly observations)
- **Post-2022 Evaluation Window:** November 2022 – October 2026 ($T_{\text{test}} = 48$ target; local `fred_data.csv` snapshot currently extends through September 2026, $T_{\text{test}} = 47$, with `HWI` trailing by one month due to release lag)

---

## Planned Model Comparison

1. **Unrestricted VAR($p$)** — Stationary baseline with BIC lag selection.
2. **Bayesian VAR (BVAR)** — Minnesota prior shrinkage, Cholesky-identified IRFs to verify monetary transmission, and **Waggoner & Zha (1999) conditional forecasting** fixing realized post-2022 paths of `FEDFUNDS`, `M2SL`, and `BAA10YM`.
3. **Deep Learning (TFT / LSTM)** — Trained on pre-2022 multivariate series treating post-2022 monetary/financial variables as known future covariates, tuned via time-series cross-validation.
4. **Zero-Shot Foundation Model (Chronos / TimesFM)** — Pretrained time-series foundation model prompted zero-shot on the pre-2022 context to benchmark structural-break handling without parametric restrictions.

---

## Repository Structure

```text
.
├── README.md
├── .gitignore
├── requirements.txt
├── fred_data.csv                           # Local FRED-MD raw dataset
├── propejct_proposal/                      # Submitted STAT 563 proposal & guidelines PDFs
├── src/
│   ├── data/
│   │   └── prepare_web_data.py             # Reproducible script exporting fred_data.csv -> docs/assets/data/
│   └── models/
│       └── README.md                       # Specification guide for VAR, BVAR, TFT/LSTM, Chronos/TimesFM
├── notebooks/                              # Exploratory & modeling notebooks
└── docs/                                   # GitHub Pages website (deploy from main -> /docs)
    ├── .nojekyll
    ├── index.html                          # Interactive editorial academic project website
    └── assets/
        ├── css/style.css
        ├── js/app.js
        └── data/
            ├── fred_macro.json             # Processed dataset, stats, & regime correlation matrices
            └── fred_macro.js               # Offline-compatible JS bundle for local file:// preview
```

---

## Regenerating Web Data & Previewing Locally

To refresh the interactive GitHub Pages datasets after updating `fred_data.csv`:

```bash
python3 src/data/prepare_web_data.py
```

To preview the GitHub Pages site locally:

```bash
python3 -m http.server 8000 --directory docs
# Open http://localhost:8000
```

You can also open `docs/index.html` directly in your browser.

---

## Deploying to GitHub Pages

This repository is deployed via **GitHub Pages** from branch `main` and folder `/docs`:

- **Repository:** https://github.com/esrasekerci/generative_ai_labor_shock
- **Live Site:** https://esrasekerci.github.io/generative_ai_labor_shock/

Any push to `main` (via GitHub Desktop or CLI) automatically updates the live site.

---

## References

Key references include Ansari et al. (2024), Bańbura, Giannone, & Lenza (2015), Caldara et al. (2023), Crump et al. (2021), Debortoli et al. (2023), Eksi & Tas (2022), Korobilis (2025), McCracken & Ng (2016), Pater (2017), Primiceri (2005), Sznajderska & Haug (2023), and Waggoner & Zha (1999). Full bibliography is available on the [project website](docs/index.html#references).
