# Isolating the Generative AI Labor Shock from Monetary Tightening

**Authors:** Eda Kal Arıca & Esra Şekerci
**Affiliation:** Department of Statistics, Middle East Technical University (STAT 563 — Track A Empirical Project)
**Project Website:** https://esrasekerci.github.io/generative_ai_labor_shock/

---

## Abstract

This project investigates whether systematic post-2022 forecast errors in U.S. nonfarm employment (`PAYEMS`) and job openings (`HWI`) reveal a structural break attributable to generative AI once contemporaneous monetary tightening (`FEDFUNDS`, `M2SL`, `BAA10YM`) and industrial production (`INDPRO`) are accounted for. Using monthly FRED-MD data from January 2010 onward, we compare classical econometric models (Unrestricted VAR and Minnesota-prior Bayesian VAR with Cholesky IRFs and Waggoner–Zha conditional forecasting) against learning-based architectures (Temporal Fusion Transformer / LSTM with known future covariates and zero-shot foundation models such as Chronos or TimesFM).

Full data documentation, interactive figures, econometric specifications, and references are published on the [project website](https://esrasekerci.github.io/generative_ai_labor_shock/).

---

## Repository Contents

```text
.
├── fred_data.csv                                                # FRED-MD raw dataset (St. Louis Fed)
├── requirements.txt                                             # Python environment dependencies
├── propejct_proposal/project_proposal_edakalarica_esrasekerci.pdf
├── src/data/prepare_web_data.py                                 # Compiles fred_data.csv -> docs/assets/data/fred_macro.js
└── docs/                                                        # GitHub Pages static site (main -> /docs)
    ├── .nojekyll
    ├── index.html
    └── assets/
        ├── css/style.css
        ├── js/app.js
        └── data/fred_macro.js
```

---

## Reproducing Web Data & Local Preview

To regenerate the processed dataset after updating `fred_data.csv`:

```bash
python3 src/data/prepare_web_data.py
```

To preview the site locally:

```bash
python3 -m http.server 8000 --directory docs
```
