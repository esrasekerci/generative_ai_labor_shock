#!/usr/bin/env python3
"""
Prepare structured JSON & JS datasets for the STAT 563 GitHub Pages site
directly from the local FRED-MD raw CSV (`fred_data.csv`).

Uses only the Python standard library so it runs reproducibly anywhere.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev

REPO_ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = REPO_ROOT / "fred_data.csv"
OUT_DIR = REPO_ROOT / "docs" / "assets" / "data"

SERIES_META = [
    {
        "id": "INDPRO",
        "label": "Industrial Production Index",
        "shortLabel": "INDPRO",
        "role": "Real Economic Activity",
        "block": "Real Activity",
        "units": "Index (2017 = 100)",
        "tcode": 5,
        "transformFormula": "Δ ln(xₜ)",
        "transformLabel": "Monthly log-difference (% growth × 100)",
        "color": "#1D4E89",
    },
    {
        "id": "PAYEMS",
        "label": "All Employees, Total Nonfarm",
        "shortLabel": "PAYEMS",
        "role": "Aggregate Employment (Target)",
        "block": "Labor Market Target",
        "units": "Thousands of Persons",
        "tcode": 5,
        "transformFormula": "Δ ln(xₜ)",
        "transformLabel": "Monthly log-difference (% growth × 100)",
        "color": "#D97706",
    },
    {
        "id": "HWI",
        "label": "JOLTS Total Nonfarm Job Openings (HWI / JTSJOL)",
        "shortLabel": "Job Openings (HWI)",
        "role": "Labor Demand / Vacancies (Target)",
        "block": "Labor Market Target",
        "units": "Thousands of Openings",
        "tcode": 2,
        "transformFormula": "Δ xₜ",
        "transformLabel": "Monthly first difference",
        "color": "#15803D",
    },
    {
        "id": "FEDFUNDS",
        "label": "Effective Federal Funds Rate",
        "shortLabel": "FEDFUNDS",
        "role": "Monetary Policy Instrument (Conditioning)",
        "block": "Policy & Financial Conditioning",
        "units": "% per annum",
        "tcode": 2,
        "transformFormula": "Δ xₜ",
        "transformLabel": "Monthly first difference (pp)",
        "color": "#B5121B",
    },
    {
        "id": "M2SL",
        "label": "M2 Money Stock",
        "shortLabel": "M2SL",
        "role": "Broad Liquidity Aggregate (Conditioning)",
        "block": "Policy & Financial Conditioning",
        "units": "Billions of Dollars",
        "tcode": 6,
        "transformFormula": "Δ² ln(xₜ)",
        "transformLabel": "Second log-difference (× 100)",
        "color": "#6D28D9",
    },
    {
        "id": "BAA10YM",
        "label": "Moody's Baa – 10Y Treasury Spread",
        "shortLabel": "BAA10YM",
        "role": "Credit Risk / Financial Stress (Conditioning)",
        "block": "Policy & Financial Conditioning",
        "units": "% per annum (BAA − GS10)",
        "tcode": 1,
        "transformFormula": "xₜ",
        "transformLabel": "Level (stationary yield spread, pp)",
        "color": "#7C2D12",
    },
]


def parse_float(val: str) -> float | None:
    s = val.strip()
    if not s:
        return None
    return float(s)


def parse_date(s: str) -> str:
    dt = datetime.strptime(s.strip(), "%m/%d/%Y")
    return dt.strftime("%Y-%m-%d")


def compute_transform(values: list[float | None], tcode: int) -> list[float | None]:
    """Compute FRED-MD stationary transformation (McCracken & Ng, 2016)."""
    n = len(values)
    out: list[float | None] = [None] * n
    if tcode == 1:
        return list(values)
    if tcode == 2:
        for i in range(1, n):
            if values[i] is not None and values[i - 1] is not None:
                out[i] = values[i] - values[i - 1]
        return out
    if tcode == 5:
        for i in range(1, n):
            if (
                values[i] is not None
                and values[i - 1] is not None
                and values[i] > 0
                and values[i - 1] > 0
            ):
                out[i] = 100.0 * (math.log(values[i]) - math.log(values[i - 1]))
        return out
    if tcode == 6:
        for i in range(2, n):
            if (
                values[i] is not None
                and values[i - 1] is not None
                and values[i - 2] is not None
                and values[i] > 0
                and values[i - 1] > 0
                and values[i - 2] > 0
            ):
                d1 = math.log(values[i]) - math.log(values[i - 1])
                d0 = math.log(values[i - 1]) - math.log(values[i - 2])
                out[i] = 100.0 * (d1 - d0)
        return out
    raise ValueError(f"Unsupported tcode: {tcode}")


def pearson_corr(x: list[float | None], y: list[float | None]) -> float | None:
    pairs = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 3:
        return None
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    mx = mean(xs)
    my = mean(ys)
    num = sum((a - mx) * (b - my) for a, b in pairs)
    den_x = math.sqrt(sum((a - mx) ** 2 for a in xs))
    den_y = math.sqrt(sum((b - my) ** 2 for b in ys))
    if den_x == 0 or den_y == 0:
        return None
    return round(num / (den_x * den_y), 3)


def summarize(vals: list[float | None], dates: list[str]) -> dict:
    valid = [(d, v) for d, v in zip(dates, vals) if v is not None]
    nums = [v for _, v in valid]
    missing_dates = [d for d, v in zip(dates, vals) if v is None]
    return {
        "nValid": len(nums),
        "nMissing": len(missing_dates),
        "missingDates": missing_dates,
        "mean": round(mean(nums), 3),
        "sd": round(stdev(nums), 3),
        "min": round(min(nums), 3),
        "max": round(max(nums), 3),
        "startVal": round(valid[0][1], 3),
        "endVal": round(valid[-1][1], 3),
        "lastDate": valid[-1][0],
    }


def main() -> None:
    with CSV_PATH.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        _tcode_row = next(reader)
        raw_rows = [row for row in reader if row and row[0].strip()]

    col_idx = {name: header.index(name) for name in ["sasdate", "INDPRO", "PAYEMS", "HWI", "FEDFUNDS", "M2SL", "BAA", "GS10"]}

    # Parse full series first so lag-1 / lag-2 transforms at 2010-01-01 use 2009-11/12 cleanly
    all_dates: list[str] = []
    full_levels: dict[str, list[float | None]] = {m["id"]: [] for m in SERIES_META}

    for row in raw_rows:
        dt = parse_date(row[col_idx["sasdate"]])
        all_dates.append(dt)
        for sid in ["INDPRO", "PAYEMS", "HWI", "FEDFUNDS", "M2SL"]:
            full_levels[sid].append(parse_float(row[col_idx[sid]]))
        baa = parse_float(row[col_idx["BAA"]])
        gs10 = parse_float(row[col_idx["GS10"]])
        full_levels["BAA10YM"].append(round(baa - gs10, 4) if (baa is not None and gs10 is not None) else None)

    full_transformed: dict[str, list[float | None]] = {}
    for m in SERIES_META:
        full_transformed[m["id"]] = compute_transform(full_levels[m["id"]], m["tcode"])

    # Slice Jan 2010 onward
    start_idx = next(i for i, d in enumerate(all_dates) if d >= "2010-01-01")
    dates = all_dates[start_idx:]

    train_cutoff = "2022-10-01"
    eval_start = "2022-11-01"

    n_total = len(dates)
    n_train = sum(1 for d in dates if d <= train_cutoff)
    n_test = sum(1 for d in dates if d >= eval_start)

    levels: dict[str, list[float | None]] = {}
    zscores: dict[str, list[float | None]] = {}
    transformed: dict[str, list[float | None]] = {}
    stats: dict[str, dict] = {}
    regime_stats: dict[str, dict] = {}

    for m in SERIES_META:
        sid = m["id"]
        lv = full_levels[sid][start_idx:]
        tr = full_transformed[sid][start_idx:]
        valid_lv = [v for v in lv if v is not None]
        mu = mean(valid_lv)
        sd = stdev(valid_lv)
        zs = [round((v - mu) / sd, 4) if v is not None else None for v in lv]

        levels[sid] = [round(v, 4) if v is not None else None for v in lv]
        zscores[sid] = zs
        transformed[sid] = [round(v, 4) if v is not None else None for v in tr]
        stats[sid] = summarize(lv, dates)

        pre_vals = [v for d, v in zip(dates, lv) if d <= train_cutoff and v is not None]
        post_vals = [v for d, v in zip(dates, lv) if d >= eval_start and v is not None]
        pre_tr = [v for d, v in zip(dates, tr) if d <= train_cutoff and v is not None]
        post_tr = [v for d, v in zip(dates, tr) if d >= eval_start and v is not None]
        regime_stats[sid] = {
            "preMean": round(mean(pre_vals), 3),
            "postMean": round(mean(post_vals), 3),
            "preTransformMean": round(mean(pre_tr), 3),
            "postTransformMean": round(mean(post_tr), 3),
        }

    ids = [m["id"] for m in SERIES_META]
    corr_levels = [[pearson_corr(levels[a], levels[b]) for b in ids] for a in ids]
    corr_transformed = [[pearson_corr(transformed[a], transformed[b]) for b in ids] for a in ids]

    # Pre-AI vs Post-AI correlations on transformed series
    pre_idx = [i for i, d in enumerate(dates) if d <= train_cutoff]
    post_idx = [i for i, d in enumerate(dates) if d >= eval_start]
    corr_pre_tr = [
        [
            pearson_corr([transformed[a][i] for i in pre_idx], [transformed[b][i] for i in pre_idx])
            for b in ids
        ]
        for a in ids
    ]
    corr_post_tr = [
        [
            pearson_corr([transformed[a][i] for i in post_idx], [transformed[b][i] for i in post_idx])
            for b in ids
        ]
        for a in ids
    ]

    payload = {
        "metadata": {
            "source": "FRED-MD (Federal Reserve Bank of St. Louis)",
            "sourceUrl": "https://www.stlouisfed.org/research/economists/mccracken/fred-databases",
            "generatedFrom": "fred_data.csv",
            "frequency": "Monthly",
            "k": len(SERIES_META),
            "sampleStart": dates[0],
            "sampleEnd": dates[-1],
            "nTotal": n_total,
            "trainWindow": {"start": dates[0], "end": train_cutoff, "n": n_train},
            "evalWindow": {"start": eval_start, "end": dates[-1], "n": n_test},
            "proposalTarget": {"end": "2026-10-01", "nTotal": 202, "nTrain": 154, "nTest": 48},
            "covidWindow": {"start": "2020-03-01", "end": "2020-12-01"},
            "dataNotes": [
                "BAA10YM is constructed directly from FRED-MD columns as BAA minus GS10 (zero missing observations).",
                "The local fred_data.csv snapshot runs through September 2026 (201 monthly observations from Jan 2010), one month prior to the proposal's Oct 2026 target (T = 202), and HWI has a single trailing release lag (missing 2026-09-01).",
            ],
        },
        "series": SERIES_META,
        "dates": dates,
        "levels": levels,
        "zscores": zscores,
        "transformed": transformed,
        "stats": stats,
        "regimeStats": regime_stats,
        "correlations": {
            "labels": ids,
            "levels": corr_levels,
            "transformed": corr_transformed,
            "preAITransformed": corr_pre_tr,
            "postAITransformed": corr_post_tr,
        },
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    js_path = OUT_DIR / "fred_macro.js"

    json_str = json.dumps(payload, indent=2)
    js_path.write_text(f"window.FRED_MACRO_DATA = {json_str};\n", encoding="utf-8")

    print(f"Wrote {js_path.relative_to(REPO_ROOT)}")
    print(f"Sample: {dates[0]} to {dates[-1]} (T={n_total}; train={n_train}, eval={n_test})")
    for sid in ids:
        st = stats[sid]
        rs = regime_stats[sid]
        print(
            f"  {sid:8s} | n={st['nValid']:3d} (miss={st['nMissing']}) | "
            f"mean={st['mean']:9.3f} | sd={st['sd']:8.3f} | "
            f"pre_mean={rs['preMean']:9.3f} | post_mean={rs['postMean']:9.3f}"
        )


if __name__ == "__main__":
    main()
