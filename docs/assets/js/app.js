/* ==========================================================================
   STAT 563 Editorial Research Site — Interactive Data & Visualization Engine
   Reads pre-compiled FRED-MD payload from `window.FRED_MACRO_DATA`
   (generated via `src/data/prepare_web_data.py` from `fred_data.csv`).
   ========================================================================== */

(function () {
  "use strict";

  const DATA = window.FRED_MACRO_DATA;
  if (!DATA) {
    console.error("FRED_MACRO_DATA not found. Run `python3 src/data/prepare_web_data.py`.");
    return;
  }

  const FONT_SANS = '"Source Sans 3", -apple-system, BlinkMacSystemFont, sans-serif';
  const FONT_MONO = '"JetBrains Mono", monospace';
  const PLOTLY_CONFIG = {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
  };

  const state = {
    fig1Scale: "zscore",   // "zscore" | "indexed"
    fig1Range: "full",     // "full" | "covid" | "postai"
    fig1Subset: "all",     // "all" | "labor" | "cond"
    fig2Series: "HWI",     // default to HWI (flagship series in hypothesis)
    fig3Mode: "preAITransformed",
  };

  /* ------------------------- Formatting Helpers ------------------------- */
  function fmtNumber(val, id) {
    if (val === null || val === undefined || Number.isNaN(val)) return "—";
    if (id === "PAYEMS" || id === "M2SL" || id === "HWI") {
      return val.toLocaleString("en-US", { maximumFractionDigits: 1, minimumFractionDigits: 1 });
    }
    return val.toLocaleString("en-US", { maximumFractionDigits: 2, minimumFractionDigits: 2 });
  }

  function fmtDateShort(iso) {
    const [y, m] = iso.split("-");
    const months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    return `${months[parseInt(m, 10) - 1]} ${y}`;
  }

  function getRegimeShapes(yRef = "paper") {
    const covid = DATA.metadata.covidWindow;
    const evalWin = DATA.metadata.evalWindow;
    return [
      {
        type: "rect",
        xref: "x",
        yref: yRef,
        x0: covid.start,
        x1: covid.end,
        y0: 0,
        y1: 1,
        fillcolor: "rgba(148, 163, 184, 0.22)",
        line: { width: 0 },
        layer: "below",
      },
      {
        type: "rect",
        xref: "x",
        yref: yRef,
        x0: evalWin.start,
        x1: evalWin.end,
        y0: 0,
        y1: 1,
        fillcolor: "rgba(245, 158, 11, 0.12)",
        line: { width: 0 },
        layer: "below",
      },
      {
        type: "line",
        xref: "x",
        yref: yRef,
        x0: evalWin.start,
        x1: evalWin.start,
        y0: 0,
        y1: 1,
        line: { color: "#B5121B", width: 1.5, dash: "dot" },
      },
    ];
  }

  /* --------------------- Populate Summary Table & KPIs ------------------- */
  function populateTableAndKPIs() {
    const meta = DATA.metadata;
    document.getElementById("kpi-k").textContent = `k = ${meta.k}`;
    document.getElementById("kpi-total").textContent = `T = ${meta.nTotal}`;
    document.getElementById("kpi-span").textContent =
      `${fmtDateShort(meta.sampleStart)} – ${fmtDateShort(meta.sampleEnd)} (${meta.proposalTarget.nTotal} target)`;
    document.getElementById("kpi-train").innerHTML = `T<sub>train</sub> = ${meta.trainWindow.n}`;
    document.getElementById("kpi-test").innerHTML = `T<sub>test</sub> = ${meta.evalWindow.n}`;
    document.getElementById("kpi-test-sub").textContent =
      `${fmtDateShort(meta.evalWindow.start)} – ${fmtDateShort(meta.evalWindow.end)} (${meta.proposalTarget.nTest} target)`;

    const tbody = document.getElementById("series-table-body");
    tbody.innerHTML = "";

    DATA.series.forEach((s) => {
      const st = DATA.stats[s.id];
      const rs = DATA.regimeStats[s.id];
      let roleClass = "real";
      let roleShort = "Real Activity";
      if (s.block === "Labor Market Target") {
        roleClass = "target";
        roleShort = "Labor Target";
      } else if (s.block === "Policy & Financial Conditioning") {
        roleClass = "cond";
        roleShort = "Conditioning";
      }

      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>
          <span class="series-dot" style="background:${s.color};"></span>
          <code>${s.id}</code>
        </td>
        <td>
          <strong>${s.label}</strong><br/>
          <span style="color:var(--ink-muted); font-size:0.8rem;">${s.units}</span>
        </td>
        <td><span class="role-badge ${roleClass}">${roleShort}</span></td>
        <td>
          <code>tcode ${s.tcode}</code> · ${s.transformFormula}
        </td>
        <td class="num">${fmtNumber(st.mean, s.id)}</td>
        <td class="num">${fmtNumber(st.sd, s.id)}</td>
        <td class="num">${fmtNumber(rs.preMean, s.id)}</td>
        <td class="num">${fmtNumber(rs.postMean, s.id)}</td>
      `;
      tbody.appendChild(tr);
    });
  }

  /* --------------------- Figure 1: Multivariate Co-Movement -------------- */
  function renderFigure1() {
    const dates = DATA.dates;
    const baseIdx = dates.indexOf("2022-10-01");

    let activeIds = DATA.series.map((s) => s.id);
    if (state.fig1Subset === "labor") {
      activeIds = ["PAYEMS", "HWI", "FEDFUNDS", "INDPRO"];
    } else if (state.fig1Subset === "cond") {
      activeIds = ["FEDFUNDS", "M2SL", "BAA10YM"];
    }

    const traces = DATA.series.filter((s) => activeIds.includes(s.id)).map((s) => {
      let yValues;
      if (state.fig1Scale === "zscore") {
        yValues = DATA.zscores[s.id];
      } else {
        const baseVal = DATA.levels[s.id][baseIdx];
        yValues = DATA.levels[s.id].map((v) =>
          v !== null && baseVal ? Number(((v / baseVal) * 100).toFixed(2)) : null
        );
      }

      return {
        type: "scatter",
        mode: "lines",
        name: s.shortLabel,
        x: dates,
        y: yValues,
        line: {
          color: s.color,
          width: s.id === "HWI" || s.id === "PAYEMS" ? 2.4 : 1.75,
        },
        hovertemplate:
          state.fig1Scale === "zscore"
            ? `<b>${s.shortLabel}</b>: %{y:+.2f} σ<extra></extra>`
            : `<b>${s.shortLabel}</b>: %{y:.1f} (Oct 2022 = 100)<extra></extra>`,
      };
    });

    // Determine x-axis range
    let xRange = [dates[0], dates[dates.length - 1]];
    if (state.fig1Range === "covid") {
      xRange = ["2019-01-01", "2022-12-01"];
    } else if (state.fig1Range === "postai") {
      xRange = ["2022-01-01", dates[dates.length - 1]];
    }

    const shapes = getRegimeShapes("paper");
    shapes.push({
      type: "line",
      xref: "paper",
      yref: "y",
      x0: 0,
      x1: 1,
      y0: state.fig1Scale === "zscore" ? 0 : 100,
      y1: state.fig1Scale === "zscore" ? 0 : 100,
      line: { color: "#94A3B8", width: 1, dash: "dash" },
      layer: "below",
    });

    const annotations = [
      {
        x: "2020-07-15",
        y: 0.99,
        xref: "x",
        yref: "paper",
        text: "COVID-19 Outlier",
        showarrow: false,
        font: { family: FONT_MONO, size: 10.5, color: "#475569" },
        yanchor: "top",
      },
      {
        x: "2024-08-01",
        y: 0.99,
        xref: "x",
        yref: "paper",
        text: "Post-AI & Tightening Window",
        showarrow: false,
        font: { family: FONT_MONO, size: 10.5, color: "#92400E" },
        yanchor: "top",
      },
    ];

    const layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#FFFFFF",
      font: { family: FONT_SANS, color: "#17202A", size: 12.5 },
      margin: { l: 58, r: 20, t: 22, b: 44 },
      hovermode: "x unified",
      xaxis: {
        range: xRange,
        showgrid: true,
        gridcolor: "#EEF2F6",
        zeroline: false,
        tickfont: { family: FONT_MONO, size: 11 },
      },
      yaxis: {
        title: {
          text:
            state.fig1Scale === "zscore"
              ? "Standard Deviations from Mean (Z-Score)"
              : "Index Level (October 2022 = 100)",
          font: { size: 12.5, color: "#434E5E" },
        },
        showgrid: true,
        gridcolor: "#EEF2F6",
        zeroline: false,
        tickfont: { family: FONT_MONO, size: 11 },
      },
      legend: {
        orientation: "h",
        yanchor: "bottom",
        y: 1.02,
        xanchor: "left",
        x: 0,
        font: { size: 12 },
      },
      shapes,
      annotations:
        state.fig1Range === "postai"
          ? [annotations[1]]
          : state.fig1Range === "covid"
            ? [annotations[0]]
            : annotations,
    };

    Plotly.react("plot-comovement", traces, layout, PLOTLY_CONFIG);
  }

  /* --------------------- Figure 2: Single-Series Detail ------------------ */
  function buildFigure2Pills() {
    const container = document.getElementById("fig2-series-pills");
    container.innerHTML = "";

    DATA.series.forEach((s) => {
      const btn = document.createElement("button");
      btn.className = `series-pill${s.id === state.fig2Series ? " active" : ""}`;
      btn.dataset.series = s.id;
      btn.innerHTML = `<span class="series-dot" style="background:${s.color}; margin-right:0;"></span>${s.id}`;
      btn.addEventListener("click", () => {
        state.fig2Series = s.id;
        container.querySelectorAll(".series-pill").forEach((el) => el.classList.remove("active"));
        btn.classList.add("active");
        renderFigure2();
      });
      container.appendChild(btn);
    });
  }

  function renderFigure2() {
    const sid = state.fig2Series;
    const meta = DATA.series.find((s) => s.id === sid);
    const st = DATA.stats[sid];
    const rs = DATA.regimeStats[sid];
    const dates = DATA.dates;

    const oct22Idx = dates.indexOf("2022-10-01");
    const oct22Val = DATA.levels[sid][oct22Idx];
    const latestVal = st.endVal;
    let deltaText = "—";
    if (oct22Val !== null && latestVal !== null) {
      if (sid === "FEDFUNDS" || sid === "BAA10YM") {
        const diff = latestVal - oct22Val;
        deltaText = `${diff >= 0 ? "+" : ""}${diff.toFixed(2)} pp`;
      } else {
        const pct = ((latestVal - oct22Val) / oct22Val) * 100;
        deltaText = `${pct >= 0 ? "+" : ""}${pct.toFixed(1)}%`;
      }
    }

    const kpiContainer = document.getElementById("fig2-kpis");
    kpiContainer.innerHTML = `
      <div class="mini-kpi">
        <div class="mini-kpi-label">FRED-MD Transform</div>
        <div class="mini-kpi-val">tcode ${meta.tcode} (${meta.transformFormula})</div>
      </div>
      <div class="mini-kpi">
        <div class="mini-kpi-label">Pre-AI Mean (2010–22)</div>
        <div class="mini-kpi-val">${fmtNumber(rs.preMean, sid)}</div>
      </div>
      <div class="mini-kpi">
        <div class="mini-kpi-label">Post-AI Mean (2022–26)</div>
        <div class="mini-kpi-val">${fmtNumber(rs.postMean, sid)}</div>
      </div>
      <div class="mini-kpi">
        <div class="mini-kpi-label">Oct 2022 → ${fmtDateShort(st.lastDate)}</div>
        <div class="mini-kpi-val">${deltaText}</div>
      </div>
    `;

    const traceLevel = {
      type: "scatter",
      mode: "lines",
      name: `${meta.id} Level`,
      x: dates,
      y: DATA.levels[sid],
      xaxis: "x",
      yaxis: "y",
      line: { color: meta.color, width: 2.3 },
      fill: "tozeroy",
      fillcolor: "rgba(16, 42, 67, 0.04)",
      hovertemplate: `<b>Level</b>: %{y:.2f}<extra></extra>`,
    };

    const traceTrans = {
      type: meta.tcode === 1 ? "scatter" : "bar",
      mode: meta.tcode === 1 ? "lines" : undefined,
      name: `Stationary (${meta.transformFormula})`,
      x: dates,
      y: DATA.transformed[sid],
      xaxis: "x",
      yaxis: "y2",
      marker:
        meta.tcode === 1
          ? undefined
          : {
              color: DATA.transformed[sid].map((v) =>
                v === null ? "#CBD5E1" : v >= 0 ? "rgba(21, 128, 61, 0.72)" : "rgba(181, 18, 27, 0.72)"
              ),
            },
      line: meta.tcode === 1 ? { color: meta.color, width: 2 } : undefined,
      hovertemplate: `<b>${meta.transformFormula}</b>: %{y:+.3f}<extra></extra>`,
    };

    const layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#FFFFFF",
      font: { family: FONT_SANS, color: "#17202A", size: 12 },
      margin: { l: 68, r: 20, t: 18, b: 40 },
      hovermode: "x unified",
      showlegend: false,
      xaxis: {
        anchor: "y2",
        showgrid: true,
        gridcolor: "#EEF2F6",
        tickfont: { family: FONT_MONO, size: 11 },
      },
      yaxis: {
        domain: [0.54, 1.0],
        anchor: "x",
        title: { text: `Raw Level (${meta.units})`, font: { size: 11.5, color: "#434E5E" } },
        showgrid: true,
        gridcolor: "#EEF2F6",
        tickfont: { family: FONT_MONO, size: 10.5 },
      },
      yaxis2: {
        domain: [0.0, 0.42],
        anchor: "x",
        title: { text: meta.transformLabel, font: { size: 11.5, color: "#434E5E" } },
        showgrid: true,
        gridcolor: "#EEF2F6",
        zeroline: true,
        zerolinecolor: "#94A3B8",
        tickfont: { family: FONT_MONO, size: 10.5 },
      },
      shapes: getRegimeShapes("paper"),
    };

    Plotly.react("plot-series-detail", [traceLevel, traceTrans], layout, PLOTLY_CONFIG);
  }

  /* --------------------- Figure 3: Correlation Heatmap ------------------- */
  const CORR_INSIGHTS = {
    preAITransformed:
      "<strong>Pre-AI Training Regime (Jan 2010 – Oct 2022, Stationary Transforms):</strong> Output growth (<code>INDPRO</code>) and payroll growth (<code>PAYEMS</code>) co-move strongly (<strong>r = +0.870</strong>), while job openings (<code>HWI</code>) correlate positively with both payrolls (<strong>r = +0.412</strong>) and output (<strong>r = +0.413</strong>)—confirming tightly coupled multivariate dynamics.",
    postAITransformed:
      "<strong>Post-AI &amp; Tightening Regime (Nov 2022 – Sep 2026, Stationary Transforms):</strong> Strikingly, monthly changes in job openings (<code>HWI</code>) decouple from both payroll growth (<strong>r = −0.103</strong>) and industrial production (<strong>r = −0.048</strong>), while showing near-zero contemporaneous correlation with <code>FEDFUNDS</code> (<strong>r = +0.008</strong>) and <code>BAA10YM</code> (<strong>r = −0.034</strong>).",
    transformed:
      "<strong>Full Sample Stationary Transforms (2010–2026):</strong> Across the full 16-year horizon, real activity and payrolls remain tightly linked (<strong>r = +0.851</strong>), with negative liquidity-output co-movement (<strong>r = −0.553</strong>) heavily influenced by the counter-cyclical 2020 COVID-19 monetary surge.",
    levels:
      "<strong>Raw Levels (2010–2026):</strong> Strong level correlations (e.g., <code>PAYEMS</code>–<code>M2SL</code> <strong>r = +0.870</strong>, <code>PAYEMS</code>–<code>BAA10YM</code> <strong>r = −0.809</strong>) reflect shared secular trends and long-run macro-financial linkages, underscoring the necessity of stationarity transforms (or shrinkage/cointegration modeling) in VAR/BVAR estimation.",
  };

  function renderFigure3() {
    const labels = DATA.correlations.labels;
    const z = DATA.correlations[state.fig3Mode];
    const annotations = [];

    for (let i = 0; i < labels.length; i++) {
      for (let j = 0; j < labels.length; j++) {
        const val = z[i][j];
        const absVal = Math.abs(val);
        const rawText = val >= 0 ? `+${val.toFixed(2)}` : `${val.toFixed(2)}`;
        annotations.push({
          x: labels[j],
          y: labels[i],
          text: i === j ? `<b>${rawText}</b>` : rawText,
          showarrow: false,
          font: {
            family: FONT_MONO,
            size: 12,
            color: absVal > 0.55 ? "#FFFFFF" : "#17202A",
          },
        });
      }
    }

    const trace = {
      type: "heatmap",
      z,
      x: labels,
      y: labels,
      zmin: -1,
      zmax: 1,
      colorscale: [
        [0.0, "#B5121B"],
        [0.25, "#E88D8F"],
        [0.5, "#F8FAFC"],
        [0.75, "#6B9AC4"],
        [1.0, "#102A43"],
      ],
      colorbar: {
        title: { text: "Pearson r", font: { family: FONT_MONO, size: 11 } },
        tickfont: { family: FONT_MONO, size: 10.5 },
        thickness: 14,
      },
      hovertemplate: "<b>%{y} × %{x}</b><br>r = %{z:+.3f}<extra></extra>",
    };

    const layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#FFFFFF",
      font: { family: FONT_SANS, color: "#17202A" },
      margin: { l: 88, r: 30, t: 18, b: 56 },
      xaxis: {
        side: "bottom",
        tickfont: { family: FONT_MONO, size: 11.5 },
      },
      yaxis: {
        autorange: "reversed",
        tickfont: { family: FONT_MONO, size: 11.5 },
      },
      annotations,
    };

    Plotly.react("plot-correlation", [trace], layout, PLOTLY_CONFIG);
    document.getElementById("fig3-insight").innerHTML = CORR_INSIGHTS[state.fig3Mode];
  }

  /* ----------------------- UI Control Event Bindings --------------------- */
  function bindControlGroup(groupId, onSelect) {
    const group = document.getElementById(groupId);
    if (!group) return;
    group.querySelectorAll(".seg-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        group.querySelectorAll(".seg-btn").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        onSelect(btn);
      });
    });
  }

  function bindReferenceFilters() {
    const filterGroup = document.getElementById("ref-filters");
    const items = document.querySelectorAll("#ref-list .ref-item");
    if (!filterGroup) return;

    filterGroup.querySelectorAll(".series-pill").forEach((btn) => {
      btn.addEventListener("click", () => {
        filterGroup.querySelectorAll(".series-pill").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        const tag = btn.dataset.tag;
        items.forEach((li) => {
          if (tag === "all" || li.dataset.tag === tag) {
            li.classList.remove("hidden");
          } else {
            li.classList.add("hidden");
          }
        });
      });
    });
  }

  function bindBibtexCopy() {
    const btn = document.getElementById("copy-bibtex-btn");
    const pre = document.getElementById("bibtex-content");
    if (!btn || !pre) return;

    btn.addEventListener("click", async () => {
      try {
        await navigator.clipboard.writeText(pre.textContent.trim());
        btn.textContent = "Copied!";
        setTimeout(() => {
          btn.textContent = "Copy BibTeX";
        }, 1800);
      } catch (_err) {
        const range = document.createRange();
        range.selectNodeContents(pre);
        const sel = window.getSelection();
        sel.removeAllRanges();
        sel.addRange(range);
        btn.textContent = "Selected";
        setTimeout(() => {
          btn.textContent = "Copy BibTeX";
        }, 1800);
      }
    });
  }

  function initScrollSpy() {
    const sections = Array.from(document.querySelectorAll(".article-section[id]"));
    const tocLinks = Array.from(document.querySelectorAll("#toc-list a"));
    if (!sections.length || !tocLinks.length) return;

    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const id = entry.target.id;
            tocLinks.forEach((a) => {
              a.classList.toggle("active", a.getAttribute("href") === `#${id}`);
            });
          }
        });
      },
      { rootMargin: "-20% 0px -68% 0px", threshold: 0.05 }
    );

    sections.forEach((sec) => observer.observe(sec));
  }

  function renderMath() {
    if (typeof window.renderMathInElement === "function") {
      window.renderMathInElement(document.body, {
        delimiters: [
          { left: "$$", right: "$$", display: true },
          { left: "\\(", right: "\\)", display: false },
        ],
        throwOnError: false,
      });
    }
  }

  /* ------------------------------ Boot ---------------------------------- */
  function init() {
    populateTableAndKPIs();

    bindControlGroup("fig1-scale-controls", (btn) => {
      state.fig1Scale = btn.dataset.scale;
      renderFigure1();
    });
    bindControlGroup("fig1-range-controls", (btn) => {
      state.fig1Range = btn.dataset.range;
      renderFigure1();
    });
    bindControlGroup("fig1-subset-controls", (btn) => {
      state.fig1Subset = btn.dataset.subset;
      renderFigure1();
    });
    bindControlGroup("fig3-corr-controls", (btn) => {
      state.fig3Mode = btn.dataset.mode;
      renderFigure3();
    });

    renderFigure1();
    buildFigure2Pills();
    renderFigure2();
    renderFigure3();
    bindReferenceFilters();
    bindBibtexCopy();
    initScrollSpy();

    // Render KaTeX equations once deferred scripts finish loading
    if (document.readyState === "complete") {
      renderMath();
    } else {
      window.addEventListener("load", renderMath);
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
