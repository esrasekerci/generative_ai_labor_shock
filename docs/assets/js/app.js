/* ==========================================================================
   Interactive Figures & Table Renderer
   Reads pre-compiled FRED-MD payload from `window.FRED_MACRO_DATA`
   (generated via `src/data/prepare_web_data.py`).
   ========================================================================== */

(function () {
  "use strict";

  const DATA = window.FRED_MACRO_DATA;
  if (!DATA) {
    console.error("FRED_MACRO_DATA missing. Run `python3 src/data/prepare_web_data.py`.");
    return;
  }

  const FONT_SANS = '"Source Sans 3", -apple-system, BlinkMacSystemFont, sans-serif';
  const FONT_MONO = '"JetBrains Mono", monospace';

  function isMobileViewport() {
    return window.innerWidth <= 680;
  }

  function getPlotlyConfig() {
    const mobile = isMobileViewport();
    return {
      responsive: true,
      displaylogo: false,
      displayModeBar: !mobile,
      scrollZoom: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
    };
  }

  const state = {
    fig1Scale: "zscore",
    fig1Range: "full",
    fig1Subset: "all",
    fig2Series: "HWI",
    fig3Mode: "preAITransformed",
  };

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
        fillcolor: "rgba(245, 158, 11, 0.11)",
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
        line: { color: "#B5121B", width: 1.4, dash: "dot" },
      },
    ];
  }

  /* --------------------------- Table 1 ----------------------------------- */
  function populateTable1() {
    const tbody = document.getElementById("series-table-body");
    if (!tbody) return;
    tbody.innerHTML = "";

    DATA.series.forEach((s) => {
      const st = DATA.stats[s.id];
      const rs = DATA.regimeStats[s.id];
      const tr = document.createElement("tr");
      tr.innerHTML = `
        <td>
          <span class="series-dot" style="background:${s.color};"></span>
          <code>${s.id}</code>
        </td>
        <td>
          ${s.label}<br/>
          <span style="color:var(--ink-muted); font-size:0.8rem;">${s.units}</span>
        </td>
        <td>${s.block}</td>
        <td><code>${s.tcode}</code> · ${s.transformFormula}</td>
        <td class="num">${fmtNumber(st.mean, s.id)}</td>
        <td class="num">${fmtNumber(st.sd, s.id)}</td>
        <td class="num">${fmtNumber(rs.preMean, s.id)}</td>
        <td class="num">${fmtNumber(rs.postMean, s.id)}</td>
      `;
      tbody.appendChild(tr);
    });
  }

  /* ------------------ Figure 2: Multivariate Co-Movement ----------------- */
  function renderComovement() {
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
          width: s.id === "HWI" || s.id === "PAYEMS" ? 2.3 : 1.65,
        },
        hovertemplate:
          state.fig1Scale === "zscore"
            ? `<b>${s.shortLabel}</b>: %{y:+.2f} σ<extra></extra>`
            : `<b>${s.shortLabel}</b>: %{y:.1f} (Oct 2022 = 100)<extra></extra>`,
      };
    });

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
        text: "COVID-19",
        showarrow: false,
        font: { family: FONT_MONO, size: 10.5, color: "#475569" },
        yanchor: "top",
      },
      {
        x: "2024-08-01",
        y: 0.99,
        xref: "x",
        yref: "paper",
        text: "Post-AI & Tightening",
        showarrow: false,
        font: { family: FONT_MONO, size: 10.5, color: "#92400E" },
        yanchor: "top",
      },
    ];

    const mobile = isMobileViewport();

    const layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#FFFFFF",
      font: { family: FONT_SANS, color: "#17202A", size: mobile ? 11 : 12.5 },
      margin: mobile ? { l: 42, r: 8, t: 16, b: 86 } : { l: 58, r: 18, t: 20, b: 42 },
      hovermode: mobile ? "closest" : "x unified",
      dragmode: mobile ? false : "zoom",
      xaxis: {
        range: xRange,
        showgrid: true,
        gridcolor: "#EEF2F6",
        zeroline: false,
        fixedrange: mobile,
        tickfont: { family: FONT_MONO, size: mobile ? 9.5 : 11 },
      },
      yaxis: {
        title: {
          text: mobile
            ? state.fig1Scale === "zscore"
              ? "Z-Score (σ)"
              : "Oct 2022 = 100"
            : state.fig1Scale === "zscore"
              ? "Standard Deviations from Mean (Z-Score)"
              : "Index (October 2022 = 100)",
          font: { size: mobile ? 10.5 : 12, color: "#434E5E" },
        },
        showgrid: true,
        gridcolor: "#EEF2F6",
        zeroline: false,
        fixedrange: mobile,
        tickfont: { family: FONT_MONO, size: mobile ? 9.5 : 11 },
      },
      legend: mobile
        ? {
            orientation: "h",
            yanchor: "top",
            y: -0.18,
            xanchor: "center",
            x: 0.5,
            font: { size: 10.5 },
          }
        : {
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

    Plotly.react("plot-comovement", traces, layout, getPlotlyConfig());
  }

  /* ------------------ Figure 3: Series Level & Transform ----------------- */
  function buildSeriesPills() {
    const container = document.getElementById("fig2-series-pills");
    if (!container) return;
    container.innerHTML = "";

    DATA.series.forEach((s) => {
      const btn = document.createElement("button");
      btn.className = `series-pill${s.id === state.fig2Series ? " active" : ""}`;
      btn.innerHTML = `<span class="series-dot" style="background:${s.color}; margin-right:0;"></span>${s.id}`;
      btn.addEventListener("click", () => {
        state.fig2Series = s.id;
        container.querySelectorAll(".series-pill").forEach((el) => el.classList.remove("active"));
        btn.classList.add("active");
        renderSeriesDetail();
      });
      container.appendChild(btn);
    });
  }

  function renderSeriesDetail() {
    const sid = state.fig2Series;
    const meta = DATA.series.find((s) => s.id === sid);
    const st = DATA.stats[sid];
    const dates = DATA.dates;

    const oct22Val = DATA.levels[sid][dates.indexOf("2022-10-01")];
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

    const traceLevel = {
      type: "scatter",
      mode: "lines",
      name: `${meta.id} Level`,
      x: dates,
      y: DATA.levels[sid],
      xaxis: "x",
      yaxis: "y",
      line: { color: meta.color, width: 2.2 },
      fill: "tozeroy",
      fillcolor: "rgba(16, 42, 67, 0.04)",
      hovertemplate: `<b>Level</b>: %{y:.2f}<extra></extra>`,
    };

    const traceTrans = {
      type: meta.tcode === 1 ? "scatter" : "bar",
      mode: meta.tcode === 1 ? "lines" : undefined,
      name: meta.transformFormula,
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

    const mobile = isMobileViewport();

    const layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#FFFFFF",
      font: { family: FONT_SANS, color: "#17202A", size: mobile ? 10.5 : 12 },
      margin: mobile ? { l: 48, r: 8, t: 12, b: 32 } : { l: 68, r: 18, t: 16, b: 38 },
      hovermode: mobile ? "closest" : "x unified",
      dragmode: mobile ? false : "zoom",
      showlegend: false,
      xaxis: {
        anchor: "y2",
        showgrid: true,
        gridcolor: "#EEF2F6",
        fixedrange: mobile,
        tickfont: { family: FONT_MONO, size: mobile ? 9.5 : 11 },
      },
      yaxis: {
        domain: [0.54, 1.0],
        anchor: "x",
        title: { text: mobile ? "Level" : `Level (${meta.units})`, font: { size: mobile ? 10.5 : 11.5, color: "#434E5E" } },
        showgrid: true,
        gridcolor: "#EEF2F6",
        fixedrange: mobile,
        tickfont: { family: FONT_MONO, size: mobile ? 9 : 10.5 },
      },
      yaxis2: {
        domain: [0.0, 0.42],
        anchor: "x",
        title: { text: mobile ? meta.transformFormula : meta.transformLabel, font: { size: mobile ? 10.5 : 11.5, color: "#434E5E" } },
        showgrid: true,
        gridcolor: "#EEF2F6",
        zeroline: true,
        zerolinecolor: "#94A3B8",
        fixedrange: mobile,
        tickfont: { family: FONT_MONO, size: mobile ? 9 : 10.5 },
      },
      shapes: getRegimeShapes("paper"),
    };

    Plotly.react("plot-series-detail", [traceLevel, traceTrans], layout, getPlotlyConfig());

    const captionEl = document.getElementById("fig2-caption");
    if (captionEl) {
      captionEl.innerHTML =
        `<strong>${meta.id} (${meta.label}):</strong> FRED-MD <code>tcode ${meta.tcode}</code> (${meta.transformFormula}). Change from Oct 2022 to ${fmtDateShort(st.lastDate)}: <strong>${deltaText}</strong>.`;
    }
  }

  /* ------------------ Figure 4: Correlation Heatmap ---------------------- */
  const CORR_CAPTIONS = {
    preAITransformed:
      "<strong>Pre-AI Stationary Transforms (Jan 2010 – Oct 2022):</strong> Output growth (<code>INDPRO</code>) and payroll growth (<code>PAYEMS</code>) co-move strongly (r = +0.87), while job openings (<code>HWI</code>) correlate positively with both payrolls (r = +0.41) and output (r = +0.41).",
    postAITransformed:
      "<strong>Post-AI Stationary Transforms (Nov 2022 – Sep 2026):</strong> Monthly changes in job openings (<code>HWI</code>) decouple from payroll growth (r = −0.10), output growth (r = −0.05), and contemporaneous policy/credit indicators (<code>FEDFUNDS</code> r = +0.01, <code>BAA10YM</code> r = −0.03).",
    transformed:
      "<strong>Full-Sample Stationary Transforms (2010–2026):</strong> Real activity and payroll growth maintain strong joint co-movement (r = +0.85), while the negative <code>M2SL</code>–<code>INDPRO</code> correlation (r = −0.55) reflects counter-cyclical 2020 liquidity expansion.",
    levels:
      "<strong>Raw Levels (2010–2026):</strong> High level correlations (e.g., <code>PAYEMS</code>–<code>M2SL</code> r = +0.87; <code>PAYEMS</code>–<code>BAA10YM</code> r = −0.81) reflect shared secular trends and motivate stationarity transformations prior to linear VAR estimation.",
  };

  function renderCorrelation() {
    const mobile = isMobileViewport();
    const labels = DATA.correlations.labels;
    const z = DATA.correlations[state.fig3Mode];
    const annotations = [];

    for (let i = 0; i < labels.length; i++) {
      for (let j = 0; j < labels.length; j++) {
        const val = z[i][j];
        const rawText = val >= 0 ? `+${val.toFixed(2)}` : `${val.toFixed(2)}`;
        annotations.push({
          x: labels[j],
          y: labels[i],
          text: i === j ? `<b>${rawText}</b>` : rawText,
          showarrow: false,
          font: {
            family: FONT_MONO,
            size: mobile ? 9 : 12,
            color: Math.abs(val) > 0.55 ? "#FFFFFF" : "#17202A",
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
      showscale: !mobile,
      colorbar: {
        title: { text: "r", font: { family: FONT_MONO, size: 11 } },
        tickfont: { family: FONT_MONO, size: 10.5 },
        thickness: 14,
      },
      hovertemplate: "<b>%{y} × %{x}</b><br>r = %{z:+.3f}<extra></extra>",
    };

    const layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#FFFFFF",
      font: { family: FONT_SANS, color: "#17202A" },
      margin: mobile ? { l: 62, r: 8, t: 10, b: 56 } : { l: 84, r: 24, t: 16, b: 52 },
      dragmode: false,
      xaxis: {
        side: "bottom",
        fixedrange: true,
        tickangle: mobile ? -35 : 0,
        tickfont: { family: FONT_MONO, size: mobile ? 9.5 : 11.5 },
      },
      yaxis: {
        autorange: "reversed",
        fixedrange: true,
        tickfont: { family: FONT_MONO, size: mobile ? 9.5 : 11.5 },
      },
      annotations,
    };

    Plotly.react("plot-correlation", [trace], layout, getPlotlyConfig());
    const cap = document.getElementById("fig3-caption");
    if (cap) cap.innerHTML = CORR_CAPTIONS[state.fig3Mode];
  }

  /* ----------------------- Controls & Utilities -------------------------- */
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

  function bindBibtexCopy() {
    const btn = document.getElementById("copy-bibtex-btn");
    const pre = document.getElementById("bibtex-content");
    if (!btn || !pre) return;

    btn.addEventListener("click", async () => {
      try {
        await navigator.clipboard.writeText(pre.textContent.trim());
        btn.textContent = "Copied";
        setTimeout(() => {
          btn.textContent = "Copy";
        }, 1600);
      } catch (_err) {
        const range = document.createRange();
        range.selectNodeContents(pre);
        const sel = window.getSelection();
        sel.removeAllRanges();
        sel.addRange(range);
      }
    });
  }

  function initScrollSpy() {
    const targets = Array.from(document.querySelectorAll(".article-section[id], h3[id]"));
    const tocLinks = Array.from(document.querySelectorAll("#toc-list a, .nav-links a"));
    if (!targets.length || !tocLinks.length) return;

    let activeId = "";

    function updateActive() {
      const scrollPos = window.scrollY + 110;
      let current = "";
      for (let i = 0; i < targets.length; i++) {
        const top = targets[i].offsetTop;
        if (scrollPos >= top) {
          current = targets[i].id;
        }
      }
      if (current && current !== activeId) {
        activeId = current;
        tocLinks.forEach((a) => {
          const href = a.getAttribute("href");
          if (href === `#${current}`) {
            a.classList.add("active");
          } else {
            const parentSection = document.getElementById(current)?.closest(".article-section");
            if (parentSection && href === `#${parentSection.id}` && !a.closest(".toc-sublist")) {
              a.classList.add("active");
            } else {
              a.classList.remove("active");
            }
          }
        });
      }
    }

    window.addEventListener("scroll", updateActive, { passive: true });
    setTimeout(updateActive, 100);
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

  function bindMobileNav() {
    const toggle = document.getElementById("mobile-nav-toggle");
    const drawer = document.getElementById("mobile-drawer");
    if (!toggle || !drawer) return;

    function closeDrawer() {
      drawer.classList.remove("open");
      toggle.setAttribute("aria-expanded", "false");
    }

    toggle.addEventListener("click", () => {
      const isOpen = drawer.classList.toggle("open");
      toggle.setAttribute("aria-expanded", String(isOpen));
    });

    drawer.querySelectorAll("a").forEach((link) => {
      link.addEventListener("click", closeDrawer);
    });

    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") closeDrawer();
    });
  }

  function injectMobileTableHints() {
    document.querySelectorAll(".table-wrap").forEach((wrap) => {
      if (wrap.querySelector(".mobile-swipe-hint")) return;
      const title = wrap.querySelector(".table-title");
      const hint = document.createElement("div");
      hint.className = "mobile-swipe-hint";
      hint.style.position = "sticky";
      hint.style.left = "0";
      hint.textContent = "Swipe table horizontally →";
      if (title) {
        title.insertAdjacentElement("afterend", hint);
      } else {
        wrap.prepend(hint);
      }
    });
  }

  function bindResponsivePlots() {
    let resizeTimer = null;
    let lastMobile = isMobileViewport();
    function onResize() {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => {
        const nowMobile = isMobileViewport();
        if (nowMobile !== lastMobile) {
          lastMobile = nowMobile;
          renderComovement();
          renderSeriesDetail();
          renderCorrelation();
        }
      }, 180);
    }
    window.addEventListener("resize", onResize, { passive: true });
    window.addEventListener("orientationchange", onResize, { passive: true });
  }

  function init() {
    populateTable1();
    injectMobileTableHints();
    bindMobileNav();
    bindResponsivePlots();

    bindControlGroup("fig1-scale-controls", (btn) => {
      state.fig1Scale = btn.dataset.scale;
      renderComovement();
    });
    bindControlGroup("fig1-range-controls", (btn) => {
      state.fig1Range = btn.dataset.range;
      renderComovement();
    });
    bindControlGroup("fig1-subset-controls", (btn) => {
      state.fig1Subset = btn.dataset.subset;
      renderComovement();
    });
    bindControlGroup("fig3-corr-controls", (btn) => {
      state.fig3Mode = btn.dataset.mode;
      renderCorrelation();
    });

    renderComovement();
    buildSeriesPills();
    renderSeriesDetail();
    renderCorrelation();
    bindBibtexCopy();
    initScrollSpy();

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
