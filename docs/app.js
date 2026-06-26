const MODULE_LABELS = [
  ["whole_brain", "Whole Brain"],
  ["microglia", "Microglia"],
  ["astrocyte", "Astrocyte"],
  ["oligo_opc", "Oligodendrocyte-OPC"],
  ["neuron", "Neuron"],
  ["scarlink", "SCARlink links"],
  ["reference_mapping", "Reference mapping"],
];

const state = {
  manifest: null,
  atlasCache: new Map(),
  scarlinkManifest: null,
  scarlinkCache: new Map(),
  currentModule: "whole_brain",
  currentFeatureType: "rna",
  currentFeature: null,
  colorBy: "subtype",
};

const palette = [
  "#1f77b4", "#e07a5f", "#81b29a", "#c1121f", "#6d597a", "#457b9d", "#2a9d8f",
  "#8d99ae", "#ef476f", "#bc6c25", "#264653", "#f4a261", "#7f5539", "#588157",
];

function getJSON(path) {
  return fetch(path).then((res) => {
    if (!res.ok) throw new Error(`Failed to load ${path}`);
    return res.json();
  });
}

function setText(id, text) {
  document.getElementById(id).textContent = text;
}

function escHtml(x) {
  return String(x ?? "").replace(/[&<>"']/g, (s) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[s]));
}

function renderNav() {
  const nav = document.getElementById("module-nav");
  nav.innerHTML = MODULE_LABELS.map(([key, label]) =>
    `<button data-module="${key}" class="${state.currentModule === key ? "active" : ""}">${label}</button>`
  ).join("");
  nav.querySelectorAll("button").forEach((btn) => {
    btn.addEventListener("click", () => switchModule(btn.dataset.module));
  });
}

async function loadAtlasModule(moduleKey) {
  if (state.atlasCache.has(moduleKey)) return state.atlasCache.get(moduleKey);
  const base = `data/${moduleKey}`;
  const payload = await Promise.all([
    getJSON(`${base}/cells.json`),
    getJSON(`${base}/categories.json`),
    getJSON(`${base}/metadata_summary.json`),
    getJSON(`${base}/rna_features.json`),
    getJSON(`${base}/atac_features.json`),
    getJSON(`${base}/markers.json`),
  ]);
  const out = {
    cells: payload[0],
    categories: payload[1],
    summary: payload[2],
    rna: payload[3],
    atac: payload[4],
    markers: payload[5],
  };
  state.atlasCache.set(moduleKey, out);
  return out;
}

async function loadScarlinkManifest() {
  if (!state.scarlinkManifest) {
    state.scarlinkManifest = await getJSON("data/scarlink/scarlink_manifest.json");
  }
  return state.scarlinkManifest;
}

async function loadScarlinkGene(gene) {
  if (state.scarlinkCache.has(gene)) return state.scarlinkCache.get(gene);
  const payload = await getJSON(`data/scarlink/${gene}.json`);
  state.scarlinkCache.set(gene, payload);
  return payload;
}

async function loadReference() {
  const [summary, example] = await Promise.all([
    getJSON("data/reference_mapping/summary.json"),
    getJSON("data/reference_mapping/example_mapping.json"),
  ]);
  return {summary, example};
}

function updateControlsForAtlas(data) {
  document.getElementById("atlas-controls").classList.remove("hidden");
  document.getElementById("scarlink-controls").classList.add("hidden");
  const colorBy = document.getElementById("color-by");
  colorBy.innerHTML = data.categories.available_color_by.map((key) => `<option value="${key}">${key}</option>`).join("");
  colorBy.value = state.colorBy;
  const featureType = document.getElementById("feature-type");
  featureType.value = state.currentFeatureType;
  populateFeatureSelect(data);
}

function populateFeatureSelect(data) {
  const search = document.getElementById("feature-search").value.trim().toUpperCase();
  const featureBlock = state.currentFeatureType === "atac" ? data.atac.features : data.rna.features;
  const options = Object.keys(featureBlock).filter((name) => !search || name.toUpperCase().includes(search));
  const select = document.getElementById("feature-select");
  select.innerHTML = options.map((name) => `<option value="${name}">${name}</option>`).join("");
  if (!options.length) {
    state.currentFeature = null;
    return;
  }
  if (!options.includes(state.currentFeature)) state.currentFeature = options[0];
  select.value = state.currentFeature;
}

function categoricalColors(values) {
  const levels = [...new Set(values)];
  const cmap = new Map(levels.map((level, idx) => [level, palette[idx % palette.length]]));
  return {
    colors: values.map((value) => cmap.get(value)),
    legend: levels.map((level) => ({label: level, color: cmap.get(level)})),
  };
}

function renderSummary(summary, rnaCount, atacCount) {
  const cards = [
    ["Exported cells", summary.n_exported_cells],
    ["Source cells", summary.n_total_source_cells],
    ["Subtypes", summary.n_subtypes],
    ["Diseases", summary.n_diseases],
    ["RNA features", rnaCount],
    ["ATAC features", atacCount],
  ];
  document.getElementById("summary-cards").innerHTML = cards.map(([label, value]) =>
    `<div class="summary-card"><span>${label}</span><strong>${value}</strong></div>`
  ).join("");
}

function renderMarkers(markers) {
  const panel = document.getElementById("markers-panel");
  const groups = markers.subtype_markers || {};
  panel.innerHTML = Object.entries(groups).slice(0, 18).map(([label, genes]) => `
    <div class="marker-group">
      <strong>${escHtml(label)}</strong>
      <div>${genes.map((gene) => `<span class="marker-chip" data-feature="${escHtml(gene)}">${escHtml(gene)}</span>`).join("")}</div>
    </div>
  `).join("");
  panel.querySelectorAll(".marker-chip").forEach((chip) => {
    chip.addEventListener("click", async () => {
      state.currentFeatureType = "rna";
      state.currentFeature = chip.dataset.feature;
      const data = await loadAtlasModule(state.currentModule);
      updateControlsForAtlas(data);
      renderAtlas(data);
    });
  });
}

function renderAtlas(data) {
  const cells = data.cells.cells;
  const x = cells.map((row) => row[0]);
  const y = cells.map((row) => row[1]);
  const subtype = cells.map((row) => row[2]);
  const disease = cells.map((row) => row[3]);
  const sample = cells.map((row) => row[4]);
  const dataset = cells.map((row) => row[5]);
  const pointSize = Number(document.getElementById("point-size").value);
  const opacity = Number(document.getElementById("point-opacity").value);
  const featureBlock = state.currentFeatureType === "atac" ? data.atac.features : data.rna.features;
  const selectedFeature = state.currentFeature && featureBlock[state.currentFeature] ? featureBlock[state.currentFeature].values : null;

  let color;
  let showscale = false;
  let legend = [];
  if (state.colorBy === "selected_feature" && selectedFeature) {
    color = selectedFeature;
    showscale = true;
  } else {
    const mapping = {
      subtype,
      disease,
      sample,
      dataset,
    };
    const mapped = categoricalColors(mapping[state.colorBy] || subtype);
    color = mapped.colors;
    legend = mapped.legend;
  }

  Plotly.newPlot("plot", [{
    type: "scattergl",
    mode: "markers",
    x,
    y,
    marker: {
      size: pointSize,
      opacity,
      color,
      colorscale: showscale ? "Viridis" : undefined,
      showscale,
      line: {width: 0},
      colorbar: showscale ? {title: state.currentFeature} : undefined,
    },
    text: cells.map((row, idx) =>
      `Subtype: ${escHtml(row[2])}<br>Disease: ${escHtml(row[3])}<br>Sample: ${escHtml(row[4])}<br>Dataset: ${escHtml(row[5])}${selectedFeature ? `<br>${escHtml(state.currentFeature)}: ${selectedFeature[idx]}` : ""}`
    ),
    hovertemplate: "%{text}<extra></extra>",
    showlegend: false,
  }], {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,.92)",
    margin: {l: 36, r: 20, t: 22, b: 42},
    xaxis: {title: "UMAP1", zeroline: false, gridcolor: "rgba(217,209,194,.45)"},
    yaxis: {title: "UMAP2", zeroline: false, gridcolor: "rgba(217,209,194,.45)"},
    annotations: legend.slice(0, 18).map((item, idx) => ({
      xref: "paper", yref: "paper", x: 1.01, y: 1 - idx * 0.05,
      text: `<span style="color:${item.color}">●</span> ${escHtml(item.label)}`,
      showarrow: false, align: "left",
    })),
  }, {responsive: true, displayModeBar: false});

  setText("view-title", MODULE_LABELS.find(([key]) => key === state.currentModule)?.[1] || state.currentModule);
  setText("view-subtitle", data.summary.blind_review_note);
  document.getElementById("module-note").textContent = (data.summary.warnings || []).join(" ") || "Downsampled static view with RNA and ATAC demo features.";
  renderSummary(data.summary, Object.keys(data.rna.features).length, Object.keys(data.atac.features).length);
  renderMarkers(data.markers);
}

function svgPoint(cx, cy, r, angleDeg) {
  const rad = (Math.PI / 180) * angleDeg;
  return {x: cx + Math.cos(rad) * r, y: cy + Math.sin(rad) * r};
}

function svgArc(cx, cy, r, a1, a2) {
  const p1 = svgPoint(cx, cy, r, a1);
  const p2 = svgPoint(cx, cy, r, a2);
  const large = Math.abs(a2 - a1) > 180 ? 1 : 0;
  const sweep = a2 > a1 ? 1 : 0;
  return `M ${p1.x.toFixed(2)} ${p1.y.toFixed(2)} A ${r} ${r} 0 ${large} ${sweep} ${p2.x.toFixed(2)} ${p2.y.toFixed(2)}`;
}

function mbLabel(pos) {
  return `${(pos / 1e6).toFixed(2)}Mb`;
}

function median(values) {
  const arr = values.filter((v) => Number.isFinite(v)).sort((a, b) => a - b);
  if (!arr.length) return 0;
  const mid = Math.floor(arr.length / 2);
  return arr.length % 2 ? arr[mid] : (arr[mid - 1] + arr[mid]) / 2;
}

function drawScarlinkCircle(payload) {
  const el = document.getElementById("scarlink-circle");
  const links = (payload.links || []).filter((x) =>
    x.enhancer_start !== null && x.enhancer_end !== null && x.tss !== null
  ).slice(0, 72);
  if (!links.length) {
    el.innerHTML = `<div class="note-box">No enhancer-gene links available.</div>`;
    return;
  }
  const W = 760, H = 560, cx = 330, cy = 285, R = 205, trackR = 175, enhancerR = 190, promoterR = 160, chordR = 126;
  const q = payload.query || {};
  const coords = [];
  links.forEach((l) => coords.push(l.enhancer_start, l.enhancer_end, l.tss));
  let minPos = Math.min(...coords), maxPos = Math.max(...coords);
  const span0 = Math.max(1, maxPos - minPos);
  const pad = Math.max(25000, span0 * 0.12);
  minPos = Math.max(1, Math.floor(minPos - pad));
  maxPos = Math.ceil(maxPos + pad);
  const span = Math.max(1, maxPos - minPos);
  const startAngle = -225, endAngle = 135;
  const angleFor = (pos) => startAngle + ((pos - minPos) / span) * (endAngle - startAngle);
  const chr = q.chr || links[0].enhancer_chr || "";
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="100%"><rect width="${W}" height="${H}" fill="transparent"/><circle cx="${cx}" cy="${cy}" r="${R + 28}" fill="rgba(255,255,255,.72)"/><text x="${cx}" y="34" text-anchor="middle" font-size="16" font-weight="700">${escHtml(chr)} SCARlink enhancer-promoter map</text><text x="${cx}" y="56" text-anchor="middle" font-size="12" fill="#6b6f76">${escHtml(chr)}:${minPos.toLocaleString()}-${maxPos.toLocaleString()}</text>`;
  svg += `<path d="${svgArc(cx, cy, R, startAngle, endAngle)}" fill="none" stroke="#111" stroke-width="2"/>`;
  for (let i = 0; i <= 18; i++) {
    const pos = minPos + span * i / 18;
    const a = angleFor(pos);
    const p1 = svgPoint(cx, cy, R - 3, a);
    const p2 = svgPoint(cx, cy, R + (i % 3 === 0 ? 13 : 8), a);
    const lab = svgPoint(cx, cy, R + 31, a);
    svg += `<line x1="${p1.x}" y1="${p1.y}" x2="${p2.x}" y2="${p2.y}" stroke="#222" stroke-width="${i % 3 === 0 ? 1.6 : 1}"/>`;
    if (i % 2 === 0) svg += `<text x="${lab.x}" y="${lab.y}" font-size="11" text-anchor="middle" transform="rotate(${a < 0 ? a + 90 : a - 90} ${lab.x} ${lab.y})">${mbLabel(pos)}</text>`;
  }
  for (let i = 0; i < 10; i++) {
    if (i % 3 === 2) continue;
    const a1 = startAngle + i * (endAngle - startAngle) / 10 + 2;
    const a2 = startAngle + (i + 0.72) * (endAngle - startAngle) / 10;
    svg += `<path d="${svgArc(cx, cy, trackR, a1, a2)}" fill="none" stroke="#1f77b4" stroke-width="15"/>`;
  }
  if (Number.isFinite(q.start) && Number.isFinite(q.end)) {
    svg += `<path d="${svgArc(cx, cy, R + 5, angleFor(q.start), angleFor(q.end))}" fill="none" stroke="#264bff" stroke-width="4"/>`;
  }
  links.slice(0, 46).forEach((l, idx) => {
    const ea = angleFor((l.enhancer_start + l.enhancer_end) / 2);
    const pa = angleFor(l.tss);
    const p1 = svgPoint(cx, cy, chordR, ea);
    const p2 = svgPoint(cx, cy, chordR, pa);
    const c1 = svgPoint(cx, cy, 42, ea);
    const c2 = svgPoint(cx, cy, 42, pa);
    const width = idx === 0 ? 5.5 : Math.max(1.1, Math.min(4.2, 1.3 + Math.abs(Number(l.regression_coef || 0)) * 18000));
    const color = idx === 0 ? "#ffea00" : (l.effect === "repression" ? "#187d8a" : "#d62728");
    svg += `<path d="M ${p1.x.toFixed(2)} ${p1.y.toFixed(2)} C ${c1.x.toFixed(2)} ${c1.y.toFixed(2)}, ${c2.x.toFixed(2)} ${c2.y.toFixed(2)}, ${p2.x.toFixed(2)} ${p2.y.toFixed(2)}" fill="none" stroke="${color}" stroke-width="${width}" stroke-opacity=".78"/>`;
  });
  links.slice(0, 40).forEach((l) => {
    const ea = angleFor((l.enhancer_start + l.enhancer_end) / 2);
    const pa = angleFor(l.tss);
    const ep = svgPoint(cx, cy, enhancerR, ea);
    const pp = svgPoint(cx, cy, promoterR, pa);
    svg += `<circle cx="${ep.x}" cy="${ep.y}" r="2.4" fill="#666"/><circle cx="${pp.x}" cy="${pp.y}" r="3" fill="#2ca02c"/>`;
  });
  svg += `</svg>`;
  el.innerHTML = svg;
}

function drawScarlinkBoxplot(payload) {
  const groups = (payload.boxplot?.groups || []).filter((g) => Array.isArray(g.values) && g.values.length);
  const traces = groups.map((g, idx) => ({
    type: "box",
    y: g.values,
    name: g.label,
    boxpoints: "all",
    jitter: 0.35,
    pointpos: 0,
    marker: {size: 4, opacity: 0.48, color: palette[idx % palette.length]},
    line: {width: 1.2},
    hovertemplate: `${g.label}<br>z-score=%{y:.4f}<extra></extra>`,
  }));
  Plotly.newPlot("scarlink-boxplot", traces, {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,.92)",
    title: {text: "SCARlink z-score by celltype_r2", font: {size: 15}},
    margin: {l: 58, r: 20, t: 42, b: 110},
    xaxis: {tickangle: -38, automargin: true},
    yaxis: {title: "z-score", zeroline: true, gridcolor: "rgba(217,208,194,.55)"},
    showlegend: false,
  }, {responsive: true, displayModeBar: false});
}

async function renderScarlink() {
  document.getElementById("atlas-controls").classList.add("hidden");
  document.getElementById("scarlink-controls").classList.remove("hidden");
  document.getElementById("plot").classList.add("hidden");
  document.getElementById("scarlink-layout").classList.remove("hidden");
  document.getElementById("reference-layout").classList.add("hidden");
  const manifest = await loadScarlinkManifest();
  const select = document.getElementById("scarlink-gene");
  select.innerHTML = manifest.genes.map((x) => `<option value="${x.gene}">${x.gene}</option>`).join("");
  const gene = select.value || manifest.genes[0]?.gene;
  const payload = await loadScarlinkGene(gene);
  setText("view-title", "SCARlink links");
  setText("view-subtitle", payload.message || "Static SCARlink examples.");
  document.getElementById("module-note").textContent = "Example enhancer-gene links exported from SCARlink result tables.";
  drawScarlinkCircle(payload);
  drawScarlinkBoxplot(payload);
  document.getElementById("scarlink-head").innerHTML = "<tr><th>Disease</th><th>Cell type</th><th>Peak</th><th>Coef</th><th>FDR</th></tr>";
  document.getElementById("scarlink-table").innerHTML = payload.table.slice(0, 120).map((row) =>
    `<tr><td>${escHtml(row.disease)}</td><td>${escHtml(row.celltype_r2)}</td><td>${escHtml(row.peak)}</td><td>${row.regression_coef}</td><td>${row.fdr}</td></tr>`
  ).join("");
  document.getElementById("scarlink-examples").innerHTML = manifest.genes.map((row) =>
    `<div class="example-card" data-gene="${row.gene}"><strong>${row.gene}</strong><div>${row.n_rows} links</div><small>${escHtml(row.query_region)}</small></div>`
  ).join("");
  document.querySelectorAll(".example-card").forEach((card) => {
    card.addEventListener("click", async () => {
      document.getElementById("scarlink-gene").value = card.dataset.gene;
      const next = await loadScarlinkGene(card.dataset.gene);
      drawScarlinkCircle(next);
      drawScarlinkBoxplot(next);
      setText("view-subtitle", next.message || "Static SCARlink examples.");
      document.getElementById("scarlink-table").innerHTML = next.table.slice(0, 120).map((row) =>
        `<tr><td>${escHtml(row.disease)}</td><td>${escHtml(row.celltype_r2)}</td><td>${escHtml(row.peak)}</td><td>${row.regression_coef}</td><td>${row.fdr}</td></tr>`
      ).join("");
    });
  });
  renderSummary({n_exported_cells: manifest.genes.length, n_total_source_cells: "Multi-disease", n_subtypes: "celltype_r2", n_diseases: "SCARlink",}, "Example", "Links");
  document.getElementById("markers-panel").innerHTML = `<div class="marker-group">Load one of the example genes to inspect circle, boxplot, and table views.</div>`;
}

async function renderReference() {
  document.getElementById("atlas-controls").classList.add("hidden");
  document.getElementById("scarlink-controls").classList.add("hidden");
  document.getElementById("plot").classList.add("hidden");
  document.getElementById("scarlink-layout").classList.add("hidden");
  document.getElementById("reference-layout").classList.remove("hidden");
  const {summary, example} = await loadReference();
  setText("view-title", "Reference mapping");
  setText("view-subtitle", summary.description);
  document.getElementById("module-note").textContent = "Static workflow-only page for blind-review deployment.";
  document.getElementById("reference-layout").innerHTML = `
    <div class="panel inner-panel">
      <h3>Workflow</h3>
      <ol>${summary.workflow.map((step) => `<li>${escHtml(step)}</li>`).join("")}</ol>
    </div>
    <div class="panel inner-panel">
      <h3>Example mapping table</h3>
      <div class="table-wrap">
        <table>
          <thead><tr>${example.columns.map((c) => `<th>${escHtml(c)}</th>`).join("")}</tr></thead>
          <tbody>${example.rows.map((row) => `<tr>${row.map((v) => `<td>${escHtml(v)}</td>`).join("")}</tr>`).join("")}</tbody>
        </table>
      </div>
    </div>`;
  renderSummary({n_exported_cells: summary.modules.length, n_total_source_cells: "Reference", n_subtypes: "Static", n_diseases: "Static"}, "Workflow", "Table");
  document.getElementById("markers-panel").innerHTML = `<div class="marker-group">Reference mapping is shown as a compact static summary in this demo.</div>`;
}

async function switchModule(moduleKey) {
  state.currentModule = moduleKey;
  renderNav();
  if (moduleKey === "scarlink") {
    return renderScarlink();
  }
  if (moduleKey === "reference_mapping") {
    return renderReference();
  }
  document.getElementById("plot").classList.remove("hidden");
  document.getElementById("scarlink-layout").classList.add("hidden");
  document.getElementById("reference-layout").classList.add("hidden");
  const data = await loadAtlasModule(moduleKey);
  state.colorBy = "subtype";
  state.currentFeatureType = "rna";
  state.currentFeature = Object.keys(data.rna.features)[0] || null;
  updateControlsForAtlas(data);
  renderAtlas(data);
}

async function init() {
  state.manifest = await getJSON("data/manifest.json");
  renderNav();
  document.getElementById("color-by").addEventListener("change", async (e) => {
    state.colorBy = e.target.value;
    const data = await loadAtlasModule(state.currentModule);
    renderAtlas(data);
  });
  document.getElementById("feature-type").addEventListener("change", async (e) => {
    state.currentFeatureType = e.target.value;
    const data = await loadAtlasModule(state.currentModule);
    populateFeatureSelect(data);
    renderAtlas(data);
  });
  document.getElementById("feature-select").addEventListener("change", async (e) => {
    state.currentFeature = e.target.value;
    const data = await loadAtlasModule(state.currentModule);
    renderAtlas(data);
  });
  document.getElementById("feature-search").addEventListener("input", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    populateFeatureSelect(data);
  });
  document.getElementById("point-size").addEventListener("input", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    renderAtlas(data);
  });
  document.getElementById("point-opacity").addEventListener("input", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    renderAtlas(data);
  });
  document.getElementById("reset-view").addEventListener("click", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    renderAtlas(data);
  });
  document.getElementById("scarlink-load").addEventListener("click", async () => {
    const gene = document.getElementById("scarlink-gene").value;
    const payload = await loadScarlinkGene(gene);
    drawScarlinkCircle(payload);
    drawScarlinkBoxplot(payload);
    setText("view-subtitle", payload.message || "Static SCARlink examples.");
  });
  switchModule("whole_brain");
}

init().catch((err) => {
  setText("view-title", "Load error");
  setText("view-subtitle", err.message);
});
