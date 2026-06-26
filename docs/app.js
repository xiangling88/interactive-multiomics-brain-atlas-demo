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
  moduleCache: new Map(),
  featureCache: new Map(),
  scarlinkManifest: null,
  scarlinkCache: new Map(),
  currentModule: "whole_brain",
  currentFeatureType: "rna",
  currentFeature: null,
  currentColorBy: "subtype",
  currentScarlinkGene: null,
  currentScarlinkDisease: null,
};

const palette = [
  "#8f2d2a", "#d16f5b", "#bfa239", "#3c7d67", "#4f8797", "#7a6eb4", "#9f5378",
  "#cb8b2f", "#4d5a68", "#87915b", "#ba5b44", "#6b89c6", "#ad6b92", "#597b80", "#9a8f7a",
];

function getJSON(path) {
  return fetch(path).then((res) => {
    if (!res.ok) throw new Error(`Failed to load ${path}`);
    return res.json();
  });
}

function escHtml(x) {
  return String(x ?? "").replace(/[&<>"']/g, (s) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[s]));
}

function setText(id, text) {
  const el = document.getElementById(id);
  if (el) el.textContent = text;
}

function renderNav() {
  const nav = document.getElementById("module-nav");
  nav.innerHTML = MODULE_LABELS.map(([key, label]) => `<button data-module="${key}" class="${state.currentModule === key ? "active" : ""}">${label}</button>`).join("");
  nav.querySelectorAll("button").forEach((btn) => btn.addEventListener("click", () => switchModule(btn.dataset.module)));
}

function moduleEntry(moduleKey) {
  return (state.manifest?.modules || []).find((m) => m.module === moduleKey);
}

function renderSummaryCards(summary, rnaCount, atacCount) {
  const cards = [
    ["Exported cells", summary.n_exported_cells],
    ["Source cells", summary.n_total_source_cells],
    ["Subtypes", summary.n_subtypes ?? "-"],
    ["Diseases", summary.n_diseases ?? "-"],
    ["RNA features", rnaCount],
    ["ATAC features", atacCount],
  ];
  document.getElementById("summary-cards").innerHTML = cards.map(([label, value]) => `<div class="summary-card"><span>${label}</span><strong>${value}</strong></div>`).join("");
}

function categoricalColors(values) {
  const levels = [...new Set(values)];
  const cmap = new Map(levels.map((level, idx) => [level, palette[idx % palette.length]]));
  return {colors: values.map((value) => cmap.get(value)), cmap, levels};
}

function emptyPlot(divId, title, note) {
  Plotly.react(divId, [{x:[0], y:[0], mode:"markers", marker:{size:0, opacity:0}}], {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "white",
    margin: {l: 42, r: 18, t: 50, b: 48},
    xaxis: {visible: false},
    yaxis: {visible: false},
    annotations: [{text: `${escHtml(title)}<br><span style="font-size:12px;color:#667085">${escHtml(note)}</span>`, x: .5, y: .55, xref: "paper", yref: "paper", showarrow: false}],
  }, {responsive: true, displayModeBar: false});
}

async function loadAtlasModule(moduleKey) {
  if (state.moduleCache.has(moduleKey)) return state.moduleCache.get(moduleKey);
  const entry = moduleEntry(moduleKey);
  const base = `data/${moduleKey}`;
  const [categories, summary, markers] = await Promise.all([
    getJSON(`${base}/categories.json`),
    getJSON(`${base}/metadata_summary.json`),
    getJSON(`${base}/markers.json`),
  ]);
  const cells = [];
  const parts = entry.cell_parts || summary.cell_parts || [];
  for (let i = 0; i < parts.length; i++) {
    setText("load-progress", `Loading cells ${i + 1}/${parts.length}`);
    const part = await getJSON(parts[i]);
    cells.push(...part.cells);
  }
  const payload = {entry, categories, summary, markers, cells, columns: ["x", "y", ...(entry.display_fields || ["subtype", "disease", "sample", "RL6", "RL_4", "RL3_2", "RL3_1", "RL_3", "RL_2"]).slice(2)]};
  state.moduleCache.set(moduleKey, payload);
  setText("load-progress", `${cells.length.toLocaleString()} cells loaded`);
  return payload;
}

async function loadFeature(moduleKey, modality, feature) {
  const cacheKey = `${moduleKey}::${modality}::${feature}`;
  if (state.featureCache.has(cacheKey)) return state.featureCache.get(cacheKey);
  const entry = moduleEntry(moduleKey);
  const path = entry?.feature_files?.[modality]?.[feature];
  if (!path) return null;
  setText("load-progress", `Loading ${feature}`);
  const payload = await getJSON(path);
  state.featureCache.set(cacheKey, payload);
  setText("load-progress", `${moduleEntry(moduleKey).n_exported_cells.toLocaleString()} cells loaded`);
  return payload;
}

function decodeFeatureValues(feature) {
  if (!feature) return null;
  if (feature.encoding === "dense") return feature.values;
  if (feature.encoding === "quantized") {
    const scale = Number(feature.scale_max || feature.q99 || 1);
    return feature.values.map((v) => (Number(v) / 255) * scale);
  }
  if (feature.encoding === "sparse") {
    const out = new Array(feature.length).fill(0);
    feature.indices.forEach((idx, i) => { out[idx] = feature.values[i]; });
    return out;
  }
  return feature.values || null;
}

function configureAtlasControls(data) {
  document.getElementById("atlas-controls").classList.remove("hidden");
  document.getElementById("scarlink-controls").classList.add("hidden");
  document.getElementById("atlas-layout").classList.remove("hidden");
  document.getElementById("scarlink-layout").classList.add("hidden");
  document.getElementById("reference-layout").classList.add("hidden");
  const colors = data.categories.available_color_by || ["subtype", "disease", "sample"];
  const colorBy = document.getElementById("color-by");
  colorBy.innerHTML = colors.map((x) => `<option value="${x}">${x}</option>`).join("");
  if (!colors.includes(state.currentColorBy)) state.currentColorBy = colors[0];
  colorBy.value = state.currentColorBy;
  populateFeatureSelect(data);
}

function populateFeatureSelect(data) {
  const featureType = state.currentFeatureType;
  const search = document.getElementById("feature-search").value.trim().toUpperCase();
  const names = (data.entry.features?.[featureType] || []).filter((name) => !search || name.toUpperCase().includes(search));
  const select = document.getElementById("feature-select");
  select.innerHTML = `<option value="">None</option>` + names.map((x) => `<option value="${x}">${x}</option>`).join("");
  if (state.currentFeature && !names.includes(state.currentFeature)) state.currentFeature = null;
  select.value = state.currentFeature || "";
}

function renderMarkers(markers) {
  const groups = markers.subtype_markers || {};
  document.getElementById("markers-panel").innerHTML = Object.entries(groups).slice(0, 22).map(([label, genes]) => `
    <div class="marker-group">
      <strong>${escHtml(label)}</strong>
      <div>${genes.map((gene) => `<span class="marker-chip" data-feature="${escHtml(gene)}">${escHtml(gene)}</span>`).join("")}</div>
    </div>
  `).join("");
  document.querySelectorAll(".marker-chip").forEach((chip) => chip.addEventListener("click", async () => {
    state.currentFeatureType = "rna";
    state.currentFeature = chip.dataset.feature;
    const data = await loadAtlasModule(state.currentModule);
    populateFeatureSelect(data);
    await renderAtlasModule(data);
  }));
}

async function renderAtlasModule(data) {
  setText("view-title", MODULE_LABELS.find(([k]) => k === state.currentModule)?.[1] || data.entry.label);
  setText("view-subtitle", `${data.entry.n_exported_cells.toLocaleString()} exported cells from ${data.entry.n_total_source_cells.toLocaleString()} source cells.`);
  document.getElementById("module-note").textContent = (data.summary.warnings || []).join(" ") || "Chunked UMAP cells and lazy-loaded selected features.";
  renderSummaryCards(data.summary, (data.entry.features?.rna || []).length, (data.entry.features?.atac || []).length);
  renderMarkers(data.markers);

  const cells = data.cells;
  const x = cells.map((row) => row[0]);
  const y = cells.map((row) => row[1]);
  const fieldIndex = {subtype: 2, disease: 3, sample: 4, RL6: 5, RL_4: 6, RL3_2: 7, RL3_1: 8, RL_3: 9, RL_2: 10};
  const pointSize = Number(document.getElementById("point-size").value);
  const opacity = Number(document.getElementById("point-opacity").value);
  let markerColor;
  let showscale = false;
  let hoverValues = null;
  let colorLegend = null;
  let detailFeature = null;

  if (state.currentColorBy === "selected_feature" && state.currentFeature) {
    detailFeature = await loadFeature(state.currentModule, state.currentFeatureType, state.currentFeature);
    hoverValues = decodeFeatureValues(detailFeature);
    markerColor = hoverValues;
    showscale = true;
  } else {
    const vals = cells.map((row) => row[fieldIndex[state.currentColorBy] ?? 2]);
    colorLegend = categoricalColors(vals);
    markerColor = colorLegend.colors;
  }

  const hoverText = cells.map((row, idx) => {
    const parts = [`Subtype: ${escHtml(row[2])}`, `Disease: ${escHtml(row[3])}`, `Sample: ${escHtml(row[4])}`];
    ["RL6", "RL_4", "RL3_2", "RL3_1", "RL_3", "RL_2"].forEach((field, offset) => {
      const v = row[5 + offset];
      if (v && v !== "NA") parts.push(`${field}: ${escHtml(v)}`);
    });
    if (hoverValues) parts.push(`${escHtml(state.currentFeature)}: ${Number(hoverValues[idx]).toFixed(4)}`);
    return parts.join("<br>");
  });

  Plotly.react("umap-plot", [{
    x, y, type: "scattergl", mode: "markers",
    marker: {
      size: pointSize,
      opacity,
      color: markerColor,
      colorscale: showscale ? (state.currentFeatureType === "rna" ? "Reds" : "Tealgrn") : undefined,
      colorbar: showscale ? {title: state.currentFeature} : undefined,
      showscale,
      line: {width: 0},
    },
    text: hoverText,
    hovertemplate: "%{text}<extra></extra>",
    showlegend: false,
  }], {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "white",
    autosize: true,
    margin: {l: 48, r: 22, t: 22, b: 44},
    xaxis: {title: "UMAP 1", zeroline: false, showgrid: false},
    yaxis: {title: "UMAP 2", zeroline: false, showgrid: false, scaleanchor: "x", scaleratio: 1},
    annotations: colorLegend ? colorLegend.levels.slice(0, 18).map((label, idx) => ({xref: "paper", yref: "paper", x: 1.02, y: 1 - idx * 0.045, text: `<span style="color:${colorLegend.cmap.get(label)}">●</span> ${escHtml(label)}`, showarrow: false, align: "left"})) : [],
  }, {responsive: true, displayModeBar: false});

  setText("umap-caption", state.currentColorBy === "selected_feature" ? `${state.currentFeatureType.toUpperCase()} feature overlay` : `${state.currentColorBy} categories`);

  if (state.currentFeature) {
    const feature = detailFeature || await loadFeature(state.currentModule, state.currentFeatureType, state.currentFeature);
    renderViolin(feature, data.categories.subtype_label || "Subtype");
  } else {
    renderRelationBar(data);
  }
}

function renderRelationBar(data) {
  setText("detail-title", "Subtype x disease relationship");
  setText("detail-caption", "No feature selected. Showing exported subtype counts.");
  const counts = Object.entries(data.summary.source_subtype_counts || {}).slice(0, 18);
  Plotly.react("detail-plot", [{type: "bar", x: counts.map((x) => x[0]), y: counts.map((x) => x[1]), marker: {color: counts.map((_, i) => palette[i % palette.length])}}], {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "white",
    margin: {l: 50, r: 20, t: 28, b: 150},
    xaxis: {tickangle: -35, automargin: true},
    yaxis: {title: "Cells", gridcolor: "rgba(217,224,232,.65)"},
    showlegend: false,
  }, {responsive: true, displayModeBar: false});
}

function renderViolin(feature, subtypeLabel) {
  const groupBy = document.getElementById("violin-group-by").value;
  const payload = feature?.violin?.[groupBy];
  if (!payload || !payload.traces?.length) {
    emptyPlot("detail-plot", "No violin data", "Selected feature has no grouped values.");
    return;
  }
  setText("detail-title", `${feature.type.toUpperCase()} violin`);
  setText("detail-caption", `${feature.label} grouped by ${groupBy}`);
  const traces = payload.traces.slice(0, 32).map((row, idx) => ({
    type: "violin",
    y: row.sample,
    name: row.name,
    box: {visible: true},
    meanline: {visible: true},
    points: row.sample.length <= 80 ? "all" : false,
    marker: {size: 3, opacity: .35, color: palette[idx % palette.length]},
    line: {color: palette[idx % palette.length]},
  }));
  Plotly.react("detail-plot", traces, {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "white",
    margin: {l: 58, r: 20, t: 24, b: 150},
    xaxis: {tickangle: -35, automargin: true},
    yaxis: {title: feature.type === "rna" ? "RNA expression" : "ATAC accessibility", gridcolor: "rgba(217,224,232,.65)"},
    showlegend: false,
  }, {responsive: true, displayModeBar: false});
}

async function loadScarlinkManifest() {
  if (!state.scarlinkManifest) state.scarlinkManifest = await getJSON("data/scarlink/scarlink_manifest.json");
  return state.scarlinkManifest;
}

async function loadScarlinkPayload(gene, diseaseSlug) {
  const cacheKey = `${gene}::${diseaseSlug}`;
  if (state.scarlinkCache.has(cacheKey)) return state.scarlinkCache.get(cacheKey);
  const payload = await getJSON(`data/scarlink/${gene}/${diseaseSlug}.json`);
  state.scarlinkCache.set(cacheKey, payload);
  return payload;
}

function svgPoint(cx, cy, r, angleDeg) {
  const rad = (Math.PI / 180) * angleDeg;
  return {x: cx + Math.cos(rad) * r, y: cy + Math.sin(rad) * r};
}

function svgArc(cx, cy, r, a1, a2) {
  const p1 = svgPoint(cx, cy, r, a1);
  const p2 = svgPoint(cx, cy, r, a2);
  const large = Math.abs(a2 - a1) > 180 ? 1 : 0;
  return `M ${p1.x.toFixed(2)} ${p1.y.toFixed(2)} A ${r} ${r} 0 ${large} 1 ${p2.x.toFixed(2)} ${p2.y.toFixed(2)}`;
}

function drawScarlinkCircle(data) {
  const el = document.getElementById("scarlink-circle");
  const links = (data.links || []).slice(0, 72);
  if (!links.length) {
    el.innerHTML = `<div class="section-title"><p>No enhancer-gene links in the selected disease layer.</p></div>`;
    return;
  }
  const W = 760, H = 560, cx = 330, cy = 285, outer = 205, inner = 128;
  const coords = [];
  links.forEach((l) => coords.push(l.enhancer_start, l.enhancer_end, l.tss));
  let minPos = Math.min(...coords), maxPos = Math.max(...coords);
  const pad = Math.max(25000, (maxPos - minPos) * .12);
  minPos -= pad; maxPos += pad;
  const span = Math.max(1, maxPos - minPos);
  const angle = (pos) => -225 + ((pos - minPos) / span) * 360;
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="100%"><rect width="${W}" height="${H}" fill="white"/><circle cx="${cx}" cy="${cy}" r="${outer + 20}" fill="rgba(255,255,255,.78)"/><text x="${cx}" y="36" text-anchor="middle" font-size="17" font-weight="700">${escHtml(data.gene)} | ${escHtml(data.disease)}</text><text x="${cx}" y="58" text-anchor="middle" font-size="12" fill="#667085">${escHtml(data.query.chr)}:${Math.round(minPos).toLocaleString()}-${Math.round(maxPos).toLocaleString()}</text>`;
  svg += `<path d="${svgArc(cx, cy, outer, -225, 135)}" fill="none" stroke="#173f5f" stroke-width="2.2"/>`;
  for (let i = 0; i <= 18; i++) {
    const pos = minPos + (span * i / 18);
    const a = angle(pos);
    const p1 = svgPoint(cx, cy, outer - 4, a);
    const p2 = svgPoint(cx, cy, outer + (i % 3 === 0 ? 12 : 7), a);
    svg += `<line x1="${p1.x}" y1="${p1.y}" x2="${p2.x}" y2="${p2.y}" stroke="#173f5f" stroke-width="${i % 3 === 0 ? 1.4 : 1}"/>`;
  }
  links.forEach((l, i) => {
    const ea = angle((l.enhancer_start + l.enhancer_end) / 2);
    const pa = angle(l.tss);
    const p1 = svgPoint(cx, cy, inner, ea);
    const p2 = svgPoint(cx, cy, inner, pa);
    const c1 = svgPoint(cx, cy, 36, ea);
    const c2 = svgPoint(cx, cy, 36, pa);
    const color = l.effect === "repression" ? "#0F766E" : "#B42318";
    const width = i < 5 ? 3.8 : 1.6;
    const opacity = i < 20 ? 0.78 : 0.35;
    svg += `<path d="M ${p1.x.toFixed(2)} ${p1.y.toFixed(2)} C ${c1.x.toFixed(2)} ${c1.y.toFixed(2)}, ${c2.x.toFixed(2)} ${c2.y.toFixed(2)}, ${p2.x.toFixed(2)} ${p2.y.toFixed(2)}" fill="none" stroke="${color}" stroke-width="${width}" opacity="${opacity}"/>`;
  });
  svg += `</svg>`;
  el.innerHTML = svg;
}

function drawScarlinkBoxplot(data) {
  const traces = Object.entries(data.box || {}).slice(0, 28).map(([name, values], idx) => ({
    type: "box", y: values, name, boxpoints: "all", jitter: .35, marker: {size: 4, opacity: .45, color: palette[idx % palette.length]},
  }));
  Plotly.react("scarlink-boxplot", traces, {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "white",
    margin: {l: 55, r: 20, t: 22, b: 110},
    yaxis: {title: "z-score", gridcolor: "rgba(217,224,232,.65)"},
    xaxis: {tickangle: -35},
    showlegend: false,
  }, {responsive: true, displayModeBar: false});
}

function renderScarlinkTable(data) {
  const filter = document.getElementById("scarlink-celltype-filter").value.trim().toLowerCase();
  const rows = (data.table || []).filter((row) => !filter || String(row.celltype_r2).toLowerCase().includes(filter));
  document.getElementById("scarlink-head").innerHTML = "<tr><th>Disease</th><th>Cell type</th><th>Peak</th><th>Coef</th><th>FDR</th><th>z-score</th></tr>";
  document.getElementById("scarlink-table").innerHTML = rows.slice(0, 160).map((row) => `<tr><td>${escHtml(row.disease)}</td><td>${escHtml(row.celltype_r2)}</td><td>${escHtml(row.peak)}</td><td>${row.regression_coef}</td><td>${row.fdr}</td><td>${row.z_score}</td></tr>`).join("");
}

async function renderScarlink() {
  const manifest = await loadScarlinkManifest();
  document.getElementById("atlas-controls").classList.add("hidden");
  document.getElementById("scarlink-controls").classList.remove("hidden");
  document.getElementById("atlas-layout").classList.add("hidden");
  document.getElementById("scarlink-layout").classList.remove("hidden");
  document.getElementById("reference-layout").classList.add("hidden");
  document.getElementById("markers-panel").innerHTML = `<div class="marker-group">SCARlink examples are organized by gene and disease. Use the disease selector to redraw the circle plot.</div>`;
  renderSummaryCards({n_exported_cells: manifest.genes.length, n_total_source_cells: manifest.diseases.length, n_subtypes: "SCARlink", n_diseases: manifest.diseases.length}, "Example", "Disease layers");
  setText("view-title", "SCARlink links");
  setText("view-subtitle", "Disease-aware static SCARlink view.");
  document.getElementById("module-note").textContent = "Circle plot, boxplot, and table are all drawn from static JSON without any backend API.";
  const geneSelect = document.getElementById("scarlink-gene");
  if (!geneSelect.options.length) geneSelect.innerHTML = manifest.genes.map((g) => `<option value="${g.gene}">${g.gene}</option>`).join("");
  const gene = geneSelect.value || manifest.genes[0]?.gene;
  const geneEntry = manifest.genes.find((g) => g.gene === gene) || manifest.genes[0];
  const diseaseSelect = document.getElementById("scarlink-disease");
  diseaseSelect.innerHTML = geneEntry.diseases.map((d) => `<option value="${d.slug}">${d.name}</option>`).join("");
  const diseaseSlug = diseaseSelect.value || geneEntry.diseases[0]?.slug;
  const payload = await loadScarlinkPayload(gene, diseaseSlug);
  state.currentScarlinkGene = gene;
  state.currentScarlinkDisease = diseaseSlug;
  setText("scarlink-caption", `${gene} in ${payload.disease}`);
  drawScarlinkCircle(payload);
  drawScarlinkBoxplot(payload);
  renderScarlinkTable(payload);
}

async function loadReference() {
  const [summary, example] = await Promise.all([getJSON("data/reference_mapping/summary.json"), getJSON("data/reference_mapping/example_mapping.json")]);
  return {summary, example};
}

async function renderReference() {
  document.getElementById("atlas-controls").classList.add("hidden");
  document.getElementById("scarlink-controls").classList.add("hidden");
  document.getElementById("atlas-layout").classList.add("hidden");
  document.getElementById("scarlink-layout").classList.add("hidden");
  document.getElementById("reference-layout").classList.remove("hidden");
  const {summary, example} = await loadReference();
  setText("view-title", "Reference mapping");
  setText("view-subtitle", summary.description);
  document.getElementById("module-note").textContent = "Static summary panel for reference mapping workflow.";
  document.getElementById("markers-panel").innerHTML = `<div class="marker-group">Reference mapping is shown as a static workflow summary in this GitHub Pages deployment.</div>`;
  renderSummaryCards({n_exported_cells: summary.modules.length, n_total_source_cells: summary.modules.length, n_subtypes: "Workflow", n_diseases: "Static"}, "Summary", "Table");
  document.getElementById("reference-layout").innerHTML = `
    <div class="plot-card reference-panel">
      <div class="plot-card-head"><h3>Workflow</h3><span>Static summary</span></div>
      <ol>${summary.workflow.map((step) => `<li>${escHtml(step)}</li>`).join("")}</ol>
    </div>
    <div class="plot-card reference-panel" style="margin-top:1rem;">
      <div class="plot-card-head"><h3>Example mapping table</h3><span>Exported modules</span></div>
      <div class="table-wrap"><table><thead><tr>${example.columns.map((c) => `<th>${escHtml(c)}</th>`).join("")}</tr></thead><tbody>${example.rows.map((row) => `<tr>${row.map((x) => `<td>${escHtml(x)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>
    </div>`;
}

async function switchModule(moduleKey) {
  state.currentModule = moduleKey;
  renderNav();
  if (moduleKey === "scarlink") return renderScarlink();
  if (moduleKey === "reference_mapping") return renderReference();
  const data = await loadAtlasModule(moduleKey);
  if (!(data.categories.available_color_by || []).includes(state.currentColorBy)) state.currentColorBy = "subtype";
  state.currentFeature = null;
  configureAtlasControls(data);
  await renderAtlasModule(data);
}

async function init() {
  state.manifest = await getJSON("data/manifest.json");
  renderNav();
  document.getElementById("color-by").addEventListener("change", async (e) => {
    state.currentColorBy = e.target.value;
    const data = await loadAtlasModule(state.currentModule);
    await renderAtlasModule(data);
  });
  document.getElementById("feature-type").addEventListener("change", async (e) => {
    state.currentFeatureType = e.target.value;
    const data = await loadAtlasModule(state.currentModule);
    populateFeatureSelect(data);
    if (state.currentFeature) await renderAtlasModule(data);
  });
  document.getElementById("feature-select").addEventListener("change", async (e) => {
    state.currentFeature = e.target.value || null;
    if (state.currentFeature) state.currentColorBy = "selected_feature";
    const data = await loadAtlasModule(state.currentModule);
    document.getElementById("color-by").value = state.currentColorBy;
    await renderAtlasModule(data);
  });
  document.getElementById("feature-search").addEventListener("input", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    populateFeatureSelect(data);
  });
  document.getElementById("violin-group-by").addEventListener("change", async () => {
    if (!state.currentFeature) return;
    const data = await loadAtlasModule(state.currentModule);
    await renderAtlasModule(data);
  });
  document.getElementById("point-size").addEventListener("input", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    await renderAtlasModule(data);
  });
  document.getElementById("point-opacity").addEventListener("input", async () => {
    if (["scarlink", "reference_mapping"].includes(state.currentModule)) return;
    const data = await loadAtlasModule(state.currentModule);
    await renderAtlasModule(data);
  });
  document.getElementById("reset-view").addEventListener("click", async () => {
    const data = await loadAtlasModule(state.currentModule);
    await renderAtlasModule(data);
  });
  document.getElementById("scarlink-load").addEventListener("click", async () => {
    await renderScarlink();
  });
  document.getElementById("scarlink-gene").addEventListener("change", async () => renderScarlink());
  document.getElementById("scarlink-disease").addEventListener("change", async () => {
    const gene = document.getElementById("scarlink-gene").value;
    const payload = await loadScarlinkPayload(gene, document.getElementById("scarlink-disease").value);
    setText("scarlink-caption", `${gene} in ${payload.disease}`);
    drawScarlinkCircle(payload);
    drawScarlinkBoxplot(payload);
    renderScarlinkTable(payload);
  });
  document.getElementById("scarlink-celltype-filter").addEventListener("input", async () => {
    if (!state.currentScarlinkGene || !state.currentScarlinkDisease) return;
    const payload = await loadScarlinkPayload(state.currentScarlinkGene, state.currentScarlinkDisease);
    renderScarlinkTable(payload);
  });
  await switchModule("whole_brain");
}

init().catch((err) => {
  setText("view-title", "Load error");
  setText("view-subtitle", err.message);
  console.error(err);
});
