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
  const payload = {entry, categories, summary, markers, cells, columns: ["x", "y", ...(entry.display_fields || ["cell_id", "subtype", "disease", "sample", "RL6", "RL_4", "RL3_2", "RL3_1", "RL_3", "RL_2"]).slice(2)]};
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
  const fieldIndex = Object.fromEntries(data.columns.map((field, idx) => [field, idx]));
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
    const vals = cells.map((row) => row[fieldIndex[state.currentColorBy] ?? fieldIndex.subtype]);
    colorLegend = categoricalColors(vals);
    markerColor = colorLegend.colors;
  }

  const hoverText = cells.map((row, idx) => {
    const parts = [];
    if (fieldIndex.cell_id !== undefined) parts.push(`Cell: ${escHtml(row[fieldIndex.cell_id])}`);
    parts.push(`Subtype: ${escHtml(row[fieldIndex.subtype])}`);
    parts.push(`Disease: ${escHtml(row[fieldIndex.disease])}`);
    parts.push(`Sample: ${escHtml(row[fieldIndex.sample])}`);
    ["RL6", "RL_4", "RL3_2", "RL3_1", "RL_3", "RL_2"].forEach((field) => {
      const v = row[fieldIndex[field]];
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

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function formatCoord(value) {
  return Number(value || 0).toLocaleString();
}

function mixHex(a, b, t) {
  const pa = a.match(/\w\w/g).map((x) => parseInt(x, 16));
  const pb = b.match(/\w\w/g).map((x) => parseInt(x, 16));
  const out = pa.map((v, i) => Math.round(v + (pb[i] - v) * t).toString(16).padStart(2, "0"));
  return `#${out.join("")}`;
}

function scarlinkLinkColor(link) {
  const sig = clamp(Number(link.significance || 0) / 6, 0, 1);
  const base = link.effect === "repression" ? "0f766e" : "b42318";
  const light = link.effect === "repression" ? "d8f3f0" : "fbe3de";
  return mixHex(light, base, 0.28 + sig * 0.72);
}

function tickLabelAttrs(angleDeg, x, y) {
  const flip = angleDeg > 90 || angleDeg < -90;
  return `x="${x.toFixed(2)}" y="${y.toFixed(2)}" font-size="11" fill="#475467" text-anchor="middle" transform="rotate(${flip ? angleDeg + 180 : angleDeg}, ${x.toFixed(2)}, ${y.toFixed(2)})"`;
}

function drawScarlinkCircle(data) {
  const el = document.getElementById("scarlink-circle");
  const links = (data.links || []).slice(0, 120);
  const topLinks = (data.top_links || []).slice(0, 5);
  if (!links.length) {
    el.innerHTML = `<div class="section-title"><p>No enhancer-gene links in the selected disease layer.</p></div>`;
    return;
  }
  const W = 920, H = 580, cx = 286, cy = 292, outer = 204, trackOuter = 176, trackInner = 160, inner = 122;
  const coords = [];
  links.forEach((l) => coords.push(l.enhancer_start, l.enhancer_end, l.tss, l.promoter_start, l.promoter_end));
  let minPos = Math.min(...coords), maxPos = Math.max(...coords);
  const pad = Math.max(25000, (maxPos - minPos) * 0.1);
  minPos -= pad;
  maxPos += pad;
  const span = Math.max(1, maxPos - minPos);
  const angle = (pos) => -225 + ((pos - minPos) / span) * 360;
  const majorTicks = 8;
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="100%"><rect width="${W}" height="${H}" rx="20" fill="white"/><circle cx="${cx}" cy="${cy}" r="${outer + 18}" fill="rgba(248,250,252,.92)" stroke="rgba(15,23,42,.05)"/>`;
  svg += `<text x="${cx}" y="38" text-anchor="middle" font-size="18" font-weight="700" fill="#1f2937">${escHtml(data.query.chr)} SCARlink enhancer-promoter map</text>`;
  svg += `<text x="${cx}" y="60" text-anchor="middle" font-size="12" fill="#667085">${escHtml(data.summary?.display_window || `${data.query.chr}:${formatCoord(minPos)}-${formatCoord(maxPos)}`)}</text>`;
  for (let i = 0; i < 6; i += 1) {
    const a1 = -220 + i * 60;
    const a2 = a1 + 36;
    svg += `<path d="${svgArc(cx, cy, trackOuter, a1, a2)}" fill="none" stroke="#2b7bbb" stroke-width="16" stroke-linecap="butt"/>`;
  }
  const promA1 = angle(data.query.promoter_start || data.query.start);
  const promA2 = angle(data.query.promoter_end || data.query.end);
  svg += `<path d="${svgArc(cx, cy, inner + 18, promA1, promA2)}" fill="none" stroke="#3da63a" stroke-width="8" stroke-linecap="round"/>`;
  topLinks.forEach((row) => {
    const a1 = angle(row.start || row.enhancer_start);
    const a2 = angle(row.end || row.enhancer_end);
    svg += `<path d="${svgArc(cx, cy, outer - 12, a1, a2)}" fill="none" stroke="#9aa0a6" stroke-width="7" stroke-linecap="round"/>`;
  });
  for (let i = 0; i <= majorTicks * 2; i += 1) {
    const pos = minPos + (span * i / (majorTicks * 2));
    const a = angle(pos);
    const p1 = svgPoint(cx, cy, outer - 4, a);
    const p2 = svgPoint(cx, cy, outer + (i % 2 === 0 ? 13 : 7), a);
    svg += `<line x1="${p1.x.toFixed(2)}" y1="${p1.y.toFixed(2)}" x2="${p2.x.toFixed(2)}" y2="${p2.y.toFixed(2)}" stroke="#111827" stroke-width="${i % 2 === 0 ? 1.3 : 0.8}"/>`;
    if (i % 2 === 0) {
      const labelPos = svgPoint(cx, cy, outer + 31, a);
      svg += `<text ${tickLabelAttrs(a, labelPos.x, labelPos.y)}>${(pos / 1e6).toFixed(2)}Mb</text>`;
    }
  }
  links.forEach((l) => {
    const ea = angle((l.enhancer_start + l.enhancer_end) / 2);
    const pa = angle(l.tss);
    const p1 = svgPoint(cx, cy, inner, ea);
    const p2 = svgPoint(cx, cy, inner, pa);
    const c1 = svgPoint(cx, cy, 28, ea);
    const c2 = svgPoint(cx, cy, 28, pa);
    const sig = clamp(Number(l.significance || 0) / 6, 0.12, 1);
    const width = l.is_top5 ? 4.8 : 0.8 + sig * 1.2;
    const opacity = l.is_top5 ? 0.96 : 0.14 + sig * 0.3;
    svg += `<path d="M ${p1.x.toFixed(2)} ${p1.y.toFixed(2)} C ${c1.x.toFixed(2)} ${c1.y.toFixed(2)}, ${c2.x.toFixed(2)} ${c2.y.toFixed(2)}, ${p2.x.toFixed(2)} ${p2.y.toFixed(2)}" fill="none" stroke="${scarlinkLinkColor(l)}" stroke-width="${width.toFixed(2)}" opacity="${opacity.toFixed(2)}"/>`;
  });
  svg += `<text x="${cx}" y="${cy - 12}" text-anchor="middle" font-size="26" font-weight="800" fill="#111827">${escHtml(data.query.chr)}</text>`;
  svg += `<text x="${cx}" y="${cy + 12}" text-anchor="middle" font-size="13" fill="#667085">${escHtml(data.gene)} promoter / TSS</text>`;
  svg += `<text x="${cx}" y="${cy + 32}" text-anchor="middle" font-size="12" fill="#475467">TSS ${formatCoord(data.query.tss)}</text>`;
  topLinks.forEach((row, idx) => {
    const mid = (row.start + row.end) / 2;
    const a = angle(mid);
    const p = svgPoint(cx, cy, outer + 55 + (idx % 2) * 12, a);
    const label = `${row.gene} ${formatCoord(row.start)}-${formatCoord(row.end)}`;
    svg += `<circle cx="${p.x.toFixed(2)}" cy="${p.y.toFixed(2)}" r="10" fill="${scarlinkLinkColor(row)}" opacity="0.92"/>`;
    svg += `<text x="${p.x.toFixed(2)}" y="${(p.y + 3.5).toFixed(2)}" text-anchor="middle" font-size="10" font-weight="700" fill="white">${idx + 1}</text>`;
    svg += `<text x="${(560).toFixed(2)}" y="${(314 + idx * 22).toFixed(2)}" font-size="12" fill="#344054">${idx + 1}. ${escHtml(label)} | FDR ${row.fdr}</text>`;
  });
  svg += `<text x="560" y="72" font-size="15" font-weight="700" fill="#1f2937">Track legend</text>`;
  [
    ["#2b7bbb", "genome track"],
    ["#3da63a", "promoter / TSS"],
    ["#9aa0a6", "enhancer tile (top links)"],
    [scarlinkLinkColor({effect: 'activation', significance: 5}), "activation link"],
    [scarlinkLinkColor({effect: 'repression', significance: 5}), "repression link"],
  ].forEach((row, idx) => {
    const y = 94 + idx * 24;
    svg += `<rect x="560" y="${y - 10}" width="16" height="10" rx="2" fill="${row[0]}"/><text x="584" y="${y - 1}" font-size="12" fill="#344054">${row[1]}</text>`;
  });
  svg += `<text x="560" y="252" font-size="15" font-weight="700" fill="#1f2937">Top 5 contacts</text>`;
  svg += `<text x="560" y="274" font-size="12" fill="#667085">Link color intensity scales with -log10(FDR).</text>`;
  svg += `<text x="560" y="300" font-size="12" fill="#667085">Showing ${Math.min(links.length, 120)} links; top 5 are labeled around the ring.</text>`;
  svg += `</svg>`;
  el.innerHTML = svg;
}

function drawScarlinkBoxplot(data) {
  const traces = Object.entries(data.box || {}).slice(0, 28).map(([name, values], idx) => ({
    type: "box",
    y: values,
    name,
    boxpoints: "all",
    jitter: 0.32,
    pointpos: -1.2,
    marker: {size: 4, opacity: 0.42, color: palette[idx % palette.length]},
    line: {color: palette[idx % palette.length]},
    fillcolor: palette[idx % palette.length],
    opacity: 0.55,
  }));
  Plotly.react("scarlink-boxplot", traces, {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "white",
    margin: {l: 55, r: 20, t: 22, b: 150},
    yaxis: {title: "z-score", gridcolor: "rgba(217,224,232,.65)"},
    xaxis: {tickangle: -35, automargin: true},
    showlegend: false,
  }, {responsive: true, displayModeBar: false});
}

function renderScarlinkTable(data) {
  const filter = document.getElementById("scarlink-celltype-filter").value.trim().toLowerCase();
  const rows = (data.table || []).filter((row) => !filter || String(row.celltype_r2).toLowerCase().includes(filter));
  document.getElementById("scarlink-head").innerHTML = "<tr><th>Rank</th><th>Disease</th><th>Cell type</th><th>Peak</th><th>Effect</th><th>Coef</th><th>FDR</th><th>-log10(FDR)</th><th>z-score</th></tr>";
  document.getElementById("scarlink-table").innerHTML = rows.slice(0, 220).map((row) => `<tr><td>${row.rank}</td><td>${escHtml(row.disease)}</td><td>${escHtml(row.celltype_r2)}</td><td>${escHtml(row.peak)}</td><td>${escHtml(row.effect)}</td><td>${row.regression_coef}</td><td>${row.fdr}</td><td>${row.significance}</td><td>${row.z_score}</td></tr>`).join("");
}

async function renderScarlink() {
  const manifest = await loadScarlinkManifest();
  document.getElementById("atlas-controls").classList.add("hidden");
  document.getElementById("scarlink-controls").classList.remove("hidden");
  document.getElementById("atlas-layout").classList.add("hidden");
  document.getElementById("scarlink-layout").classList.remove("hidden");
  document.getElementById("reference-layout").classList.add("hidden");
  document.getElementById("markers-panel").innerHTML = `<div class="marker-group">SCARlink examples are organized by gene and disease. Link color intensity reflects FDR significance, and the top 5 contacts are labeled in the circle map.</div>`;
  renderSummaryCards({n_exported_cells: manifest.genes.length, n_total_source_cells: manifest.diseases.length, n_subtypes: "SCARlink", n_diseases: manifest.diseases.length}, "Example", "Disease layers");
  setText("view-title", "SCARlink links");
  setText("view-subtitle", "Disease-aware static SCARlink view with coordinate-rich enhancer-promoter maps.");
  document.getElementById("module-note").textContent = "Circle plot, boxplot, and table are drawn from static JSON only. Top 5 contacts are highlighted and link color intensity follows FDR significance.";
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
  setText("scarlink-caption", `${gene} in ${payload.disease} | ${payload.summary?.query_region || "query"}`);
  drawScarlinkCircle(payload);
  drawScarlinkBoxplot(payload);
  renderScarlinkTable(payload);
}

async function loadReference() {
  const [summary, example, heatmap] = await Promise.all([
    getJSON("data/reference_mapping/summary.json"),
    getJSON("data/reference_mapping/example_mapping.json"),
    getJSON("data/reference_mapping/heatmap.json"),
  ]);
  return {summary, example, heatmap};
}

async function renderReference() {
  document.getElementById("atlas-controls").classList.add("hidden");
  document.getElementById("scarlink-controls").classList.add("hidden");
  document.getElementById("atlas-layout").classList.add("hidden");
  document.getElementById("scarlink-layout").classList.add("hidden");
  document.getElementById("reference-layout").classList.remove("hidden");
  const {summary, example, heatmap} = await loadReference();
  setText("view-title", "Reference mapping");
  setText("view-subtitle", summary.description);
  document.getElementById("module-note").textContent = `C5832Cd concordance ${(Number(summary.example_concordance || 0) * 100).toFixed(1)}% | ${summary.reference_label || "atlas second_label"}`;
  document.getElementById("markers-panel").innerHTML = `<div class="marker-group">Reference mapping uses the C5832Cd Huntington disease query labels, standardizes Inhibitory into Neuron, and compares the result with atlas second_label families. The layout below keeps the original query / transfer / concordance story compact for GitHub Pages.</div>`;
  renderSummaryCards({n_exported_cells: summary.modules.length, n_total_source_cells: summary.modules.length, n_subtypes: "Workflow", n_diseases: "Static"}, "Summary", "Table");
  const queryDataset = escHtml(summary.query_dataset || "C5832Cd");
  const standardization = (summary.query_standardization || []).map((rule) => `<span class="ref-chip">${escHtml(rule)}</span>`).join("");
  document.getElementById("reference-layout").innerHTML = `
    <div class="reference-grid">
      <div class="plot-card reference-panel reference-flow-panel">
        <div class="plot-card-head"><h3>Workflow</h3><span>${escHtml(summary.reference_label || "atlas second_label")} transfer overview</span></div>
        <div class="reference-flow">
          <div class="reference-node reference-node-source">
            <strong>${queryDataset}</strong>
            <span>Query label space</span>
            <small>Standardize raw Cluster labels</small>
          </div>
          <div class="reference-arrow">→</div>
          <div class="reference-node reference-node-target">
            <strong>Atlas</strong>
            <span>${escHtml(summary.reference_label || "second_label")}</span>
            <small>Reference family set</small>
          </div>
        </div>
        <div class="reference-note">
          <div class="ref-chip-row">${standardization}</div>
          <p>${escHtml(summary.workflow?.[0] || "Load the query label set and compare with atlas families.")}</p>
          <p>${escHtml(summary.workflow?.[1] || "Compare transferred labels against the reference atlas.")}</p>
          <p>${escHtml(summary.workflow?.[4] || "Potential next step: MIDAS-style joint training for missing-modality completion.")}</p>
          <ol>
            ${(summary.workflow || []).map((step) => `<li>${escHtml(step)}</li>`).join("")}
          </ol>
        </div>
      </div>
      <div class="plot-card reference-panel">
        <div class="plot-card-head"><h3>Example mapping table</h3><span>Exported modules</span></div>
        <div class="table-wrap"><table><thead><tr>${example.columns.map((c) => `<th>${escHtml(c)}</th>`).join("")}</tr></thead><tbody>${example.rows.map((row) => `<tr>${row.map((x) => `<td>${escHtml(x)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>
      </div>
    </div>`;
  if ((heatmap.query_labels || []).length && (heatmap.atlas_second_labels || []).length) {
    document.getElementById("reference-layout").innerHTML += `
      <div class="plot-card reference-panel reference-heatmap-panel">
        <div class="plot-card-head"><h3>Example concordance heatmap</h3><span>Query l1 labels vs atlas second_label families</span></div>
        <div id="reference-heatmap" class="plot reference-heatmap"></div>
      </div>
    `;
    Plotly.react("reference-heatmap", [{
      type: "heatmap",
      x: heatmap.atlas_second_labels,
      y: heatmap.query_labels,
      z: heatmap.matrix,
      text: heatmap.text,
      hovertemplate: "Query: %{y}<br>Atlas: %{x}<br>Score: %{z:.2f}<extra></extra>",
      colorscale: [
        [0, "#f7fbff"],
        [0.2, "#d7e8f8"],
        [0.45, "#9ecae1"],
        [0.7, "#4f8fc8"],
        [1, "#173f5f"],
      ],
      zmin: 0,
      zmax: 1,
      colorbar: {title: "Concordance"},
    }], {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "white",
      margin: {l: 170, r: 24, t: 20, b: 120},
      xaxis: {tickangle: -35, automargin: true},
      yaxis: {automargin: true},
      annotations: (heatmap.rows || []).map((row, idx) => ({
        xref: "paper",
        yref: "paper",
        x: 1.02,
        y: 1 - idx * 0.14,
        text: `${escHtml(row.query_label)}: ${(Number(row.best_match_score || 0) * 100).toFixed(0)}%`,
        showarrow: false,
        align: "left",
        font: {size: 12, color: "#344054"},
      })),
    }, {responsive: true, displayModeBar: false});
  }
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
    state.currentScarlinkGene = gene;
    state.currentScarlinkDisease = document.getElementById("scarlink-disease").value;
    setText("scarlink-caption", `${gene} in ${payload.disease} | ${payload.summary?.query_region || "query"}`);
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
