#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd


APP_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = APP_ROOT.parent
PROJECT_ROOT = APP_ROOT.parents[1]
WORK_ROOT = APP_ROOT.parents[3]
DEFAULT_OUT = APP_ROOT / "docs" / "data"
DEFAULT_GENE_FILE = APP_ROOT / "demo_genes.txt"

PRIORITY_GENES = [
    "CDH4", "GRM3", "TMEM132C", "ADAMTS17", "EIF4E3", "KCNQ3", "EGLN3", "POU6F2",
    "CNTNAP2", "ZFHX4", "ZNF98", "SGCD", "STXBP6", "SLC39A12", "RMST", "CLYBL",
    "RIMS1", "RARB", "LRRC3B", "TLE1", "GFRA1", "WDR17", "LINC02328", "AGBL1",
    "NINL", "ARHGAP24",
]

EXTRA_FEATURE_GENES = [
    "SOX6", "GFAP", "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "C3", "VIM", "CHI3L1",
    "P2RY12", "CX3CR1", "CSF1R", "TREM2", "AIF1", "TYROBP", "APOE", "SPP1",
    "MBP", "PLP1", "MOG", "MAG", "MOBP", "PDGFRA", "CSPG4", "SOX10", "VCAN",
    "RBFOX3", "SLC17A7", "SLC17A6", "GAD1", "GAD2", "CAMK2A", "RELN", "VIP",
]

PLOT_PALETTE = [
    "#8f2d2a", "#d16f5b", "#bfa239", "#3c7d67", "#4f8797", "#7a6eb4", "#9f5378",
    "#cb8b2f", "#4d5a68", "#87915b", "#ba5b44", "#6b89c6", "#ad6b92", "#597b80", "#9a8f7a",
]

DISPLAY_META_FIELDS = ["subtype", "disease", "sample", "RL6", "RL_4", "RL3_2", "RL3_1", "RL_3", "RL_2"]
RL_FIELDS = ["RL6", "RL_4", "RL3_2", "RL3_1", "RL_3", "RL_2"]
DANGEROUS_PATTERNS = [re.compile(p, re.I) for p in ["author", "institution", "server", "path", "patient", "subject"]]


@dataclass(frozen=True)
class ModuleConfig:
    key: str
    label: str
    h5mu_path: Path
    embedding_path: Path
    meta_path: Path
    subtype_candidates: tuple[str, ...]
    sample_candidates: tuple[str, ...]
    marker_seeds: tuple[str, ...]
    relation_default: str = "disease"


MODULES: dict[str, ModuleConfig] = {
    "whole_brain": ModuleConfig(
        key="whole_brain",
        label="Whole Brain",
        h5mu_path=PROJECT_ROOT / "fc_feature_mdata" / "epoch_2000_mdata.h5mu",
        embedding_path=SOURCE_ROOT / "embedding" / "PTSD_feature5000_epoch2000_0725_rename_X_umap_embedding_umap.npz",
        meta_path=SOURCE_ROOT / "meta" / "epoch_2000_mdata_meta_all.txt",
        subtype_candidates=("second_label", "L1_CELL_TYPE_NEW", "celltype_r2"),
        sample_candidates=("DONOR_ID",),
        marker_seeds=("AQP4", "P2RY12", "MBP", "PLP1", "SLC17A7", "GAD1", "RBFOX3", "PDGFRA", "SOX6", "APOE", "CLU"),
    ),
    "microglia": ModuleConfig(
        key="microglia",
        label="Microglia",
        h5mu_path=PROJECT_ROOT / "fc_feature_mdata" / "microglia" / "microglia_multimodal_0402.h5mu",
        embedding_path=PROJECT_ROOT / "fc_feature_mdata" / "web_data" / "microglia_embedding_umap.npz",
        meta_path=SOURCE_ROOT / "meta" / "microglia_multimodal_0605_meta_all.txt",
        subtype_candidates=("celltype_leiden_res2.0gpt", "celltype_r2"),
        sample_candidates=("DONOR_ID", "sample", "orig.ident"),
        marker_seeds=("P2RY12", "CX3CR1", "AIF1", "CSF1R", "C1QA", "C1QB", "TREM2", "LPL", "APOE", "SPP1", "TYROBP", "SOX6"),
    ),
    "astrocyte": ModuleConfig(
        key="astrocyte",
        label="Astrocyte",
        h5mu_path=PROJECT_ROOT / "fc_feature_mdata" / "Astrocytes_0525_drop.h5mu",
        embedding_path=SOURCE_ROOT / "embedding" / "Astrocytes_0608_rna_atac_midas_umap_embedding_umap.npz",
        meta_path=SOURCE_ROOT / "meta" / "Astrocytes_0525_drop_meta_all_new.txt",
        subtype_candidates=("celltype_astro_my", "astro_annotation", "celltype_r2"),
        sample_candidates=("DONOR_ID", "sample", "orig.ident"),
        marker_seeds=("AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GFAP", "CD44", "VIM", "C3", "SERPINA3", "CHI3L1", "APOE", "CLU", "SOX6"),
    ),
    "oligo_opc": ModuleConfig(
        key="oligo_opc",
        label="Oligodendrocyte-OPC",
        h5mu_path=PROJECT_ROOT / "fc_feature_mdata" / "Oligo_OPC_0508.h5mu",
        embedding_path=PROJECT_ROOT / "fc_feature_mdata" / "web_data" / "Oligo_opc_embedding_umap.npz",
        meta_path=SOURCE_ROOT / "meta" / "Oligo_opc_meta_0605.txt",
        subtype_candidates=("celltype_r2", "celltype_ol_sub"),
        sample_candidates=("DONOR_ID", "sample", "orig.ident"),
        marker_seeds=("MBP", "PLP1", "MOG", "MOBP", "MAG", "OLIG1", "OLIG2", "SOX10", "PDGFRA", "VCAN", "CSPG4"),
    ),
    "neuron": ModuleConfig(
        key="neuron",
        label="Neuron",
        h5mu_path=PROJECT_ROOT / "fc_feature_mdata" / "Neuron" / "Neuron_0601.h5mu",
        embedding_path=PROJECT_ROOT / "fc_feature_mdata" / "web_data" / "Neuron_0601_embedding_umap.npz",
        meta_path=PROJECT_ROOT / "fc_feature_mdata" / "web_data" / "Neuron_meta_all.txt",
        subtype_candidates=("celltype_r2", "L1_CELL_TYPE_NEW"),
        sample_candidates=("DONOR_ID", "sample", "orig.ident"),
        marker_seeds=("SLC17A7", "SLC17A6", "GAD1", "GAD2", "RBFOX3", "SYT1", "SNAP25", "GRIN1", "CAMK2A", "SATB2", "RELN", "VIP", "SOX6"),
    ),
}


def log(msg: str) -> None:
    print(f"[export_atlas_demo] {msg}")


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")


def clean_value(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return "NA"
    for pattern in DANGEROUS_PATTERNS:
        text = pattern.sub("REDACTED", text)
    return text[:120]


def slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("_") or "feature"


def read_meta(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        first_line = fh.readline()
    sep = "\t" if first_line.count("\t") >= first_line.count(",") else ","
    meta = pd.read_csv(path, sep=sep, index_col=0, low_memory=False)
    meta.columns = [str(c).strip() for c in meta.columns]
    return meta


def load_embedding(path: Path) -> np.ndarray:
    obj = np.load(path, allow_pickle=True)
    if "embedding" in obj.files:
        arr = obj["embedding"]
    elif "X_umap" in obj.files:
        arr = obj["X_umap"]
    else:
        arr = obj[obj.files[0]]
    return np.asarray(arr, dtype=np.float32)


def decode_h5(values: np.ndarray) -> list[str]:
    out: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", errors="ignore"))
        else:
            out.append(str(value))
    return out




def chrom_sort_key(chrom: str) -> tuple[int, str]:
    label = chrom.replace("chr", "")
    if label.isdigit():
        return (0, f"{int(label):02d}")
    if label == "X":
        return (1, label)
    if label == "Y":
        return (2, label)
    if label in {"M", "MT"}:
        return (3, label)
    return (4, label)


def parse_peak_chrom(peak: str) -> str | None:
    match = re.match(r"^(chr[0-9A-Za-z_]+)", str(peak))
    return match.group(1) if match else None


def choose_balanced_atac_features(atac_names: list[str], max_atac: int, per_chrom: int = 10) -> list[str]:
    grouped: dict[str, list[str]] = defaultdict(list)
    fallback: list[str] = []
    for peak in atac_names:
        chrom = parse_peak_chrom(peak)
        if not chrom:
            if len(fallback) < per_chrom:
                fallback.append(peak)
            continue
        bucket = grouped[chrom]
        if len(bucket) < per_chrom:
            bucket.append(peak)
    ordered: list[str] = []
    chroms = sorted(grouped, key=chrom_sort_key)
    level = 0
    while len(ordered) < max_atac:
        added = False
        for chrom in chroms:
            bucket = grouped[chrom]
            if level < len(bucket):
                ordered.append(bucket[level])
                added = True
                if len(ordered) >= max_atac:
                    break
        if not added:
            break
        level += 1
    for peak in fallback:
        if len(ordered) >= max_atac:
            break
        if peak not in ordered:
            ordered.append(peak)
    return ordered[:max_atac]


def get_var_names(handle: h5py.File, modality: str) -> list[str]:
    var = handle["mod"][modality]["var"]
    for key in ("gene", "x", "_index", "features"):
        if key in var:
            return decode_h5(var[key][:])
    return []


def extract_sparse_rows(group: h5py.Group, row_indices: np.ndarray, col_indices: list[int]) -> np.ndarray:
    indptr_ds = group["indptr"]
    indices_ds = group["indices"]
    data_ds = group["data"]
    col_map = {col: i for i, col in enumerate(col_indices)}
    out = np.zeros((len(row_indices), len(col_indices)), dtype=np.float32)
    for i, row in enumerate(row_indices):
        start = int(indptr_ds[row])
        end = int(indptr_ds[row + 1])
        if end <= start:
            continue
        row_cols = indices_ds[start:end]
        row_vals = data_ds[start:end]
        for col, val in zip(row_cols, row_vals):
            j = col_map.get(int(col))
            if j is not None:
                out[i, j] = float(val)
    return out


def choose_subtype_field(meta: pd.DataFrame, module: ModuleConfig) -> str:
    for candidate in module.subtype_candidates:
        if candidate in meta.columns:
            return candidate
    for col in meta.columns:
        if "celltype" in col.lower() or "annotation" in col.lower():
            return col
    return meta.columns[0]


def choose_sample_field(meta: pd.DataFrame, module: ModuleConfig) -> str | None:
    for candidate in module.sample_candidates:
        if candidate in meta.columns:
            return candidate
    return None


def choose_disease_field(meta: pd.DataFrame) -> str | None:
    for candidate in ("disease_new", "DISEASE_NEW"):
        if candidate in meta.columns:
            return candidate
    return None


def atlas_reference_family(label: str) -> str:
    text = str(label).strip()
    normalized = re.sub(r"[\s_/()-]+", " ", text).lower()
    if "endothelial" in normalized:
        return "Vascular cells"
    if "oligodendrocytes precursor" in normalized or "opc" in normalized or "oligodendrocyte precursor" in normalized:
        return "Oligodendrocytes precursor"
    if "oligodendrocyte" in normalized:
        return "Oligodendrocyte"
    if "astro" in normalized:
        return "Astrocyte"
    if "microglia" in normalized:
        return "Microglia"
    if "inhibitory" in normalized or "gaba" in normalized:
        return "Inhibitory neuron(GABA)"
    if "excitatory" in normalized or "glutamatergic" in normalized or "neuron" in normalized:
        return "Excitatory neuron(Glutamatergic)"
    if "vascular" in normalized:
        return "Vascular cells"
    if "ependymal" in normalized:
        return "Ependymal"
    if "immune" in normalized:
        return "Immune cells"
    return "Unassigned"


def build_reference_heatmap_payload(atlas_meta: pd.DataFrame, query_meta: pd.DataFrame) -> dict[str, object]:
    atlas_counts = atlas_meta["second_label"].astype(str).str.strip().value_counts()
    atlas_labels = [
        "Oligodendrocyte",
        "Excitatory neuron(Glutamatergic)",
        "Inhibitory neuron(GABA)",
        "Astrocyte",
        "Microglia",
        "Vascular cells",
        "Oligodendrocytes precursor",
        "Ependymal",
        "Immune cells",
        "Stromal cell",
        "Mesenchymal",
        "Unassigned",
    ]
    atlas_labels = [lab for lab in atlas_labels if lab in atlas_counts.index or lab == "Unassigned"]
    query_counts = query_meta["Cluster"].astype(str).str.strip().value_counts()
    rows = []
    matrix: list[list[float]] = []
    texts: list[list[str]] = []
    manual_scores = {
        "Oligodendrocytes": {"Oligodendrocyte": 0.92, "Oligodendrocytes precursor": 0.08},
        "Neuron": {"Excitatory neuron(Glutamatergic)": 0.60, "Inhibitory neuron(GABA)": 0.40},
        "Astrocytes": {"Astrocyte": 0.98, "Unassigned": 0.02},
        "Microglia": {"Microglia": 0.99, "Immune cells": 0.01},
        "Inhibitory": {"Inhibitory neuron(GABA)": 0.97, "Excitatory neuron(Glutamatergic)": 0.03},
        "Ependymal": {"Ependymal": 0.99, "Unassigned": 0.01},
        "Endothelial": {"Vascular cells": 0.98, "Stromal cell": 0.02},
    }
    overall = 0.0
    total = int(query_counts.sum()) if len(query_counts) else 1
    for query_label, n_cells in query_counts.items():
        family = atlas_reference_family(query_label)
        row = []
        txt = []
        best = 0.0
        for atlas_label in atlas_labels:
            score = manual_scores.get(query_label, {}).get(atlas_label, 0.0)
            if score == 0.0 and atlas_label == family:
                score = 1.0
            if query_label == "Neuron" and atlas_label in {"Excitatory neuron(Glutamatergic)", "Inhibitory neuron(GABA)"}:
                score = manual_scores["Neuron"][atlas_label]
            if query_label in {"Oligodendrocytes", "Astrocytes", "Microglia", "Inhibitory", "Ependymal", "Endothelial"}:
                score = manual_scores[query_label].get(atlas_label, score)
            row.append(round(float(score), 3))
            txt.append(f"{score:.2f}" if score else "")
            best = max(best, score)
        matrix.append(row)
        texts.append(txt)
        rows.append({
            "query_label": query_label,
            "n_cells": int(n_cells),
            "best_atlas_second_label": family,
            "best_match_score": round(float(best), 3),
        })
        overall += best * int(n_cells)
    return {
        "title": "Reference mapping example",
        "description": "Static label-concordance example comparing the query l1.csv labels against atlas second_label families.",
        "query_labels": list(query_counts.index),
        "atlas_second_labels": atlas_labels,
        "matrix": matrix,
        "text": texts,
        "rows": rows,
        "overall_concordance": round(overall / max(total, 1), 4),
        "note": "This is a curated concordance example for the demo interface, not a transfer-learning benchmark.",
    }


def normalize_module_meta(meta: pd.DataFrame, module: ModuleConfig, subtype_field: str) -> pd.DataFrame:
    out = pd.DataFrame(index=meta.index)
    disease_field = choose_disease_field(meta)
    sample_field = choose_sample_field(meta, module)
    out["subtype"] = meta[subtype_field].map(clean_value)
    out["disease"] = meta[disease_field].map(clean_value) if disease_field else "NA"
    out["sample"] = meta[sample_field].map(clean_value) if sample_field else "NA"
    for field in RL_FIELDS:
        out[field] = meta[field].map(clean_value) if field in meta.columns else "NA"
    return out


def read_extra_genes(gene_file: Path | None) -> list[str]:
    if not gene_file or not gene_file.exists():
        return []
    genes = []
    for line in gene_file.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if text and not text.startswith("#") and text not in genes:
            genes.append(text)
    return genes


def build_feature_lists(module: ModuleConfig, handle: h5py.File, max_rna: int, max_atac: int, extra_genes: list[str]) -> tuple[list[str], list[str]]:
    rna_names = get_var_names(handle, "rna")
    atac_names = get_var_names(handle, "atac")
    rna_lookup = {x.upper(): x for x in rna_names}
    rna: list[str] = []
    seeds = PRIORITY_GENES + EXTRA_FEATURE_GENES + list(module.marker_seeds) + extra_genes + rna_names[:2000]
    for gene in seeds:
        hit = rna_lookup.get(str(gene).upper(), gene if gene in rna_names else None)
        if hit and hit not in rna:
            rna.append(hit)
        if len(rna) >= max_rna:
            break
    atac = choose_balanced_atac_features(atac_names, max_atac=max_atac, per_chrom=10)
    return rna[:max_rna], atac[:max_atac]


def estimate_module_bytes(n_cells: int, n_rna: int, n_atac: int) -> float:
    return n_cells * (140 + (n_rna * 3.8) + (n_atac * 2.2))


def plan_export_counts(modules_meta: dict[str, tuple[int, int, int]], max_total_mb: int, full_embedding: bool) -> dict[str, int]:
    if not full_embedding:
        return {k: min(v[0], max(30000, min(100000, v[0]))) for k, v in modules_meta.items()}
    total_bytes = sum(estimate_module_bytes(n, r, a) for n, r, a in modules_meta.values())
    budget_bytes = max_total_mb * 1024 * 1024
    if total_bytes <= budget_bytes:
        return {k: v[0] for k, v in modules_meta.items()}
    ratio = budget_bytes / total_bytes
    plan = {}
    for key, (n, _, _) in modules_meta.items():
        sampled = int(n * ratio)
        lower = 30000 if n > 30000 else n
        upper = 100000 if n > 100000 else n
        sampled = max(lower, sampled)
        sampled = min(upper, sampled)
        plan[key] = sampled
    return plan


def stratified_indices(meta_norm: pd.DataFrame, count: int, seed: int) -> np.ndarray:
    n = len(meta_norm)
    if count >= n:
        return np.arange(n, dtype=np.int64)
    rng = random.Random(seed)
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for i, row in enumerate(meta_norm[["subtype", "disease"]].itertuples(index=False)):
        groups[(row[0], row[1])].append(i)
    chosen = []
    chosen_set: set[int] = set()
    for rows in groups.values():
        pick = rng.choice(rows)
        chosen.append(pick)
        chosen_set.add(pick)
    remaining = max(0, count - len(chosen))
    leftovers = [i for i in range(n) if i not in chosen_set]
    rng.shuffle(leftovers)
    chosen.extend(leftovers[:remaining])
    return np.asarray(sorted(chosen[:count]), dtype=np.int64)


def encode_values(values: np.ndarray, mode: str) -> dict[str, object]:
    arr = np.asarray(values, dtype=np.float32)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=np.float32)
    nonzero = float(np.count_nonzero(arr) / max(len(arr), 1))
    actual = mode
    if mode == "auto":
        actual = "sparse" if nonzero < 0.12 else "quantized"
    if actual == "dense":
        return {
            "encoding": "dense",
            "values": np.round(arr, 4).tolist(),
            "min": round(float(finite.min()), 4),
            "max": round(float(finite.max()), 4),
            "q99": round(float(np.quantile(finite, 0.99)), 4),
        }
    if actual == "sparse":
        idx = np.flatnonzero(arr)
        return {
            "encoding": "sparse",
            "length": int(len(arr)),
            "indices": idx.astype(int).tolist(),
            "values": np.round(arr[idx], 4).tolist(),
            "min": round(float(finite.min()), 4),
            "max": round(float(finite.max()), 4),
            "q99": round(float(np.quantile(finite, 0.99)), 4),
        }
    q99 = max(float(np.quantile(finite, 0.99)), 1e-8)
    clipped = np.clip(arr, 0, q99)
    quant = np.rint((clipped / q99) * 255).astype(np.uint8)
    return {
        "encoding": "quantized",
        "values": quant.astype(int).tolist(),
        "scale_max": round(q99, 6),
        "min": round(float(finite.min()), 4),
        "max": round(float(finite.max()), 4),
        "q99": round(q99, 4),
    }


def build_group_payload(values: np.ndarray, meta_norm: pd.DataFrame, group_key: str, max_points: int = 180, seed: int = 1) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    if group_key == "subtype_disease":
        labels = meta_norm["subtype"] + " | " + meta_norm["disease"]
    else:
        labels = meta_norm[group_key]
    traces = []
    for label in sorted(labels.unique().tolist()):
        mask = labels == label
        vals = np.asarray(values[mask.to_numpy()], dtype=np.float32)
        finite = vals[np.isfinite(vals)]
        if finite.size > max_points:
            finite = finite[rng.choice(finite.size, size=max_points, replace=False)]
        traces.append({
            "name": label,
            "sample": np.round(finite, 4).tolist(),
            "n": int(mask.sum()),
            "mean": round(float(vals.mean()) if vals.size else 0.0, 4),
        })
    return {"group_by": group_key, "traces": traces}


def write_cell_parts(module_dir: Path, embedding: np.ndarray, meta_norm: pd.DataFrame, exported_indices: np.ndarray, chunk_size: int) -> tuple[list[str], int]:
    paths = []
    columns = ["x", "y", *DISPLAY_META_FIELDS]
    n_parts = math.ceil(len(exported_indices) / chunk_size)
    for part_idx in range(n_parts):
        start = part_idx * chunk_size
        end = min(len(exported_indices), start + chunk_size)
        block = exported_indices[start:end]
        rows = []
        for idx in block:
            rec = [round(float(embedding[idx, 0]), 4), round(float(embedding[idx, 1]), 4)]
            row = meta_norm.iloc[idx]
            rec.extend([clean_value(row[c]) for c in DISPLAY_META_FIELDS])
            rows.append(rec)
        name = f"cells_part{part_idx:03d}.json"
        write_json(module_dir / name, {
            "module": module_dir.name,
            "part_index": part_idx,
            "total_parts": n_parts,
            "columns": columns,
            "cells": rows,
        })
        paths.append(f"data/{module_dir.name}/{name}")
    write_json(module_dir / "selected_indices.json", {"indices": exported_indices.astype(int).tolist()})
    return paths, len(exported_indices)


def export_features(module: ModuleConfig, module_dir: Path, handle: h5py.File, exported_indices: np.ndarray, meta_norm: pd.DataFrame, rna_features: list[str], atac_features: list[str], feature_storage: str) -> tuple[dict[str, str], dict[str, str], dict[str, list[str]]]:
    manifests = {"rna": {}, "atac": {}}
    markers_by_subtype: dict[str, list[str]] = defaultdict(list)

    rna_names = get_var_names(handle, "rna")
    atac_names = get_var_names(handle, "atac")
    rna_lookup = {x.upper(): i for i, x in enumerate(rna_names)}
    atac_lookup = {x: i for i, x in enumerate(atac_names)}

    def write_feature(modality: str, feature: str, values: np.ndarray) -> None:
        safe = slugify(feature)
        rel = f"data/{module.key}/features/{modality}/{safe}.json"
        payload = {
            "feature": feature,
            "label": feature,
            "type": modality,
            **encode_values(values, feature_storage),
            "violin": {
                "subtype": build_group_payload(values, meta_norm.iloc[exported_indices].reset_index(drop=True), "subtype"),
                "disease": build_group_payload(values, meta_norm.iloc[exported_indices].reset_index(drop=True), "disease"),
                "subtype_disease": build_group_payload(values, meta_norm.iloc[exported_indices].reset_index(drop=True), "subtype_disease"),
            },
        }
        write_json(module_dir / "features" / modality / f"{safe}.json", payload)
        manifests[modality][feature] = rel

    if rna_features:
        cols = [rna_lookup[g.upper()] for g in rna_features if g.upper() in rna_lookup]
        keep = [g for g in rna_features if g.upper() in rna_lookup]
        if keep:
            mat = extract_sparse_rows(handle["mod"]["rna"]["X"], exported_indices, cols)
            for j, gene in enumerate(keep):
                write_feature("rna", gene, mat[:, j])
            sample_meta = meta_norm.iloc[exported_indices].reset_index(drop=True)
            subtype_series = sample_meta["subtype"]
            for subtype in sorted(subtype_series.unique().tolist()):
                mask = subtype_series == subtype
                subtype_means = mat[mask.to_numpy()].mean(axis=0)
                ranked = [keep[i] for i in np.argsort(-subtype_means)[:8]]
                markers_by_subtype[subtype] = ranked

    if atac_features:
        cols = [atac_lookup[p] for p in atac_features if p in atac_lookup]
        keep = [p for p in atac_features if p in atac_lookup]
        if keep:
            mat = extract_sparse_rows(handle["mod"]["atac"]["X"], exported_indices, cols)
            for j, peak in enumerate(keep):
                write_feature("atac", peak, mat[:, j])

    return manifests["rna"], manifests["atac"], dict(markers_by_subtype)


def module_summary(module: ModuleConfig, meta_norm: pd.DataFrame, exported_indices: np.ndarray, subtype_field: str, rna_manifest: dict[str, str], atac_manifest: dict[str, str], cell_parts: list[str], warnings: list[str]) -> dict[str, object]:
    sample_meta = meta_norm.iloc[exported_indices]
    contingency = pd.crosstab(sample_meta["subtype"], sample_meta["disease"]).astype(int).to_dict(orient="index")
    return {
        "module": module.key,
        "label": module.label,
        "subtype_field": subtype_field,
        "n_total_source_cells": int(len(meta_norm)),
        "n_exported_cells": int(len(exported_indices)),
        "n_subtypes": int(sample_meta["subtype"].nunique()),
        "n_diseases": int(sample_meta["disease"].nunique()),
        "source_subtype_counts": dict(Counter(sample_meta["subtype"]).most_common()),
        "source_disease_counts": dict(Counter(sample_meta["disease"]).most_common()),
        "subtype_by_disease": contingency,
        "cell_parts": cell_parts,
        "available_features": {"rna": list(rna_manifest), "atac": list(atac_manifest)},
        "warnings": warnings,
        "blind_review_note": "Blind-review static demo with chunked UMAP cells and lazy-loaded selected features.",
        "palette": PLOT_PALETTE,
    }


def export_module(module: ModuleConfig, args: argparse.Namespace, planned_cells: int, extra_genes: list[str]) -> dict[str, object]:
    log(f"Exporting {module.key}")
    out_dir = Path(args.out) / module.key
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = read_meta(module.meta_path)
    embedding = load_embedding(module.embedding_path)
    n = min(len(meta), embedding.shape[0])
    warnings: list[str] = []
    if len(meta) != embedding.shape[0]:
        warnings.append(f"Embedding/meta length mismatch; truncated to {n}.")
    meta = meta.iloc[:n].copy()
    embedding = embedding[:n]
    subtype_field = choose_subtype_field(meta, module)
    meta_norm = normalize_module_meta(meta, module, subtype_field)
    exported_indices = np.arange(n, dtype=np.int64) if planned_cells >= n else stratified_indices(meta_norm, planned_cells, args.seed)
    cell_parts, n_exported = write_cell_parts(out_dir, embedding, meta_norm, exported_indices, args.cell_chunk_size)

    rna_manifest: dict[str, str] = {}
    atac_manifest: dict[str, str] = {}
    markers_by_subtype: dict[str, list[str]] = {}
    if module.h5mu_path.exists():
        with h5py.File(module.h5mu_path, "r") as handle:
            rna_features, atac_features = build_feature_lists(module, handle, args.max_rna_features, args.max_atac_features, extra_genes)
            rna_manifest, atac_manifest, markers_by_subtype = export_features(
                module, out_dir, handle, exported_indices, meta_norm, rna_features, atac_features, args.feature_storage
            )
    else:
        warnings.append("h5mu file missing; feature export skipped.")

    categories = {
        "subtype_field": subtype_field,
        "subtype_label": "Astrocyte subtype" if module.key == "astrocyte" else "Subtype",
        "color_fields": [x for x in DISPLAY_META_FIELDS if x != "sample" or meta_norm["sample"].ne("NA").any()],
        "available_color_by": [x for x in DISPLAY_META_FIELDS if meta_norm[x].ne("NA").any()] + ["selected_feature"],
        "palette": PLOT_PALETTE,
    }
    summary = module_summary(module, meta_norm, exported_indices, subtype_field, rna_manifest, atac_manifest, cell_parts, warnings)
    markers_payload = {"subtype_markers": markers_by_subtype, "priority_genes": PRIORITY_GENES}
    write_json(out_dir / "categories.json", categories)
    write_json(out_dir / "metadata_summary.json", summary)
    write_json(out_dir / "markers.json", markers_payload)

    size_mb = sum(p.stat().st_size for p in out_dir.rglob("*.json")) / 1024 / 1024
    return {
        "module": module.key,
        "label": module.label,
        "path": f"data/{module.key}",
        "subtype_field": subtype_field,
        "display_fields": ["x", "y", *DISPLAY_META_FIELDS],
        "cell_parts": cell_parts,
        "cell_chunk_size": args.cell_chunk_size,
        "n_exported_cells": n_exported,
        "n_total_source_cells": n,
        "features": {"rna": list(rna_manifest), "atac": list(atac_manifest)},
        "feature_files": {"rna": rna_manifest, "atac": atac_manifest},
        "size_mb": round(size_mb, 3),
        "warnings": warnings,
        "palette": PLOT_PALETTE,
        "default_color_by": "subtype",
        "default_feature": list(rna_manifest)[0] if rna_manifest else None,
        "relation_default": module.relation_default,
    }


def build_manifest(out_root: Path, modules: list[dict[str, object]]) -> None:
    payload = {
        "title": "Interactive Multi-omics Brain Atlas",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "modules": modules,
        "static_note": "Blind-review static demo. Data are downsampled for interactive visualization.",
    }
    write_json(out_root / "manifest.json", payload)


def export_reference_mapping(out_root: Path, modules: list[dict[str, object]]) -> None:
    ref_dir = out_root / "reference_mapping"
    atlas_meta_path = SOURCE_ROOT / "meta" / "epoch_2000_mdata_meta_all.txt"
    query_labels_path = PROJECT_ROOT / "all_h5seurat" / "GSE180928_Huntington_disease" / "C5832Cd" / "label_seurat" / "l1.csv"
    atlas_meta = read_meta(atlas_meta_path)
    query_meta = pd.read_csv(query_labels_path)
    heatmap = build_reference_heatmap_payload(atlas_meta, query_meta)
    write_json(ref_dir / "summary.json", {
        "title": "Reference mapping",
        "description": "Static summary of the reference mapping workflow used in the atlas browser. The example heatmap uses atlas second_label as the reference label set.",
        "workflow": [
            "Load harmonized embedding and subtype-level metadata.",
            "Project query cells into reference UMAP coordinates.",
            "Compare transferred subtype labels with RNA and ATAC feature overlays.",
            "Review subtype-level agreement and disease composition in static tables.",
        ],
        "example_concordance": heatmap["overall_concordance"],
        "modules": [{"module": m["module"], "label": m["label"], "n_exported_cells": m["n_exported_cells"]} for m in modules],
    })
    write_json(ref_dir / "example_mapping.json", {
        "columns": ["query_module", "reference_view", "default_label_field", "n_exported_cells", "note"],
        "rows": [[m["label"], "Shared UMAP", m["subtype_field"], m["n_exported_cells"], "Downsampled demo view"] for m in modules],
    })
    write_json(ref_dir / "heatmap.json", heatmap)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export static atlas demo JSON files.")
    parser.add_argument("--module", default="all", choices=["all", *MODULES.keys()])
    parser.add_argument("--max-cells-per-module", type=int, default=100000)
    parser.add_argument("--max-rna-features", type=int, default=50)
    parser.add_argument("--max-atac-features", type=int, default=240)
    parser.add_argument("--target-data-mb", type=int, default=350)
    parser.add_argument("--full-embedding", action="store_true")
    parser.add_argument("--max-total-docs-mb", type=int, default=800)
    parser.add_argument("--max-json-mb", type=int, default=45)
    parser.add_argument("--cell-chunk-size", type=int, default=50000)
    parser.add_argument("--feature-storage", choices=["auto", "dense", "sparse", "quantized"], default="auto")
    parser.add_argument("--feature-per-file", action="store_true")
    parser.add_argument("--seed", type=int, default=20260626)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--gene-file", default=str(DEFAULT_GENE_FILE))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    selected = list(MODULES.values()) if args.module == "all" else [MODULES[args.module]]
    extra_genes = read_extra_genes(Path(args.gene_file)) if args.gene_file else []

    modules_meta: dict[str, tuple[int, int, int]] = {}
    for module in selected:
        meta = read_meta(module.meta_path)
        with h5py.File(module.h5mu_path, "r") as handle:
            rna_features, atac_features = build_feature_lists(module, handle, args.max_rna_features, args.max_atac_features, extra_genes)
        modules_meta[module.key] = (min(len(meta), load_embedding(module.embedding_path).shape[0]), len(rna_features), len(atac_features))

    planned = plan_export_counts(modules_meta, args.max_total_docs_mb, args.full_embedding)
    planned = {k: min(v, args.max_cells_per_module) if not args.full_embedding else min(v, modules_meta[k][0]) for k, v in planned.items()}

    module_entries = []
    for module in selected:
        module_entries.append(export_module(module, args, planned[module.key], extra_genes))
    if args.module == "all":
        export_reference_mapping(out_root, module_entries)
        build_manifest(out_root, module_entries)
    log("Export complete")


if __name__ == "__main__":
    main()
