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

import h5py
import numpy as np
import pandas as pd


APP_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = APP_ROOT / "docs" / "data"
CACHE_ROOT = APP_ROOT / "atlas_light_cache"

PRIORITY_GENES = [
    "CDH4", "GRM3", "TMEM132C", "ADAMTS17", "EIF4E3", "KCNQ3", "EGLN3", "POU6F2",
    "CNTNAP2", "ZFHX4", "ZNF98", "SGCD", "STXBP6", "SLC39A12", "RMST", "CLYBL",
    "RIMS1", "RARB", "LRRC3B", "TLE1", "GFRA1", "WDR17", "LINC02328", "AGBL1",
    "NINL", "ARHGAP24",
]

DANGEROUS_COL_PATTERNS = [
    re.compile(p, re.I)
    for p in [
        "donor", "patient", "subject", "author", "institution", "server", "path",
        "unique_id", "barcode", "sample_id", "precise", "sub_region_precise",
    ]
]


@dataclass(frozen=True)
class ModuleConfig:
    key: str
    label: str
    h5mu_path: Path
    embedding_path: Path
    meta_path: Path
    subtype_field: str
    subtype_fallbacks: tuple[str, ...]
    marker_seeds: tuple[str, ...]


MODULES: dict[str, ModuleConfig] = {
    "whole_brain": ModuleConfig(
        key="whole_brain",
        label="Whole Brain",
        h5mu_path=Path("/share/home/fxl/project/Atlas_1230/fc_feature_mdata/epoch_2000_mdata.h5mu"),
        embedding_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/embedding/epoch_2000_mdata_embedding_umap.npz"),
        meta_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/meta/epoch_2000_mdata_meta_all.txt"),
        subtype_field="L1_CELL_TYPE_NEW",
        subtype_fallbacks=("celltype_r2",),
        marker_seeds=("AQP4", "P2RY12", "MBP", "PLP1", "SLC17A7", "GAD1", "RBFOX3", "PDGFRA", "SOX6", "APOE", "CLU"),
    ),
    "microglia": ModuleConfig(
        key="microglia",
        label="Microglia",
        h5mu_path=Path("/share/home/fxl/project/Atlas_1230/fc_feature_mdata/microglia/microglia_multimodal_0402.h5mu"),
        embedding_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/embedding/microglia_embedding_umap.npz"),
        meta_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/meta/microglia_multimodal_0605_meta_all.txt"),
        subtype_field="celltype_leiden_res2.0gpt",
        subtype_fallbacks=("celltype_r2",),
        marker_seeds=("P2RY12", "CX3CR1", "AIF1", "CSF1R", "C1QA", "C1QB", "TREM2", "LPL", "APOE", "SPP1"),
    ),
    "astrocyte": ModuleConfig(
        key="astrocyte",
        label="Astrocyte",
        h5mu_path=Path("/share/home/fxl/project/Atlas_1230/fc_feature_mdata/Astrocytes_0525_drop.h5mu"),
        embedding_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/embedding/Astrocytes_0525_drop_embedding_umap.npz"),
        meta_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/meta/Astrocytes_0525_drop_meta_all_0605.txt"),
        subtype_field="astro_annotation",
        subtype_fallbacks=("celltype_r2",),
        marker_seeds=("AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "GFAP", "CD44", "VIM", "C3", "SERPINA3", "CHI3L1"),
    ),
    "oligo_opc": ModuleConfig(
        key="oligo_opc",
        label="Oligodendrocyte-OPC",
        h5mu_path=Path("/share/home/fxl/project/Atlas_1230/fc_feature_mdata/Oligo_OPC_0508.h5mu"),
        embedding_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/embedding/Oligo_opc_embedding_umap.npz"),
        meta_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/meta/Oligo_opc_meta_0605.txt"),
        subtype_field="celltype_r2",
        subtype_fallbacks=("celltype_ol_sub",),
        marker_seeds=("MBP", "PLP1", "MOG", "MOBP", "MAG", "OLIG1", "OLIG2", "SOX10", "PDGFRA", "VCAN"),
    ),
    "neuron": ModuleConfig(
        key="neuron",
        label="Neuron",
        h5mu_path=Path("/share/home/fxl/project/Atlas_1230/fc_feature_mdata/Neuron/Neuron_0601.h5mu"),
        embedding_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/embedding/Neuron_0601_embedding_umap.npz"),
        meta_path=Path("/share/home/fxl/project/Atlas_1230/web_可视化/meta/Neuron_meta_all.txt"),
        subtype_field="celltype_r2",
        subtype_fallbacks=("L1_CELL_TYPE_NEW",),
        marker_seeds=("SLC17A7", "SLC17A6", "GAD1", "GAD2", "RBFOX3", "SYT1", "SNAP25", "GRIN1", "CAMK2A", "SATB2"),
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
    if not text:
        return "NA"
    text = re.sub(r"/share/home/\S+", "REDACTED", text)
    text = re.sub(r"/root/\S+", "REDACTED", text)
    text = re.sub(r"Atlas_1230", "REDACTED", text, flags=re.I)
    return text[:120]


def safe_column(name: str) -> bool:
    lname = name.lower()
    if lname.startswith("unnamed"):
        return False
    return not any(p.search(name) for p in DANGEROUS_COL_PATTERNS)


def read_meta(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
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


def choose_subtype_field(meta: pd.DataFrame, module: ModuleConfig) -> str:
    if module.subtype_field in meta.columns:
        return module.subtype_field
    for alt in module.subtype_fallbacks:
        if alt in meta.columns:
            return alt
    for col in meta.columns:
        if "celltype" in col.lower() or "annotation" in col.lower():
            return col
    return meta.columns[0]


def pick_display_columns(meta: pd.DataFrame, subtype_field: str) -> list[str]:
    candidates = [
        subtype_field, "disease_new", "DISEASE_NEW", "sample", "SAMPLE", "batch",
        "rna:batch", "dataset", "DATASET", "Publication", "REGION_L1", "REGION_L3",
    ]
    out: list[str] = []
    for col in candidates:
        if col in meta.columns and safe_column(col) and col not in out:
            out.append(col)
    if subtype_field not in out:
        out.insert(0, subtype_field)
    return out[:6]


def normalize_display_meta(meta: pd.DataFrame, subtype_field: str) -> pd.DataFrame:
    cols = pick_display_columns(meta, subtype_field)
    out = pd.DataFrame(index=meta.index)
    out["subtype"] = meta[subtype_field].map(clean_value)
    disease_col = "disease_new" if "disease_new" in meta.columns else "DISEASE_NEW" if "DISEASE_NEW" in meta.columns else None
    out["disease"] = meta[disease_col].map(clean_value) if disease_col else "NA"
    sample_col = next((c for c in ["sample", "SAMPLE", "batch", "rna:batch"] if c in meta.columns and safe_column(c)), None)
    dataset_col = next((c for c in ["dataset", "DATASET", "Publication", "DATA_MODALITY"] if c in meta.columns and safe_column(c)), None)
    if sample_col:
        out["sample"] = meta[sample_col].map(clean_value)
    else:
        out["sample"] = "NA"
    if dataset_col:
        out["dataset"] = meta[dataset_col].map(clean_value)
    else:
        out["dataset"] = "NA"
    return out


def capped_cells(max_cells: int, rna_features: int, atac_features: int, target_mb: int) -> int:
    bytes_per_cell = (2 + rna_features + atac_features) * 4 + 64
    rough_cap = max(1000, int((target_mb * 1024 * 1024 * 0.75) / max(bytes_per_cell, 1)))
    return min(max_cells, rough_cap)


def stratified_sample(meta: pd.DataFrame, subtype_field: str, max_cells: int, seed: int) -> np.ndarray:
    disease_col = "disease_new" if "disease_new" in meta.columns else "DISEASE_NEW" if "DISEASE_NEW" in meta.columns else None
    rng = random.Random(seed)
    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    subtype_vals = meta[subtype_field].map(clean_value)
    disease_vals = meta[disease_col].map(clean_value) if disease_col else pd.Series(["NA"] * len(meta), index=meta.index)
    for i, (subtype, disease) in enumerate(zip(subtype_vals, disease_vals)):
        groups[(subtype, disease)].append(i)
    if len(meta) <= max_cells:
        return np.arange(len(meta))
    picked: list[int] = []
    for rows in groups.values():
        picked.append(rng.choice(rows))
    remaining = max(0, max_cells - len(picked))
    leftovers = [idx for rows in groups.values() for idx in rows if idx not in set(picked)]
    rng.shuffle(leftovers)
    picked.extend(leftovers[:remaining])
    picked = sorted(set(picked))
    if len(picked) > max_cells:
        rng.shuffle(picked)
        picked = sorted(picked[:max_cells])
    return np.asarray(picked, dtype=np.int64)


def decode_h5(values: np.ndarray) -> list[str]:
    out = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", errors="ignore"))
        else:
            out.append(str(value))
    return out


def get_var_names(handle: h5py.File, modality: str) -> list[str]:
    var = handle["mod"][modality]["var"]
    for key in ("gene", "x", "_index", "features"):
        if key in var:
            return decode_h5(var[key][:])
    return []


def extract_sparse_rows(group: h5py.Group, row_indices: np.ndarray, col_indices: list[int]) -> np.ndarray:
    if not {"data", "indices", "indptr"} <= set(group.keys()):
        raise ValueError("Unsupported sparse matrix layout")
    indptr_ds = group["indptr"]
    indices_ds = group["indices"]
    data_ds = group["data"]
    col_map = {col: j for j, col in enumerate(col_indices)}
    out = np.zeros((len(row_indices), len(col_indices)), dtype=np.float32)
    for i, row in enumerate(row_indices):
        ptr = indptr_ds[row: row + 2]
        start = int(ptr[0])
        end = int(ptr[1])
        if end <= start:
            continue
        row_cols = indices_ds[start:end]
        row_vals = data_ds[start:end]
        for col, val in zip(row_cols, row_vals):
            j = col_map.get(int(col))
            if j is not None:
                out[i, j] = float(val)
    return out


def feature_catalog(module_key: str) -> dict[str, list[str]]:
    path = CACHE_ROOT / {
        "whole_brain": "atlas.json",
        "microglia": "microglia.json",
        "astrocyte": "astrocyte.json",
        "oligo_opc": "oligo.json",
        "neuron": "neuron.json",
    }[module_key]
    if not path.exists():
        return {"rna": [], "atac": []}
    obj = json.loads(path.read_text(encoding="utf-8"))
    return obj.get("feature_catalog", {"rna": [], "atac": []})


def select_features(module: ModuleConfig, catalog: dict[str, list[str]], max_rna: int, max_atac: int) -> tuple[list[str], list[str]]:
    rna = []
    for gene in PRIORITY_GENES + list(module.marker_seeds) + catalog.get("rna", []):
        if gene not in rna:
            rna.append(gene)
    atac = []
    for peak in catalog.get("atac", []):
        if peak not in atac:
            atac.append(peak)
    return rna[:max_rna], atac[:max_atac]


def summarize_feature(values: np.ndarray, ftype: str, gene: str | None = None) -> dict[str, object]:
    flat = values.astype(np.float32)
    finite = flat[np.isfinite(flat)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=np.float32)
    payload: dict[str, object] = {
        "type": ftype,
        "values": np.round(flat, 4).tolist(),
        "min": round(float(finite.min()), 4),
        "max": round(float(finite.max()), 4),
        "q99": round(float(np.quantile(finite, 0.99)), 4),
    }
    if gene:
        payload["gene"] = gene
    return payload


def subtype_markers(subtypes: list[str], rna_features: dict[str, dict[str, object]], min_genes: int = 5, max_genes: int = 8) -> dict[str, list[str]]:
    if not rna_features:
        return {}
    genes = list(rna_features)
    matrix = np.vstack([np.asarray(rna_features[g]["values"], dtype=np.float32) for g in genes]).T
    out: dict[str, list[str]] = {}
    subtype_arr = np.asarray(subtypes, dtype=object)
    for subtype in sorted(set(subtypes)):
        mask = subtype_arr == subtype
        if not mask.any():
            continue
        mean_vals = matrix[mask].mean(axis=0)
        ranked = [genes[i] for i in np.argsort(-mean_vals)]
        filtered = [g for g in ranked if float(mean_vals[genes.index(g)]) > 0][:max_genes]
        if len(filtered) < min_genes:
            filtered = ranked[:min_genes]
        out[subtype] = filtered[:max_genes]
    return out


def export_module(module: ModuleConfig, args: argparse.Namespace) -> dict[str, object]:
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
    catalog = feature_catalog(module.key)
    selected_rna, selected_atac = select_features(module, catalog, args.max_rna_features, args.max_atac_features)
    max_cells = capped_cells(args.max_cells_per_module, len(selected_rna), len(selected_atac), args.target_data_mb // len(MODULES))
    sample_idx = stratified_sample(meta, subtype_field, max_cells, args.seed)
    sampled_meta = meta.iloc[sample_idx].copy()
    display_meta = normalize_display_meta(sampled_meta, subtype_field)

    cells = []
    for pos, row in enumerate(sample_idx):
        rec = [round(float(embedding[row, 0]), 4), round(float(embedding[row, 1]), 4)]
        for key in ["subtype", "disease", "sample", "dataset"]:
            rec.append(clean_value(display_meta.iloc[pos][key]))
        cells.append(rec)

    rna_payload: dict[str, dict[str, object]] = {}
    atac_payload: dict[str, dict[str, object]] = {}
    if module.h5mu_path.exists():
        try:
            with h5py.File(module.h5mu_path, "r") as handle:
                rna_names = get_var_names(handle, "rna")
                atac_names = get_var_names(handle, "atac")
                rna_lookup = {name.upper(): i for i, name in enumerate(rna_names)}
                atac_lookup = {name: i for i, name in enumerate(atac_names)}
                kept_rna = [g for g in selected_rna if g.upper() in rna_lookup]
                kept_atac = [p for p in selected_atac if p in atac_lookup]
                if kept_rna:
                    rna_mat = extract_sparse_rows(handle["mod"]["rna"]["X"], sample_idx, [rna_lookup[g.upper()] for g in kept_rna])
                    for j, gene in enumerate(kept_rna):
                        rna_payload[gene] = summarize_feature(rna_mat[:, j], "rna")
                if kept_atac:
                    atac_mat = extract_sparse_rows(handle["mod"]["atac"]["X"], sample_idx, [atac_lookup[p] for p in kept_atac])
                    for j, peak in enumerate(kept_atac):
                        atac_payload[peak] = summarize_feature(atac_mat[:, j], "atac")
        except Exception as exc:
            warnings.append(f"Feature export skipped: {type(exc).__name__}: {exc}")
    else:
        warnings.append("h5mu file missing; only UMAP and metadata exported.")

    subtype_list = [clean_value(x) for x in display_meta["subtype"].tolist()]
    disease_list = [clean_value(x) for x in display_meta["disease"].tolist()]
    markers = {
        "subtype_markers": subtype_markers(subtype_list, rna_payload),
        "priority_genes": PRIORITY_GENES,
    }
    categories = {
        "subtype_field": subtype_field,
        "available_color_by": ["subtype", "disease", "sample", "dataset", "selected_feature"],
        "subtype_levels": sorted(set(subtype_list)),
        "disease_levels": sorted(set(disease_list)),
        "sample_levels": sorted(set(display_meta["sample"].tolist())),
        "dataset_levels": sorted(set(display_meta["dataset"].tolist())),
    }

    source_subtype_counts = Counter(meta[subtype_field].map(clean_value))
    disease_col = "disease_new" if "disease_new" in meta.columns else "DISEASE_NEW" if "DISEASE_NEW" in meta.columns else None
    disease_series = meta[disease_col].map(clean_value) if disease_col else pd.Series(["NA"] * len(meta), index=meta.index)
    source_disease_counts = Counter(disease_series)
    contingency = (
        pd.crosstab(meta[subtype_field].map(clean_value), disease_series)
        .astype(int)
        .to_dict(orient="index")
    )
    summary = {
        "module": module.key,
        "label": module.label,
        "subtype_field": subtype_field,
        "n_total_source_cells": int(n),
        "n_exported_cells": int(len(sample_idx)),
        "n_subtypes": int(len(source_subtype_counts)),
        "n_diseases": int(len(source_disease_counts)),
        "source_subtype_counts": dict(source_subtype_counts.most_common()),
        "source_disease_counts": dict(source_disease_counts.most_common()),
        "subtype_by_disease": contingency,
        "export_parameters": {
            "max_cells_per_module": args.max_cells_per_module,
            "max_rna_features": args.max_rna_features,
            "max_atac_features": args.max_atac_features,
            "target_data_mb": args.target_data_mb,
            "seed": args.seed,
        },
        "export_time_utc": datetime.now(timezone.utc).isoformat(),
        "blind_review_note": "Static blind-review demo with downsampled cells and selected RNA/ATAC features only.",
        "warnings": warnings,
    }

    write_json(out_dir / "cells.json", {
        "module": module.key,
        "subtype_field": subtype_field,
        "n_total_source_cells": int(n),
        "n_exported_cells": int(len(sample_idx)),
        "columns": ["x", "y", "subtype", "disease", "sample", "dataset"],
        "cells": cells,
    })
    write_json(out_dir / "categories.json", categories)
    write_json(out_dir / "metadata_summary.json", summary)
    write_json(out_dir / "rna_features.json", {"features": rna_payload})
    write_json(out_dir / "atac_features.json", {"features": atac_payload})
    write_json(out_dir / "markers.json", markers)

    size_mb = sum(p.stat().st_size for p in out_dir.glob("*.json")) / 1024 / 1024
    return {
        "module": module.key,
        "label": module.label,
        "path": f"data/{module.key}",
        "default_subtype_field": subtype_field,
        "default_color_by": "subtype",
        "n_exported_cells": int(len(sample_idx)),
        "n_total_source_cells": int(n),
        "rna_feature_count": len(rna_payload),
        "atac_feature_count": len(atac_payload),
        "size_mb": round(size_mb, 3),
        "warnings": warnings,
    }


def export_reference_mapping(out_root: Path, manifest_modules: list[dict[str, object]]) -> None:
    ref_dir = out_root / "reference_mapping"
    summary = {
        "title": "Reference Mapping Demo",
        "description": "Static summary of the reference mapping workflow used for blind-review visualization.",
        "workflow": [
            "Load harmonized embedding and subtype-level metadata.",
            "Project query cells into reference UMAP coordinates.",
            "Inspect transferred labels together with RNA and ATAC demo features.",
            "Review subtype-level agreement and disease composition in static tables.",
        ],
        "modules": [{k: m[k] for k in ["module", "label", "n_exported_cells", "n_total_source_cells"]} for m in manifest_modules],
    }
    example = {
        "columns": ["query_module", "reference_view", "default_label_field", "n_exported_cells", "note"],
        "rows": [
            [m["label"], "Shared UMAP", m["default_subtype_field"], m["n_exported_cells"], "Downsampled demo view"]
            for m in manifest_modules
        ],
    }
    write_json(ref_dir / "summary.json", summary)
    write_json(ref_dir / "example_mapping.json", example)


def build_manifest(out_root: Path, modules: list[dict[str, object]]) -> None:
    manifest = {
        "title": "Interactive Multi-omics Brain Atlas Demo",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "modules": modules,
        "static_note": "Blind-review static demo. Data are downsampled for interactive visualization.",
    }
    write_json(out_root / "manifest.json", manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export static atlas demo JSON files for GitHub Pages.")
    parser.add_argument("--module", default="all", choices=["all", *MODULES.keys()])
    parser.add_argument("--max-cells-per-module", type=int, default=25000)
    parser.add_argument("--max-rna-features", type=int, default=48)
    parser.add_argument("--max-atac-features", type=int, default=20)
    parser.add_argument("--target-data-mb", type=int, default=350)
    parser.add_argument("--seed", type=int, default=20260626)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    selected = MODULES.values() if args.module == "all" else [MODULES[args.module]]
    manifest_modules = []
    for module in selected:
        manifest_modules.append(export_module(module, args))
    if args.module == "all":
        export_reference_mapping(out_root, manifest_modules)
        build_manifest(out_root, manifest_modules)
    log("Export complete")


if __name__ == "__main__":
    main()
