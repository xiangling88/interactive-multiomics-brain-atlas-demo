#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import re
import shutil
from collections import Counter
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd


APP_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = APP_ROOT.parent
ATLAS_ROOT = APP_ROOT.parents[1]

WHOLE_H5MU = ATLAS_ROOT / "Downstream_analysis_atlas_1230" / "Analysis" / "PTSD_feature5000_epoch2000_0725_rename.h5mu"
TMP_META = ATLAS_ROOT / "fc_feature_mdata" / "tmp.feather_new.txt"
ASTRO_RDS = ATLAS_ROOT / "fc_feature_mdata" / "Astrocytes" / "Astrocytes_0608_rna_atac.rds"
ASTRO_META_IN = SOURCE_ROOT / "meta" / "Astrocytes_0525_drop_meta_all_new.txt"

WHOLE_NPZ = SOURCE_ROOT / "embedding" / "PTSD_feature5000_epoch2000_0725_rename_X_umap_embedding_umap.npz"
WHOLE_META = SOURCE_ROOT / "meta" / "epoch_2000_mdata_meta_all.txt"
ASTRO_NPZ = SOURCE_ROOT / "embedding" / "Astrocytes_0608_rna_atac_midas_umap_embedding_umap.npz"
ASTRO_TSV = SOURCE_ROOT / "embedding" / "Astrocytes_0608_rna_atac_midas_umap_embedding.tsv"
ASTRO_CELLS = SOURCE_ROOT / "embedding" / "Astrocytes_0608_rna_atac_midas_umap_cells.txt"
ASTRO_META_OUT = ASTRO_META_IN

REPORT = APP_ROOT / "docs" / "data" / "alignment_qc_report.json"
BACKUP_DIR = APP_ROOT / "backup_20260819"

TMP_ADD_COLS = [
    "SUBSET_ID",
    "DONOR_ID",
    "RL6",
    "RL_5",
    "RL_4",
    "RL3_2",
    "RL3_1",
    "RL_3",
    "RL_2",
    "C_NAME",
    "FUNCTION",
    "SUB_REGION_PRECISE",
    "REGION_CLASS",
    "REGION_SPECIFIC",
    "REGION_L4",
    "REGION_L1",
    "REGION_L2",
    "REGION_L3",
    "FUNCTION_L1",
    "UNIQUE_ID",
    "ORDER_ID",
    "L1_CELL_TYPE_NEW",
]


def decode_h5(values: np.ndarray) -> list[str]:
    out = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8", errors="ignore"))
        else:
            out.append(str(value))
    return out


def read_obs_column(obs: h5py.Group, key: str) -> list[object]:
    obj = obs[key]
    if isinstance(obj, h5py.Dataset):
        values = obj[:]
        if values.dtype.kind in {"O", "S"}:
            return decode_h5(values)
        return values.tolist()
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = obj["codes"][:]
        categories = decode_h5(obj["categories"][:])
        return [categories[int(code)] if int(code) >= 0 else "" for code in codes]
    return [""] * len(obs["_index"])


def read_h5mu_obs_and_umap(path: Path) -> tuple[pd.DataFrame, np.ndarray, dict[str, object]]:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        cell_names = np.asarray(decode_h5(obs["_index"][:]), dtype=object)
        columns = [key for key in obs.keys() if key != "_index"]
        meta = pd.DataFrame({key: read_obs_column(obs, key) for key in columns}, index=cell_names)
        embedding_key = "X_umap"
        embedding = np.asarray(handle["obsm"][embedding_key][:], dtype=np.float32)
        structure = {
            "obs_shape": [int(len(cell_names)), int(len(columns))],
            "obs_columns": columns,
            "obs_names_first20": cell_names[:20].tolist(),
            "obsm_keys": list(handle["obsm"].keys()),
            "obsm_shapes": {key: list(handle["obsm"][key].shape) for key in handle["obsm"].keys()},
            "mod_keys": list(handle["mod"].keys()),
            "mod_obsm_keys": {
                mod: list(handle["mod"][mod].get("obsm", {}).keys()) for mod in handle["mod"].keys()
            },
            "embedding_source": f'mdata.obsm["{embedding_key}"]',
        }
    meta.index = meta.index.astype(str)
    return meta, embedding, structure


def non_missing(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    text = str(value).strip()
    return text != "" and text.lower() not in {"nan", "na", "none", "<na>"}


def modal_non_missing(values: Iterable[object]) -> object:
    cleaned = [str(v).strip() for v in values if non_missing(v)]
    if not cleaned:
        return ""
    return Counter(cleaned).most_common(1)[0][0]


def backup_once(paths: list[Path]) -> None:
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            dest = BACKUP_DIR / f"{path.name}.backup_before_20260819"
            if not dest.exists():
                shutil.copy2(path, dest)


def write_compatible_npz(old_npz: Path, out_npz: Path, embedding: np.ndarray, cell_names: np.ndarray) -> list[str]:
    with np.load(old_npz, allow_pickle=True) as old:
        keys = list(old.files)
    payload: dict[str, np.ndarray] = {}
    if "embedding" in keys:
        payload["embedding"] = embedding.astype(np.float32, copy=False)
    else:
        payload[keys[0]] = embedding.astype(np.float32, copy=False)
    if "cell_names" in keys:
        payload["cell_names"] = cell_names.astype(str)
    np.savez_compressed(out_npz, **payload)
    return keys


def build_tmp_subset(tmp: pd.DataFrame, fields: list[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    tmp = tmp.copy()
    tmp["SUBSET_ID_merge"] = tmp["SUBSET_ID"].astype(str).str.strip()
    conflict_counts = {}
    for field in fields:
        if field == "SUBSET_ID":
            continue
        per_subset = tmp.groupby("SUBSET_ID_merge", sort=False)[field].nunique(dropna=True)
        conflict_counts[field] = int((per_subset > 1).sum())
    tmp_subset = tmp.groupby("SUBSET_ID_merge", sort=False)[fields].agg(modal_non_missing).reset_index()
    tmp_subset["SUBSET_ID"] = tmp_subset["SUBSET_ID_merge"]
    return tmp_subset, {
        "subset_rows": int(len(tmp_subset)),
        "tmp_rows": int(len(tmp)),
        "tmp_subset_id_unique": int(tmp["SUBSET_ID_merge"].nunique()),
        "conflicting_subset_counts_by_field": conflict_counts,
    }


def prepare_whole() -> dict[str, object]:
    meta, embedding, structure = read_h5mu_obs_and_umap(WHOLE_H5MU)
    cell_names = meta.index.to_numpy(dtype=str)
    if embedding.shape[0] != len(cell_names):
        raise ValueError(f"Whole embedding rows {embedding.shape[0]} != obs names {len(cell_names)}")
    if embedding.shape[1] != 2:
        raise ValueError(f"Whole embedding must be 2D, got {embedding.shape}")
    if not np.isfinite(embedding).all():
        raise ValueError("Whole embedding contains NaN or Inf")

    header = pd.read_csv(TMP_META, sep="\t", nrows=0).columns.tolist()
    fields = [col for col in TMP_ADD_COLS if col in header]
    tmp = pd.read_csv(TMP_META, sep="\t", usecols=fields, low_memory=False)
    tmp_subset, tmp_report = build_tmp_subset(tmp, fields)

    meta_out = meta.copy()
    meta_out.insert(0, "cell_id", cell_names)
    meta_out["_original_order"] = np.arange(meta_out.shape[0])
    meta_out["dataset_merge"] = meta_out["dataset"].astype(str).str.strip()

    add_cols = ["SUBSET_ID"] + [col for col in fields if col != "SUBSET_ID" and col not in meta_out.columns]
    tmp_add = tmp_subset[["SUBSET_ID_merge", *add_cols]].drop_duplicates("SUBSET_ID_merge")
    meta_new = meta_out.merge(
        tmp_add,
        how="left",
        left_on="dataset_merge",
        right_on="SUBSET_ID_merge",
        validate="many_to_one",
    )
    meta_new = meta_new.sort_values("_original_order")
    if not np.array_equal(meta_new["cell_id"].astype(str).to_numpy(), cell_names):
        raise ValueError("Whole meta cell_id order does not match h5mu obs_names")
    meta_new = meta_new.drop(columns=["_original_order", "dataset_merge", "SUBSET_ID_merge"])

    backup_once([WHOLE_NPZ, WHOLE_META])
    old_npz_keys = write_compatible_npz(WHOLE_NPZ, WHOLE_NPZ, embedding, cell_names)
    meta_new.to_csv(WHOLE_META, sep="\t", index=False)

    dataset_values = meta["dataset"].astype(str).str.strip()
    tmp_values = tmp["SUBSET_ID"].astype(str).str.strip()
    matched_datasets = sorted(set(dataset_values.dropna()) & set(tmp_values.dropna()))
    unmatched_datasets = sorted(set(dataset_values.dropna()) - set(tmp_values.dropna()))
    rl6 = meta_new["RL6"] if "RL6" in meta_new.columns else pd.Series([""] * len(meta_new))

    return {
        "h5mu": str(WHOLE_H5MU),
        "cells": int(len(cell_names)),
        "embedding_source": structure["embedding_source"],
        "embedding_shape": list(map(int, embedding.shape)),
        "embedding_min": np.nanmin(embedding, axis=0).astype(float).tolist(),
        "embedding_max": np.nanmax(embedding, axis=0).astype(float).tolist(),
        "embedding_nan_count": int(np.isnan(embedding).sum()),
        "embedding_inf_count": int(np.isinf(embedding).sum()),
        "meta_cells": int(meta_new.shape[0]),
        "unique_cell_names": int(pd.Index(cell_names).nunique()),
        "duplicated_cell_names": int(pd.Index(cell_names).duplicated().sum()),
        "dataset_unique": int(dataset_values.nunique()),
        "subset_id_matched": int(len(matched_datasets)),
        "subset_id_unmatched": int(len(unmatched_datasets)),
        "unmatched_datasets_first20": unmatched_datasets[:20],
        "rl6_non_null_cells": int(rl6.map(non_missing).sum()),
        "rl6_missing_cells": int((~rl6.map(non_missing)).sum()),
        "old_npz_keys": old_npz_keys,
        "final_embedding": str(WHOLE_NPZ),
        "final_meta": str(WHOLE_META),
        "alignment_percent": 100.0,
        "structure": structure,
        "tmp_subset_report": tmp_report,
    }


def prepare_astro_from_existing_r_exports() -> dict[str, object]:
    umap = pd.read_csv(ASTRO_TSV, sep="\t")
    if not {"cell_id", "UMAP1", "UMAP2"}.issubset(umap.columns):
        raise ValueError(f"Astro UMAP TSV missing required columns: {ASTRO_TSV}")
    meta = pd.read_csv(ASTRO_META_OUT, sep="\t", dtype=str, low_memory=False)
    if "cell_id" not in meta.columns:
        raise ValueError("Astro meta must include cell_id from the R alignment export")
    target = umap["cell_id"].astype(str).to_numpy()
    if not np.array_equal(meta["cell_id"].astype(str).to_numpy(), target):
        raise ValueError("Astro meta cell_id order does not match RDS midas.umap rownames")
    embedding = umap[["UMAP1", "UMAP2"]].to_numpy(dtype=np.float32)
    if not np.isfinite(embedding).all():
        raise ValueError("Astro embedding contains NaN or Inf")

    backup_once([ASTRO_NPZ])
    old_npz_keys = write_compatible_npz(ASTRO_NPZ, ASTRO_NPZ, embedding, target)
    counts = meta["celltype_astro_my"].value_counts(dropna=False).to_dict()
    return {
        "rds": str(ASTRO_RDS),
        "reduction": "midas.umap",
        "rds_cells": int(len(target)),
        "astro_meta_cells": int(len(meta)),
        "matched": int(len(target)),
        "unmatched": 0,
        "celltype_column": "celltype_astro_my",
        "celltype_counts": {str(k): int(v) for k, v in counts.items()},
        "embedding_shape": list(map(int, embedding.shape)),
        "embedding_min": np.nanmin(embedding, axis=0).astype(float).tolist(),
        "embedding_max": np.nanmax(embedding, axis=0).astype(float).tolist(),
        "embedding_nan_count": int(np.isnan(embedding).sum()),
        "embedding_inf_count": int(np.isinf(embedding).sum()),
        "old_npz_keys": old_npz_keys,
        "final_embedding": str(ASTRO_NPZ),
        "final_meta": str(ASTRO_META_OUT),
        "alignment_percent": 100.0,
    }


def main() -> None:
    report = {
        "whole_atlas": prepare_whole(),
        "astrocytes": prepare_astro_from_existing_r_exports(),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print("[PASS] Whole atlas embedding/meta alignment")
    print("[PASS] Astrocyte embedding/meta alignment")


if __name__ == "__main__":
    main()
