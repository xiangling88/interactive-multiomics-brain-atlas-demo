#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


APP_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = APP_ROOT.parent
ATLAS_ROOT = APP_ROOT.parents[1]

WHOLE_NPZ = SOURCE_ROOT / "embedding" / "fc_cellid_my_umap_embedding_umap.npz"
WHOLE_META = SOURCE_ROOT / "meta" / "fc_cellid_my_umap_embedding_meta.tsv"

ASTRO_NPZ = SOURCE_ROOT / "embedding" / "Astrocytes_0608_rna_atac_midas_umap_embedding_umap.npz"
ASTRO_META = SOURCE_ROOT / "meta" / "Astrocytes_0525_drop_meta_all_new.txt"
ASTRO_CELLS = SOURCE_ROOT / "embedding" / "Astrocytes_0608_rna_atac_midas_umap_embedding.tsv"


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([
        x.decode("utf-8", errors="ignore") if isinstance(x, bytes) else str(x)
        for x in values
    ], dtype=object)


def check_finite(name: str, embedding: np.ndarray) -> None:
    print(f"{name} embedding shape: {embedding.shape}")
    print(f"{name} embedding min: {np.nanmin(embedding, axis=0)}")
    print(f"{name} embedding max: {np.nanmax(embedding, axis=0)}")
    print(f"{name} embedding NaN count: {int(np.isnan(embedding).sum())}")
    print(f"{name} embedding Inf count: {int(np.isinf(embedding).sum())}")
    if not np.isfinite(embedding).all():
        raise SystemExit(f"{name} embedding contains NaN or Inf")


def is_missing_label(value: object) -> bool:
    if value is None:
        return True
    text = str(value).strip().lower()
    return text in {"", "na", "nan", "none", "<na>"}


def validate_whole() -> dict[str, object]:
    with np.load(WHOLE_NPZ, allow_pickle=True) as npz:
        embedding = np.asarray(npz["embedding"], dtype=np.float32)
        npz_cellid = np.asarray(npz["cellid"]) if "cellid" in npz.files else None
        npz_cells = decode(npz["cell_names"][:]) if "cell_names" in npz.files else None
        npz_keys = list(npz.files)
    meta = pd.read_csv(WHOLE_META, sep="\t", dtype=str, low_memory=False)
    if "cell" not in meta.columns or "cellid" not in meta.columns:
        raise SystemExit("Whole meta missing required cell/cellid columns")
    meta_cells = meta["cell"].astype(str).to_numpy()
    if embedding.shape[0] != len(meta_cells):
        raise SystemExit("Whole embedding/meta length mismatch")
    if npz_cellid is None:
        raise SystemExit("Whole npz missing cellid array")
    if not np.array_equal(npz_cellid.astype(str), meta["cellid"].astype(str).to_numpy()):
        raise SystemExit("Whole npz cellid differs from meta cellid")
    if pd.Index(meta_cells).duplicated().any():
        raise SystemExit("Whole meta cell column contains duplicated cell names")
    if npz_cells is not None and not np.array_equal(npz_cells.astype(str), meta_cells):
        raise SystemExit("Whole npz cell_names differ from meta cell")
    if "second_label" not in meta.columns:
        raise SystemExit("Whole meta missing second_label")
    missing_second_label = int(meta["second_label"].map(is_missing_label).sum())
    if missing_second_label:
        counts = meta.loc[meta["second_label"].map(is_missing_label), "celltype_r2"].value_counts(dropna=False).to_dict()
        raise SystemExit(f"Whole second_label still has {missing_second_label} missing values: {counts}")
    check_finite("Whole atlas", embedding)
    print(f"Whole atlas second_label missing count: {missing_second_label}")
    print("[PASS] Whole atlas embedding/meta alignment")
    return {
        "cells": int(len(meta_cells)),
        "embedding_shape": list(map(int, embedding.shape)),
        "meta_cells": int(len(meta_cells)),
        "npz_keys": npz_keys,
        "alignment_percent": 100.0,
    }


def validate_astro() -> dict[str, object]:
    with np.load(ASTRO_NPZ, allow_pickle=True) as npz:
        embedding = np.asarray(npz["embedding"], dtype=np.float32)
        npz_cells = decode(npz["cell_names"][:]) if "cell_names" in npz.files else None
        npz_keys = list(npz.files)
    meta = pd.read_csv(ASTRO_META, sep="\t", dtype=str, low_memory=False)
    umap = pd.read_csv(ASTRO_CELLS, sep="\t", dtype={"cell_id": str})
    if "cell_id" not in meta.columns:
        raise SystemExit("Astrocyte meta missing required cell_id column")
    if "celltype_astro_my" not in meta.columns:
        raise SystemExit("Astrocyte meta missing celltype_astro_my")
    meta_cells = meta["cell_id"].astype(str).to_numpy()
    rds_cells = umap["cell_id"].astype(str).to_numpy()
    if embedding.shape[0] != len(meta_cells):
        raise SystemExit("Astrocyte embedding/meta length mismatch")
    if not np.array_equal(meta_cells, rds_cells):
        raise SystemExit("Astrocyte meta cell_id order differs from RDS midas.umap export")
    if npz_cells is not None and not np.array_equal(npz_cells.astype(str), meta_cells):
        raise SystemExit("Astrocyte npz cell_names differ from meta cell_id")
    check_finite("Astrocyte", embedding)
    print(meta["celltype_astro_my"].value_counts(dropna=False))
    print("[PASS] Astrocyte embedding/meta alignment")
    return {
        "cells": int(len(meta_cells)),
        "embedding_shape": list(map(int, embedding.shape)),
        "meta_cells": int(len(meta_cells)),
        "npz_keys": npz_keys,
        "alignment_percent": 100.0,
    }


def main() -> None:
    report = {
        "whole_atlas": validate_whole(),
        "astrocytes": validate_astro(),
    }
    out = APP_ROOT / "docs" / "data" / "embedding_meta_alignment_validation.json"
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
