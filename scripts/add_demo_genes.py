#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

import export_atlas_demo as atlas
import export_scarlink_demo as scar

APP_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = APP_ROOT / "docs" / "data"


def log(msg: str) -> None:
    print(f"[add_demo_genes] {msg}")


def read_genes(args: argparse.Namespace) -> list[str]:
    genes = []
    for gene in args.genes or []:
        if gene not in genes:
            genes.append(gene)
    if args.gene_file:
        path = Path(args.gene_file)
        if path.exists():
            for line in path.read_text(encoding="utf-8").splitlines():
                text = line.strip()
                if text and not text.startswith("#") and text not in genes:
                    genes.append(text)
    return genes


def update_module_features(module_key: str, genes: list[str], out_root: Path, update_rna: bool, update_atac: bool, feature_storage: str, dry_run: bool) -> None:
    module = atlas.MODULES[module_key]
    module_dir = out_root / module_key
    idx_path = module_dir / "selected_indices.json"
    manifest_path = out_root / "manifest.json"
    if not idx_path.exists():
        log(f"warning: {module_key} has no selected_indices.json; skip")
        return
    exported_indices = np.asarray(json.loads(idx_path.read_text())["indices"], dtype=np.int64)
    meta = atlas.read_meta(module.meta_path)
    embedding = atlas.load_embedding(module.embedding_path)
    n = min(len(meta), embedding.shape[0])
    meta = meta.iloc[:n].copy()
    subtype_field = atlas.choose_subtype_field(meta, module)
    meta_norm = atlas.normalize_module_meta(meta, module, subtype_field)
    if dry_run:
        log(f"dry-run: would update {module_key} for {', '.join(genes)}")
        return
    with h5py.File(module.h5mu_path, "r") as handle:
        rna_names = atlas.get_var_names(handle, "rna")
        atac_names = atlas.get_var_names(handle, "atac")
        rna_lookup = {x.upper(): i for i, x in enumerate(rna_names)}
        atac_lookup = {x: i for i, x in enumerate(atac_names)}
        feature_files = {"rna": {}, "atac": {}}
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            for mod in manifest.get("modules", []):
                if mod.get("module") == module_key:
                    feature_files = mod.get("feature_files", feature_files)
                    break
        if update_rna:
            keep = [g for g in genes if g.upper() in rna_lookup]
            if keep:
                mat = atlas.extract_sparse_rows(handle["mod"]["rna"]["X"], exported_indices, [rna_lookup[g.upper()] for g in keep])
                for j, gene in enumerate(keep):
                    payload = {
                        "feature": gene,
                        "label": gene,
                        "type": "rna",
                        **atlas.encode_values(mat[:, j], feature_storage),
                        "violin": {
                            "subtype": atlas.build_group_payload(mat[:, j], meta_norm.iloc[exported_indices].reset_index(drop=True), "subtype"),
                            "disease": atlas.build_group_payload(mat[:, j], meta_norm.iloc[exported_indices].reset_index(drop=True), "disease"),
                            "subtype_disease": atlas.build_group_payload(mat[:, j], meta_norm.iloc[exported_indices].reset_index(drop=True), "subtype_disease"),
                        },
                    }
                    safe = atlas.slugify(gene)
                    rel = f"data/{module_key}/features/rna/{safe}.json"
                    atlas.write_json(module_dir / "features" / "rna" / f"{safe}.json", payload)
                    feature_files.setdefault("rna", {})[gene] = rel
            missing = [g for g in genes if g.upper() not in rna_lookup]
            for gene in missing:
                log(f"warning: {gene} not found in RNA for {module_key}")
        if update_atac:
            for gene in genes:
                peak = next((x for x in atac_names if gene.upper() in x.upper()), None)
                if not peak:
                    log(f"warning: no ATAC-like feature found for {gene} in {module_key}")
                    continue
                mat = atlas.extract_sparse_rows(handle["mod"]["atac"]["X"], exported_indices, [atac_lookup[peak]])
                payload = {
                    "feature": peak,
                    "label": peak,
                    "type": "atac",
                    **atlas.encode_values(mat[:, 0], feature_storage),
                    "violin": {
                        "subtype": atlas.build_group_payload(mat[:, 0], meta_norm.iloc[exported_indices].reset_index(drop=True), "subtype"),
                        "disease": atlas.build_group_payload(mat[:, 0], meta_norm.iloc[exported_indices].reset_index(drop=True), "disease"),
                        "subtype_disease": atlas.build_group_payload(mat[:, 0], meta_norm.iloc[exported_indices].reset_index(drop=True), "subtype_disease"),
                    },
                }
                safe = atlas.slugify(peak)
                rel = f"data/{module_key}/features/atac/{safe}.json"
                atlas.write_json(module_dir / "features" / "atac" / f"{safe}.json", payload)
                feature_files.setdefault("atac", {})[peak] = rel
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            for mod in manifest.get("modules", []):
                if mod.get("module") == module_key:
                    mod["feature_files"] = feature_files
                    mod["features"] = {k: list(v) for k, v in feature_files.items()}
            atlas.write_json(manifest_path, manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Append selected genes to exported static demo features.")
    parser.add_argument("--genes", nargs="*")
    parser.add_argument("--gene-file")
    parser.add_argument("--modules", default="all")
    parser.add_argument("--update-scarlink", action="store_true")
    parser.add_argument("--update-rna", action="store_true")
    parser.add_argument("--update-atac", action="store_true")
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--feature-storage", choices=["auto", "dense", "sparse", "quantized"], default="auto")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genes = read_genes(args)
    if not genes:
        raise SystemExit("No genes provided")
    modules = list(atlas.MODULES) if args.modules == "all" else [x.strip() for x in args.modules.split(",") if x.strip()]
    out_root = Path(args.out)
    for module_key in modules:
        update_module_features(module_key, genes, out_root, args.update_rna, args.update_atac, args.feature_storage, args.dry_run)
    if args.update_scarlink:
        if args.dry_run:
            log(f"dry-run: would update SCARlink for {', '.join(genes)}")
        else:
            frames = scar.load_frames(scar.DEFAULT_SCARLINK_DIR)
            chosen = scar.choose_genes(frames, genes, max_genes=len(genes))
            missing = [g for g in genes if g.upper() not in {x.upper() for x in chosen}]
            for gene in missing:
                log(f"warning: no SCARlink data for {gene}")
            scar.export_selected_genes(chosen, scar.DEFAULT_SCARLINK_DIR, out_root / "scarlink")
    log("Done")


if __name__ == "__main__":
    main()
