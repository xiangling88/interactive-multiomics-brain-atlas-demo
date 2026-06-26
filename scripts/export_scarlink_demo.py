#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


APP_ROOT = Path(__file__).resolve().parents[1]
GENE_REF = APP_ROOT / "hg38_new_web.txt"
DEFAULT_OUT = APP_ROOT / "docs" / "data" / "scarlink"
DEFAULT_SCARLINK_DIR = Path("/share/home/fxl/SCARlink/web_add_24_genes")
DEFAULT_CACHE_DIR = APP_ROOT / "atlas_light_cache"
DEFAULT_GENES = ["CDH4", "GRM3", "ADAMTS17", "CNTNAP2", "ZFHX4"]


def log(msg: str) -> None:
    print(f"[export_scarlink_demo] {msg}")


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")


def clean_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return str(value).strip()[:120] or "NA"


def load_gene_reference(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", low_memory=False)
    cols = {c.lower(): c for c in df.columns}
    gene_col = next((cols[c] for c in cols if c in {"gene", "gene_name", "symbol"}), None)
    chr_col = next((cols[c] for c in cols if c in {"chr", "chrom", "chromosome"}), None)
    start_col = next((cols[c] for c in cols if "start" in c), None)
    end_col = next((cols[c] for c in cols if "end" in c), None)
    strand_col = next((cols[c] for c in cols if "strand" in c), None)
    if not all([gene_col, chr_col, start_col, end_col]):
        return {}
    ref = {}
    for _, row in df.iterrows():
        gene = clean_text(row[gene_col]).upper()
        strand = clean_text(row[strand_col]) if strand_col else "+"
        start = int(row[start_col])
        end = int(row[end_col])
        tss = start if strand == "+" else end
        ref[gene] = {"chr": clean_text(row[chr_col]), "start": start, "end": end, "tss": tss, "strand": strand}
    return ref


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename = {}
    for col in df.columns:
        lower = col.lower().replace(" ", "_").replace("-", "_")
        rename[col] = lower
    df = df.rename(columns=rename)
    return df


def iter_scarlink_frames(scarlink_dir: Path) -> list[tuple[str, pd.DataFrame]]:
    suffix = "_gene_linked_tiles_celltype_r2.csv.gz"
    frames = []
    for path in sorted(scarlink_dir.glob(f"*{suffix}")):
        disease = path.name.removesuffix(suffix).replace("_", " ")
        df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        frames.append((disease, standardize_columns(df)))
    return frames


def choose_genes(frames: list[tuple[str, pd.DataFrame]], requested: list[str], max_genes: int = 5) -> list[str]:
    present = Counter()
    for _, df in frames:
        if "gene" not in df.columns:
            continue
        present.update(g.upper() for g in df["gene"].dropna().astype(str).unique().tolist())
    out = [g for g in requested if g.upper() in present]
    if len(out) < max_genes:
        extras = [g for g, _ in present.most_common() if g not in {x.upper() for x in out}]
        out.extend(extras[: max_genes - len(out)])
    return out[:max_genes]


def build_gene_payload(gene: str, frames: list[tuple[str, pd.DataFrame]], gene_ref: dict[str, dict[str, object]]) -> dict[str, object]:
    rows = []
    for disease, df in frames:
        if "gene" not in df.columns:
            continue
        sub = df[df["gene"].astype(str).str.upper() == gene.upper()].copy()
        if sub.empty:
            continue
        sub["disease"] = disease
        rows.append(sub)
    merged = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if merged.empty:
        return {"gene": gene, "summary": {"n_rows": 0}, "circle": {"nodes": [], "links": [], "tracks": []}, "boxplot": {"groups": [], "values": []}, "table": []}
    merged = merged.sort_values(["fdr", "p_value"], na_position="last")
    top = merged.head(200).copy()
    ref = gene_ref.get(gene.upper(), {})
    query_row = top.iloc[0]
    chr_name = clean_text(query_row.get("chr", ref.get("chr", "chrNA")))
    q_start = int(query_row.get("start", ref.get("start", 0)))
    q_end = int(query_row.get("end", ref.get("end", 0)))
    tss = int(ref.get("tss", q_start))
    table = []
    for _, row in top.iterrows():
        table.append({
            "disease": clean_text(row.get("disease")),
            "celltype_r2": clean_text(row.get("celltype_r2")),
            "region": f"{clean_text(row.get('chr'))}:{int(row.get('start'))}-{int(row.get('end'))}",
            "peak": f"{clean_text(row.get('chr'))}:{int(row.get('start'))}-{int(row.get('end'))}",
            "gene": gene,
            "regression_coef": round(float(row.get("regression_coef", 0.0)), 6),
            "pval": round(float(row.get("p_value", 1.0)), 6),
            "fdr": round(float(row.get("fdr", 1.0)), 6),
            "z_score": round(float(row.get("z_score", 0.0)), 6),
            "spearman_corr": round(float(row.get("spearman_corr", 0.0)), 6),
        })
    links = []
    for row in table[:72]:
        links.append({
            "gene": gene,
            "enhancer_chr": chr_name,
            "enhancer_start": int(row["region"].split(":")[1].split("-")[0]),
            "enhancer_end": int(row["region"].split("-")[1]),
            "promoter_start": max(1, tss - 2000),
            "promoter_end": tss + 2000,
            "tss": tss,
            "regression_coef": row["regression_coef"],
            "fdr": row["fdr"],
            "effect": "activation" if row["regression_coef"] >= 0 else "repression",
        })
    groups = []
    grouped = {}
    for row in table:
        grouped.setdefault(row["celltype_r2"], []).append(row["z_score"])
    for label, values in sorted(grouped.items(), key=lambda x: -abs(sum(x[1]) / max(len(x[1]), 1)))[:28]:
        groups.append({"label": label, "values": [round(float(v), 6) for v in values]})
    return {
        "gene": gene,
        "summary": {
            "n_rows": len(table),
            "n_diseases": len(set(r["disease"] for r in table)),
            "top_disease": table[0]["disease"] if table else "NA",
            "query_region": f"{chr_name}:{q_start}-{q_end}",
        },
        "message": f"Example SCARlink links for {gene} across {len(set(r['disease'] for r in table))} disease result sets.",
        "query": {"chr": chr_name, "start": q_start, "end": q_end},
        "links": links,
        "rows": table,
        "circle": {"nodes": [], "links": links, "tracks": []},
        "boxplot": {"groups": groups, "values": [g["values"] for g in groups]},
        "table": table,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export SCARlink demo JSON files for the static atlas browser.")
    parser.add_argument("--scarlink-dir", default=str(DEFAULT_SCARLINK_DIR))
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR))
    parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = iter_scarlink_frames(Path(args.scarlink_dir))
    gene_ref = load_gene_reference(GENE_REF)
    genes = choose_genes(frames, args.genes, max_genes=5)
    manifest = []
    table_summary = []
    for gene in genes:
        payload = build_gene_payload(gene, frames, gene_ref)
        write_json(out_dir / f"{gene}.json", payload)
        manifest.append({
            "gene": gene,
            "path": f"data/scarlink/{gene}.json",
            "n_rows": payload["summary"]["n_rows"],
            "query_region": payload["summary"]["query_region"],
        })
        table_summary.append({
            "gene": gene,
            "n_rows": payload["summary"]["n_rows"],
            "n_diseases": payload["summary"]["n_diseases"],
            "top_disease": payload["summary"]["top_disease"],
        })
    write_json(out_dir / "scarlink_manifest.json", {
        "title": "SCARlink links",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "genes": manifest,
    })
    write_json(out_dir / "table_summary.json", {
        "columns": ["gene", "n_rows", "n_diseases", "top_disease"],
        "rows": [[row["gene"], row["n_rows"], row["n_diseases"], row["top_disease"]] for row in table_summary],
    })
    log(f"Exported {len(manifest)} SCARlink examples")


if __name__ == "__main__":
    main()
