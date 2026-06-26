#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

APP_ROOT = Path(__file__).resolve().parents[1]
WORK_ROOT = APP_ROOT.parents[2]
GENE_REF = APP_ROOT / "hg38_new_web.txt"
DEFAULT_OUT = APP_ROOT / "docs" / "data" / "scarlink"
DEFAULT_SCARLINK_DIR = WORK_ROOT / "SCARlink" / "web_add_24_genes"
DEFAULT_GENES = ["CDH4", "GRM3", "ADAMTS17", "CNTNAP2", "ZFHX4"]


def log(msg: str) -> None:
    print(f"[export_scarlink_demo] {msg}")


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")


def clean_text(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    text = str(value).strip()
    return text if text else "NA"


def slugify(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "NA"


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename = {c: c.lower().replace(" ", "_").replace("-", "_") for c in df.columns}
    return df.rename(columns=rename)


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
    out = {}
    for _, row in df.iterrows():
        gene = clean_text(row[gene_col]).upper()
        strand = clean_text(row[strand_col]) if strand_col else "+"
        start = int(row[start_col])
        end = int(row[end_col])
        out[gene] = {"chr": clean_text(row[chr_col]), "tss": start if strand == "+" else end, "start": start, "end": end}
    return out


def load_frames(scarlink_dir: Path) -> list[tuple[str, pd.DataFrame]]:
    suffix = "_gene_linked_tiles_celltype_r2.csv.gz"
    frames = []
    for path in sorted(scarlink_dir.glob(f"*{suffix}")):
        disease = path.name.removesuffix(suffix).replace("_", " ")
        frames.append((disease, standardize_columns(pd.read_csv(path, sep="\t", compression="gzip", low_memory=False))))
    return frames


def choose_genes(frames: list[tuple[str, pd.DataFrame]], requested: list[str], max_genes: int = 5) -> list[str]:
    present = Counter()
    for _, df in frames:
        if "gene" in df.columns:
            present.update(g.upper() for g in df["gene"].dropna().astype(str).unique().tolist())
    selected = [g for g in requested if g.upper() in present]
    if len(selected) < max_genes:
        selected.extend([g for g, _ in present.most_common() if g not in {x.upper() for x in selected}][: max_genes - len(selected)])
    return selected[:max_genes]


def build_payload(gene: str, disease: str, df: pd.DataFrame, gene_ref: dict[str, dict[str, object]]) -> dict[str, object]:
    sub = df[df["gene"].astype(str).str.upper() == gene.upper()].copy()
    if sub.empty:
        return {"gene": gene, "disease": disease, "links": [], "table": [], "box": {}, "message": f"No rows matched {gene} in {disease}."}
    sub = sub.sort_values(["fdr", "p_value"], na_position="last").head(220)
    ref = gene_ref.get(gene.upper(), {})
    q = sub.iloc[0]
    chr_name = clean_text(q.get("chr", ref.get("chr", "chrNA")))
    q_start = int(q.get("start", ref.get("start", 0)))
    q_end = int(q.get("end", ref.get("end", 0)))
    tss = int(ref.get("tss", q_start))
    rows = []
    grouped = defaultdict(list)
    for _, row in sub.iterrows():
        peak = f"{clean_text(row.get('chr'))}:{int(row.get('start'))}-{int(row.get('end'))}"
        z = round(float(row.get("z_score", 0.0)), 6)
        grouped[clean_text(row.get("celltype_r2"))].append(z)
        rows.append({
            "disease": disease,
            "celltype_r2": clean_text(row.get("celltype_r2")),
            "gene": gene,
            "peak": peak,
            "region": peak,
            "regression_coef": round(float(row.get("regression_coef", 0.0)), 6),
            "pval": round(float(row.get("p_value", 1.0)), 6),
            "fdr": round(float(row.get("fdr", 1.0)), 6),
            "z_score": z,
            "spearman_corr": round(float(row.get("spearman_corr", 0.0)), 6),
        })
    links = []
    for row in rows[:72]:
        start = int(row["peak"].split(":")[1].split("-")[0])
        end = int(row["peak"].split("-")[1])
        links.append({
            "gene": gene,
            "enhancer_chr": chr_name,
            "enhancer_start": start,
            "enhancer_end": end,
            "promoter_start": max(1, tss - 2000),
            "promoter_end": tss + 2000,
            "tss": tss,
            "regression_coef": row["regression_coef"],
            "fdr": row["fdr"],
            "effect": "activation" if row["regression_coef"] >= 0 else "repression",
        })
    box = {k: v[:500] for k, v in sorted(grouped.items(), key=lambda kv: -abs(sum(kv[1]) / max(len(kv[1]), 1)))[:28]}
    return {
        "gene": gene,
        "disease": disease,
        "summary": {"n_rows": len(rows), "query_region": f"{chr_name}:{q_start}-{q_end}"},
        "query": {"chr": chr_name, "start": q_start, "end": q_end},
        "links": links,
        "table": rows,
        "box": box,
        "message": f"SCARlink links for {gene} in {disease}.",
    }


def export_selected_genes(genes: list[str], scarlink_dir: Path, out_dir: Path) -> dict[str, object]:
    frames = load_frames(scarlink_dir)
    gene_ref = load_gene_reference(GENE_REF)
    manifest = []
    disease_names = [d for d, _ in frames]
    for gene in genes:
        gene_entry = {"gene": gene, "diseases": []}
        for disease, df in frames:
            payload = build_payload(gene, disease, df, gene_ref)
            disease_slug = slugify(disease)
            rel = f"data/scarlink/{gene}/{disease_slug}.json"
            write_json(out_dir / gene / f"{disease_slug}.json", payload)
            gene_entry["diseases"].append({"name": disease, "slug": disease_slug, "path": rel, "n_rows": payload.get("summary", {}).get("n_rows", 0)})
        manifest.append(gene_entry)
    write_json(out_dir / "scarlink_manifest.json", {
        "title": "SCARlink links",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "genes": manifest,
        "diseases": disease_names,
    })
    write_json(out_dir / "table_summary.json", {
        "columns": ["gene", "disease", "n_rows"],
        "rows": [[g["gene"], d["name"], d["n_rows"]] for g in manifest for d in g["diseases"]],
    })
    return {"genes": manifest, "diseases": disease_names}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export SCARlink demo JSON files.")
    parser.add_argument("--scarlink-dir", default=str(DEFAULT_SCARLINK_DIR))
    parser.add_argument("--cache-dir", default=str(APP_ROOT / "atlas_light_cache"))
    parser.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = load_frames(Path(args.scarlink_dir))
    genes = choose_genes(frames, args.genes, max_genes=5)
    export_selected_genes(genes, Path(args.scarlink_dir), out_dir)
    log(f"Exported {len(genes)} SCARlink genes")


if __name__ == "__main__":
    main()
