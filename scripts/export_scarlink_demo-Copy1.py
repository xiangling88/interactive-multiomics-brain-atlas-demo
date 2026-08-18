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
DEFAULT_GENES = [
    "CDH4", "GRM3", "ADAMTS17", "CNTNAP2", "ZFHX4", "SOX6", "GFAP", "AQP4", "P2RY12", "CX3CR1",
    "APOE", "CLU", "MBP", "PLP1", "SLC17A7", "CAMK2A", "RBFOX3", "MAG", "TREM2", "C3",
]


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


def safe_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(out) or math.isinf(out):
        return default
    return out


def slugify(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "NA"


def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename = {c: c.lower().replace(" ", "_").replace("-", "_") for c in df.columns}
    return df.rename(columns=rename)


def load_gene_reference(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="	", low_memory=False)
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
        df = pd.read_csv(path, sep="	", compression="gzip", low_memory=False)
        frames.append((disease, standardize_columns(df)))
    return frames


def choose_genes(frames: list[tuple[str, pd.DataFrame]], requested: list[str], max_genes: int = 20) -> list[str]:
    present = Counter()
    for _, df in frames:
        if "gene" in df.columns:
            present.update(g.upper() for g in df["gene"].dropna().astype(str).unique().tolist())
    selected = []
    for gene in requested:
        if gene.upper() in present and gene.upper() not in {x.upper() for x in selected}:
            selected.append(gene)
    if len(selected) < max_genes:
        selected.extend([g for g, _ in present.most_common() if g not in {x.upper() for x in selected}][: max_genes - len(selected)])
    return selected[:max_genes]


def format_peak(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{start}-{end}"


def fdr_to_score(fdr: float) -> float:
    return round(-math.log10(max(float(fdr), 1e-300)), 6)


def build_payload(gene: str, disease: str, df: pd.DataFrame, gene_ref: dict[str, dict[str, object]]) -> dict[str, object]:
    sub = df[df["gene"].astype(str).str.upper() == gene.upper()].copy()
    if sub.empty:
        return {"gene": gene, "disease": disease, "links": [], "table": [], "box": {}, "message": f"No rows matched {gene} in {disease}."}

    if "fdr" not in sub.columns:
        sub["fdr"] = 1.0
    if "p_value" not in sub.columns:
        sub["p_value"] = 1.0
    if "z_score" not in sub.columns:
        sub["z_score"] = 0.0
    if "spearman_corr" not in sub.columns:
        sub["spearman_corr"] = 0.0

    sub["fdr"] = pd.to_numeric(sub["fdr"], errors="coerce").fillna(1.0)
    sub["p_value"] = pd.to_numeric(sub["p_value"], errors="coerce").fillna(1.0)
    sub["z_score"] = pd.to_numeric(sub["z_score"], errors="coerce").fillna(0.0)
    sub["spearman_corr"] = pd.to_numeric(sub["spearman_corr"], errors="coerce").fillna(0.0)
    sub["regression_coef"] = pd.to_numeric(sub.get("regression_coef", 0.0), errors="coerce").fillna(0.0)
    sub = sub.sort_values(["fdr", "p_value", "z_score"], ascending=[True, True, False], na_position="last").head(320).reset_index(drop=True)

    ref = gene_ref.get(gene.upper(), {})
    q = sub.iloc[0]
    chr_name = clean_text(q.get("chr", ref.get("chr", "chrNA")))
    q_start = int(safe_float(q.get("start", ref.get("start", 0)), 0))
    q_end = int(safe_float(q.get("end", ref.get("end", 0)), 0))
    tss = int(ref.get("tss", q_start))
    promoter_start = max(1, tss - 2000)
    promoter_end = tss + 2000

    grouped = defaultdict(list)
    rows = []
    for rank, (_, row) in enumerate(sub.iterrows(), start=1):
        chrom = clean_text(row.get("chr", chr_name))
        start = int(safe_float(row.get("start", 0), 0))
        end = int(safe_float(row.get("end", 0), 0))
        peak = format_peak(chrom, start, end)
        fdr = round(safe_float(row.get("fdr", 1.0), 1.0), 6)
        pval = round(safe_float(row.get("p_value", 1.0), 1.0), 6)
        coef = round(safe_float(row.get("regression_coef", 0.0), 0.0), 6)
        z = round(safe_float(row.get("z_score", 0.0), 0.0), 6)
        corr = round(safe_float(row.get("spearman_corr", 0.0), 0.0), 6)
        celltype = clean_text(row.get("celltype_r2"))
        grouped[celltype].append(z)
        rows.append({
            "rank": rank,
            "is_top5": rank <= 5,
            "disease": disease,
            "celltype_r2": celltype,
            "gene": gene,
            "chr": chrom,
            "start": start,
            "end": end,
            "peak": peak,
            "region": peak,
            "promoter_start": promoter_start,
            "promoter_end": promoter_end,
            "tss": tss,
            "regression_coef": coef,
            "pval": pval,
            "fdr": fdr,
            "significance": fdr_to_score(fdr),
            "z_score": z,
            "spearman_corr": corr,
            "effect": "activation" if coef >= 0 else "repression",
        })

    links = []
    for row in rows[:120]:
        links.append({
            "rank": row["rank"],
            "is_top5": row["is_top5"],
            "gene": gene,
            "enhancer_chr": row["chr"],
            "enhancer_start": row["start"],
            "enhancer_end": row["end"],
            "enhancer_label": row["peak"],
            "promoter_start": promoter_start,
            "promoter_end": promoter_end,
            "promoter_label": format_peak(chr_name, promoter_start, promoter_end),
            "tss": tss,
            "regression_coef": row["regression_coef"],
            "fdr": row["fdr"],
            "significance": row["significance"],
            "z_score": row["z_score"],
            "celltype_r2": row["celltype_r2"],
            "effect": row["effect"],
        })

    box = {k: v[:600] for k, v in sorted(grouped.items(), key=lambda kv: (-max(abs(x) for x in kv[1]) if kv[1] else 0, kv[0]))[:30]}
    top_links = [row for row in rows[:5]]
    return {
        "gene": gene,
        "disease": disease,
        "summary": {
            "n_rows": len(rows),
            "query_region": format_peak(chr_name, q_start, q_end),
            "display_window": format_peak(chr_name, min([x["start"] for x in rows] + [promoter_start]), max([x["end"] for x in rows] + [promoter_end])),
            "top_fdr": top_links[0]["fdr"] if top_links else 1.0,
            "top_contacts": len(top_links),
        },
        "query": {"chr": chr_name, "start": q_start, "end": q_end, "promoter_start": promoter_start, "promoter_end": promoter_end, "tss": tss},
        "links": links,
        "top_links": top_links,
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
        total_rows = 0
        for disease, df in frames:
            payload = build_payload(gene, disease, df, gene_ref)
            disease_slug = slugify(disease)
            rel = f"data/scarlink/{gene}/{disease_slug}.json"
            write_json(out_dir / gene / f"{disease_slug}.json", payload)
            n_rows = int(payload.get("summary", {}).get("n_rows", 0))
            total_rows += n_rows
            gene_entry["diseases"].append({"name": disease, "slug": disease_slug, "path": rel, "n_rows": n_rows})
        gene_entry["total_rows"] = total_rows
        manifest.append(gene_entry)
    write_json(out_dir / "scarlink_manifest.json", {
        "title": "SCARlink links",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "genes": manifest,
        "diseases": disease_names,
        "default_gene": manifest[0]["gene"] if manifest else None,
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
    genes = choose_genes(frames, args.genes, max_genes=20)
    export_selected_genes(genes, Path(args.scarlink_dir), out_dir)
    log(f"Exported {len(genes)} SCARlink genes")


if __name__ == "__main__":
    main()
