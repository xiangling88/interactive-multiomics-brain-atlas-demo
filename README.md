# Interactive Multi-omics Brain Atlas

This repository hosts a blind-review static atlas browser for GitHub Pages. The deployed site is built from `docs/` only and does not require any Python backend or `/api/*` endpoint.

## Static site contents

The published browser reads only:

- `docs/index.html`
- `docs/app.js`
- `docs/style.css`
- `docs/data/**/*.json`

Large raw files such as full `h5mu`, `npz`, full metadata tables, compressed caches, and original SCARlink result archives are not committed for web deployment.

## Export atlas data

Full chunked embedding export with per-feature lazy-loading:

```bash
python scripts/export_atlas_demo.py --module all --full-embedding --feature-per-file --out docs/data --overwrite
```

Useful size controls:

```bash
python scripts/export_atlas_demo.py \
  --module all \
  --full-embedding \
  --feature-per-file \
  --max-total-docs-mb 800 \
  --max-json-mb 45 \
  --cell-chunk-size 50000 \
  --out docs/data \
  --overwrite
```

If the estimated total size exceeds the configured limit, the exporter automatically falls back to larger stratified subsamples instead of writing oversized payloads.

## Export SCARlink demo data

```bash
python scripts/export_scarlink_demo.py --out docs/data/scarlink --overwrite
```

## Add more display genes later

Append selected genes to existing exported module features and optional SCARlink examples:

```bash
python scripts/add_demo_genes.py --genes KLF12 CDH4 --modules astrocyte --update-rna --update-atac --update-scarlink
```

You can also read genes from a file:

```bash
python scripts/add_demo_genes.py --gene-file genes.txt --modules all --update-rna --update-atac
```

Dry-run is supported:

```bash
python scripts/add_demo_genes.py --genes KLF12 CDH4 --modules astrocyte --update-rna --dry-run
```

## Local preview

```bash
cd docs && python -m http.server 8000
```

Then open `http://localhost:8000`.

## GitHub Pages update workflow

After regenerating `docs/data`:

```bash
git add docs scripts README.md .gitignore
git commit -m "Improve full atlas static demo"
git push origin main
```

In the repository settings, configure GitHub Pages to publish from branch `main` and folder `/docs`.

## File size policy

Do not commit:

- raw `h5mu` / `h5ad`
- full `npz`
- full metadata text tables
- compressed source matrices
- original SCARlink result directories
- large cache artifacts not needed by `docs/`

Commit only the downsampled static JSON files required by the published browser.

## Blind-review notes

- Keep page copy generic.
- Do not expose local paths, usernames, author names, institution names, or server identifiers in `docs/` or `README.md`.
- Re-run the content scans before pushing updates.
