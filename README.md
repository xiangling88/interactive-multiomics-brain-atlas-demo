# Interactive Multi-omics Brain Atlas Demo

This repository is a blind-review static demo for a multi-omics brain atlas browser. It contains only downsampled JSON files for GitHub Pages deployment and does not include raw large-scale source data.

The static site lives in `docs/` and reads only:

- `docs/index.html`
- `docs/app.js`
- `docs/style.css`
- `docs/data/**/*.json`

No Python backend, `/api/*` endpoint, `h5mu`, full `npz`, full metadata table, or large SCARlink cache is required at deploy time.

## Local export

Recommended dependencies:

- Python 3.10+
- `numpy`
- `pandas`
- `h5py`
- `anndata`
- `mudata`

Export atlas demo data:

```bash
python scripts/export_atlas_demo.py --module all --out docs/data --overwrite
```

Export SCARlink demo data:

```bash
python scripts/export_scarlink_demo.py --out docs/data/scarlink --overwrite
```

One-command build:

```bash
bash scripts/build_demo_all.sh
```

## Local preview

```bash
cd docs
python -m http.server 8000
```

Then open `http://localhost:8000`.

## GitHub Pages

1. Open repository `Settings`.
2. Go to `Pages`.
3. Under `Build and deployment`, select `Deploy from a branch`.
4. Choose branch `main`.
5. Choose folder `/docs`.

The page loads Plotly from a public CDN. If a fully offline deployment is required, vendor Plotly locally into `docs/`.

## Updating the demo

1. Re-run the export scripts.
2. Review the generated `docs/data/**/*.json`.
3. Confirm no local path or identifying metadata leaked into `docs/`, `README.md`, or commit content.
4. Commit only the small static demo assets.

## File size policy

Do not commit:

- raw `h5mu` / `h5ad`
- full `npz`
- full metadata text files
- compressed source matrices
- the original SCARlink result directory
- large cache artifacts not required by `docs/`

Only the downsampled JSON files inside `docs/data/` should be committed for the web demo.

## Blind-review notes

- Keep page copy generic.
- Do not expose author names, institution names, local servers, usernames, or absolute paths.
- Re-check `docs/`, `README.md`, and git history before pushing.
