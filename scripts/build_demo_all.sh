#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MAX_CELLS_PER_MODULE="${MAX_CELLS_PER_MODULE:-25000}"
TARGET_DATA_MB="${TARGET_DATA_MB:-350}"
OVERWRITE="${OVERWRITE:-1}"

ATLAS_ARGS=(
  --module all
  --max-cells-per-module "${MAX_CELLS_PER_MODULE}"
  --target-data-mb "${TARGET_DATA_MB}"
  --out "${ROOT}/docs/data"
)
SCARLINK_ARGS=(
  --out "${ROOT}/docs/data/scarlink"
)

if [[ "${OVERWRITE}" == "1" ]]; then
  ATLAS_ARGS+=(--overwrite)
  SCARLINK_ARGS+=(--overwrite)
fi

python "${ROOT}/scripts/export_atlas_demo.py" "${ATLAS_ARGS[@]}"
python "${ROOT}/scripts/export_scarlink_demo.py" "${SCARLINK_ARGS[@]}"
