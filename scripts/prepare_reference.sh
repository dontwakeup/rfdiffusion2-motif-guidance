#!/usr/bin/env bash
set -euo pipefail

uv run python scripts/prepare_reference.py \
  --motif-json inputs/motif.json \
  --pdb inputs/2LZM.pdb \
  --out inputs/2LZM.npz \
  --summary inputs/reference_summary.json