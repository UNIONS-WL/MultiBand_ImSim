#!/usr/bin/env bash
set -euo pipefail

cd /home/hervas/n25/SP_simu_fab/

pairs=(
  "SP_1z2p_small_low_SNR_grid SP_1z2p_grid"
  "SP_1z2m_small_low_SNR_grid SP_1z2m_grid"
  "SP_1z2z_small_low_SNR_grid SP_1z2z_grid"
  "SP_1p2z_small_low_SNR_grid SP_1p2z_grid"
  "SP_1m2z_small_low_SNR_grid SP_1m2z_grid"
)

for pair in "${pairs[@]}"; do
  set -- $pair
  src=$1
  dst=$2

  echo "==> Copying: $src -> $dst"

  if [[ ! -d "$src" ]]; then
    echo "ERROR: source directory not found: $src" >&2
    exit 1
  fi

  rm -rf -- "$dst"
  mkdir -p -- "$dst"

  rsync -a \
    --exclude='outputs/' \
    --exclude='output/' \
    --exclude='input_tiles' \
    --exclude='input_exp' \
    "$src"/ "$dst"/

  echo "    Done."
done

echo "All copies completed."

