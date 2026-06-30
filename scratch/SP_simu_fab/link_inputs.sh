#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

BASE="/home/hervas/n25"
SIMU_DIR="$BASE/SP_simu_fab"
SKILLS_OUT="$BASE/skills_out"

shears=(1z2p 1z2m 1z2z 1p2z 1m2z)
runs=(1 2 3 4 5)

cd "$SIMU_DIR"

for sh in "${shears[@]}"; do
  grid_dir="SP_${sh}"

  if [[ ! -d "$grid_dir" ]]; then
    echo "WARNING: missing $grid_dir — skipping"
    continue
  fi

  echo "==> Linking inputs for $grid_dir"

  # Recreate clean link directories
  rm -rf "$grid_dir/input_tiles" "$grid_dir/input_exp"
  mkdir -p "$grid_dir/input_tiles" "$grid_dir/input_exp"

  for r in "${runs[@]}"; do
    src_tiles="$SKILLS_OUT/${sh}_${r}/images/SP_tiles"
    src_exp="$SKILLS_OUT/${sh}_${r}/images/SP_exp"

    if [[ -d "$src_tiles" ]]; then
      for item in "$src_tiles"/*; do
        ln -sfn "$item" "$grid_dir/input_tiles/"
      done
    fi

    if [[ -d "$src_exp" ]]; then
      for item in "$src_exp"/*; do
        ln -sfn "$item" "$grid_dir/input_exp/"
      done
    fi
  done

  echo "    Done."
done

echo "All linking completed."

