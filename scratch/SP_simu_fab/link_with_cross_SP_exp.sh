#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

cd /home/hervas/n25/SP_simu_fab/



echo "All copies completed."

# ---------------------------------------------
# 2) Cross-link exposures across skills_out runs
# ---------------------------------------------
BASE="/home/hervas/n25/skills_out"
shears=(1z2p 1z2m 1z2z 1p2z 1m2z)
runs=(1 2 3 4 5)

echo "==> Cross-linking exposures across runs in: $BASE"

for sh in "${shears[@]}"; do
  for r_dst in "${runs[@]}"; do
    dst_dir="$BASE/${sh}_${r_dst}/images/SP_exp"

    if [[ ! -d "$dst_dir" ]]; then
      echo "  WARNING: missing destination: $dst_dir (skipping)"
      continue
    fi

    echo "  -> ${sh} run ${r_dst}: linking other runs into $dst_dir"

    for r_src in "${runs[@]}"; do
      [[ "$r_src" == "$r_dst" ]] && continue

      src_dir="$BASE/${sh}_grid_${r_src}/images/SP_exp"
      if [[ ! -d "$src_dir" ]]; then
        echo "     NOTE: missing source: $src_dir"
        continue
      fi

      for f in "$src_dir"/*; do
        bn="$(basename "$f")"
        # Don't overwrite anything already present in destination.
        if [[ -e "$dst_dir/$bn" || -L "$dst_dir/$bn" ]]; then
          continue
        fi
        ln -s -- "$f" "$dst_dir/$bn"
      done
    done
  done
done

echo "Exposure cross-linking completed."

