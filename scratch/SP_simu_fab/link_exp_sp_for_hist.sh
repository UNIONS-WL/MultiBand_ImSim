#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

BASE="/home/hervas/n25/skills_out"
shears=(1z2p 1z2m 1z2z 1p2z 1m2z)
runs=(1 2 3 4)

echo "==> Cross-linking exposures across runs in: $BASE"
echo "==> Also building run 5 (links-only) for each shear"

for sh in "${shears[@]}"; do
  # Ensure run-5 destination exists (links-only)
  run5_dir="$BASE/${sh}_5/images/SP_exp"
  mkdir -p "$run5_dir"

  # 1) For runs 1-4: link other runs into each run's SP_exp
  for r_dst in "${runs[@]}"; do
    dst_dir="$BASE/${sh}_${r_dst}/images/SP_exp"
    if [[ ! -d "$dst_dir" ]]; then
      echo "  WARNING: missing destination: $dst_dir (skipping)"
      continue
    fi

    echo "  -> ${sh} run ${r_dst}: linking other runs into $dst_dir"

    for r_src in "${runs[@]}"; do
      [[ "$r_src" == "$r_dst" ]] && continue
      src_dir="$BASE/${sh}_${r_src}/images/SP_exp"
      [[ -d "$src_dir" ]] || continue

      for f in "$src_dir"/*; do
        bn="$(basename "$f")"
        # Do not overwrite anything already present
        [[ -e "$dst_dir/$bn" || -L "$dst_dir/$bn" ]] && continue
        ln -s -- "$f" "$dst_dir/$bn"
      done
    done
  done

  # 2) Run 5: link *everything* from runs 1-4 into run5_dir
  echo "  -> ${sh} run 5: linking all runs into $run5_dir"
  for r_src in "${runs[@]}"; do
    src_dir="$BASE/${sh}_${r_src}/images/SP_exp"
    [[ -d "$src_dir" ]] || continue

    for f in "$src_dir"/*; do
      bn="$(basename "$f")"
      # Do not overwrite anything already present
      [[ -e "$run5_dir/$bn" || -L "$run5_dir/$bn" ]] && continue
      ln -s -- "$f" "$run5_dir/$bn"
    done
  done
done

echo "Done."

