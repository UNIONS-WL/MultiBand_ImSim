#!/usr/bin/env bash
set -euo pipefail

# --- User inputs ---
TILES_FILE="${1:-tile_numbers_real.txt}"
BASE="${2:-/n09data/hervas/SP_simu_fab}"

# The 4 shear combinations you mentioned
SHEARS=(1p2z 1m2z 1z2p 1z2m)

# --- Helpers ---
die() { echo "ERROR: $*" >&2; exit 1; }

[[ -f "$TILES_FILE" ]] || die "Tiles file not found: $TILES_FILE"
[[ -d "$BASE" ]] || die "Base directory not found: $BASE"

# Read tiles, ignoring blank lines / comments
mapfile -t TILES < <(grep -E '^[0-9]{3}-[0-9]{3}$' "$TILES_FILE" || true)
(( ${#TILES[@]} > 0 )) || die "No tiles found in $TILES_FILE (expected lines like 278-281)"

echo "BASE=$BASE"
echo "TILES_FILE=$TILES_FILE"
echo "Tiles: ${#TILES[@]}"

for shear in "${SHEARS[@]}"; do
  # Match your pattern: SP_YYYY where YYYY is shear combo; allow suffix like _grid
  # Prefer exact SP_${shear}* directory; if multiple, loop them
  shopt -s nullglob
  sp_dirs=( "$BASE/SP_${shear}"* )
  shopt -u nullglob

  if (( ${#sp_dirs[@]} == 0 )); then
    echo "WARN: No SP dirs found for shear '$shear' under $BASE (pattern: SP_${shear}*)"
    continue
  fi

  for sp_dir in "${sp_dirs[@]}"; do
    [[ -d "$sp_dir" ]] || continue

    input_tiles_dir="$sp_dir/input_tiles"
    [[ -d "$input_tiles_dir" ]] || {
      echo "WARN: Missing input_tiles dir: $input_tiles_dir (skipping this SP dir)"
      continue
    }

    echo "----"
    echo "Processing SP dir: $sp_dir"

    for tile in "${TILES[@]}"; do
      out_tile_dir="$sp_dir/outputs/output_${tile}"
      run_image_dir="$out_tile_dir/run_image1"
      links_dir="$run_image_dir/get_images_runner_run_1/output"

      # Destination links
      dest_img="$links_dir/CFIS_simu_image-${tile}.fits"
      dest_wgt="$links_dir/CFIS_simu_weight-${tile}.fits"

      # Source files (assumed naming matches; adjust here if needed)
      src_img="$input_tiles_dir/CFIS_simu_image-${tile}.fits"
      src_wgt="$input_tiles_dir/CFIS_simu_weight-${tile}.fits"

      # Create needed directories
      mkdir -p "$links_dir"

      # Create/overwrite symlinks if sources exist
      if [[ -f "$src_img" ]]; then
        ln -sfn "$src_img" "$dest_img"
      else
        echo "WARN: Missing source image: $src_img"
      fi

      if [[ -f "$src_wgt" ]]; then
        ln -sfn "$src_wgt" "$dest_wgt"
      else
        echo "WARN: Missing source weight: $src_wgt"
      fi

      # Ensure log_run.txt contains the single required line
      log_file="$out_tile_dir/log_run_sp.txt"
      mkdir -p "$out_tile_dir"

      required_line="${sp_dir}/outputs/output_${tile}/run_image1 get_images_runner"

      if [[ -f "$log_file" ]]; then
        if ! grep -Fqx "$required_line" "$log_file"; then
          echo "$required_line" >> "$log_file"
        fi
      else
        printf "%s\n" "$required_line" > "$log_file"
      fi
    done
  done
done

echo "Done."

