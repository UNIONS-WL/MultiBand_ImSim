#!/usr/bin/env bash
set -euo pipefail

DIR="/home/hervas/n25/SP_simu_fab/SP_1p2z/input_tiles"
OUTFILE="/home/hervas/n25/SP_simu_fab/tile_numbers_real.txt"

cd "$DIR"

# Extract tile numbers from image files only
for f in CFIS_simu_image-*.fits; do
    # Remove prefix and suffix
    tile=${f#CFIS_simu_image-}
    tile=${tile%.fits}
    echo "$tile"
done | sort -u > "$OUTFILE"

echo "Tile numbers written to $DIR/$OUTFILE"

