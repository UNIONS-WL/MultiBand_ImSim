#!/bin/bash

n_tiles=50 # Total number of tiles to launch
tile_count=0

echo "You cleaned up outputs right?"

for dir in SP_????_star_grid; do
  echo "$dir"
  ls -l numbers_run.txt
  if [ -d "$dir" ]; then
    cd "$dir" || exit 1
    pwd
    apptainer exec --bind /n17data,/n23data1,/n09data /home/hervas/fhervas/shapepipe-1.4 /home/hervas/fhervas/shapepipe_cail/shapepipe/scripts/sh/combine_runs.bash   
    apptainer exec --bind /n17data,/n23data1,/n09data,/automnt/n09data /home/hervas/fhervas/shapepipe-1.4 python /home/hervas/fhervas/shapepipe-1/scripts/python/merge_final_cat.py -i output/run_sp_combined_final/make_catalog_runner/output -p /home/hervas/fhervas/shapepipe-1/example/cfis/final_cat.param -v
     

    cd /home/hervas/n25/SP_simu_fab/ || exit 1
  fi
done

