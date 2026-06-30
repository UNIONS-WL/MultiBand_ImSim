#!/bin/bash
stat /automnt/n09data/hervas/SP_simu_fab >/dev/null 2>&1
stat /automnt/n09data >/dev/null 2>&1

n_tiles=1000 # Total number of tiles to launch
tile_count=0

echo "You cleaned up outputs right?"

for dir in SP_????; do
  echo "$dir"
  ls -l numbers_run.txt
  if [ -d "$dir" ]; then
    cd "$dir" || exit 1
    pwd
    ARGUMENTS_FILE="numbers_run.txt"
    ARGS=($(cat "$ARGUMENTS_FILE"))
    echo $ARGS
    for arg in "${ARGS[@]}"; do
      #if [ "$tile_count" -ge "$n_tiles" ]; then
      #  echo "Reached maximum tile limit ($n_tiles)."
      #  break 2
      #fi

      echo "Processing: $arg"
     
      # Remove output file if it exists
      output_path="outputs/output_$arg"
      #echo $output_path
      #output_path="outputs/output_$arg"
      #echo "$output_path"

      #if [ -d "$output_path" ]; then
      #  index=1
      #  new_output_path="${output_path}_$index"
      #  # Find next available version
      #  while [ -d "$new_output_path" ]; do
      #    index=$((index + 1))
      #    new_output_path="${output_path}_$index"
      #  done
     #   echo "Renaming $output_path to $new_output_path"
      #  mv "$output_path" "$new_output_path"

     #fi
      

      # Submit job
      sbatch --output="/home/hervas/n25/SP_simu_fab/output_log_job/${dir}${arg}" /n09data/hervas/SP_simu_fab/job_per_tile.job "$arg"
      
      tile_count=$((tile_count + 1))
      echo $tile_count
      sleep 15
    done

    cd /home/hervas/n25/SP_simu_fab/ || exit 1
  fi
done

