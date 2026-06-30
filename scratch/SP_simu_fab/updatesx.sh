echo "You cleaned up outputs right?"
for dir in SP_????/; do
  echo $dir	
  if [ -d "$dir" ]; then
	  cd "$dir"
	  pwd
	  python /home/hervas/n25/SP_simu_fab/update_tiles_Sx.py
	  sleep 1
	  cd /home/hervas/n25/SP_simu_fab
  fi
done
