#!/bin/bash

# Check if there are any Torque scripts
if ! compgen -G "job_*" >/dev/null; then
    echo "No Torque scripts found."
    exit 1
fi

# Iterate over Torque scripts starting with "job_"
for TORQUE_SCRIPT in job_*; do
    # Define SLURM submission script template
    SLURM_TEMPLATE='#!/bin/bash
export HDF5_USE_FILE_LOCKING="FALSE"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Your commands here'

    # Extract job name from Torque script
    JOB_NAME=$(awk '/#PBS -N/{print $3}' $TORQUE_SCRIPT)
    
    # Extract nodes and ppn from Torque script, if specified
    NODES=$(awk -F= '/#PBS -l nodes=/{print $2}' $TORQUE_SCRIPT | cut -d":" -f1)
    PPN=$(awk -F= '/ppn=/{print $3}' $TORQUE_SCRIPT | awk '{gsub(/[^0-9]/,"",$0); print}')
    echo $PPN
    echo $NODES
    # Replace PBS syntax with SLURM syntax
        SLURM_SCRIPT=$(sed -e 's/#PBS -N/#SBATCH --job-name /g' \
                      -e "s/#PBS -l nodes=.*:ppn=.*/#SBATCH --nodes=$NODES \n#SBATCH --cpus-per-task=$PPN \n#SBATCH --ntask-per-node=1/g" \
                      -e "s/#PBS -l nodes=$NODES/#SBATCH --nodes=$NODES/g" \
                      -e "s/#PBS -l ppn=$PPN/#SBATCH --cpus-per-task=$PPN/g" \
                      -e 's/#PBS -l walltime/#SBATCH --time/g' \
                      -e 's/#PBS -l mem/#SBATCH --mem/g' \
                      -e 's/#PBS -M/#SBATCH --mail-user/g' \
                      -e 's/#PBS -m ea/#SBATCH --mail-type END/g' \
                      -e 's/#PBS -o/#SBATCH --output/g' \
                      -e 's/#PBS -e/#SBATCH --error/g' \
                      -e "s/python /srun python /g" \
		      -e '/#PBS -j oe/d' \
		      -e '\|^#!/bin/bash$|d' \
		       $TORQUE_SCRIPT)
    # Substitute job name in SLURM template
    SLURM_SCRIPT="${SLURM_SCRIPT//job_updated_name/job_updated_$JOB_NAME}"

    # Output the transformed SLURM script to a new file
    echo "$SLURM_TEMPLATE" | sed "s|# Your commands here|$(echo "$SLURM_SCRIPT" | sed ':a;N;$!ba;s/\n/\\n/g')|g" > "job_updated_$TORQUE_SCRIPT"
done

