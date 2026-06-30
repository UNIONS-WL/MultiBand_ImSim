#!/bin/bash

# Check if there are any Torque scripts
if [ $# -eq 0 ]; then
    echo "Usage: $0 <torque_script>"
    exit 1
fi

TORQUE_SCRIPT=$1
echo 'modifying torque script' $TORQUE_SCRIPT
# Iterate over Torque scripts starting with "job_"

# Define SLURM submission script template
SLURM_TEMPLATE='#!/bin/bash
# Your commands here'

# Extract job name from Torque script
JOB_NAME=$(awk '/#PBS -N/{print $3}' $TORQUE_SCRIPT)

# Extract nodes and ppn from Torque script, if specified
NODES=$(awk -F= '/#PBS -l nodes=/{print $2}' $TORQUE_SCRIPT | cut -d":" -f1)
PPN=$(awk -F= '/ppn=/{print $3}' $TORQUE_SCRIPT | awk '{gsub(/[^0-9]/,"",$0); print}')
echo 'PPN in torque script=' $PPN
echo 'nodes in torque script=' $NODES
# Replace PBS syntax with SLURM syntax
    SLURM_SCRIPT=$(sed -e 's/#PBS -N/#SBATCH --job-name /g' \
                    -e "s/#PBS -l nodes=.*:ppn=.*/#SBATCH --nodes=$NODES \n#SBATCH --cpus-per-task=$PPN \n#SBATCH --ntask-per-node=1\n#SBATCH --ntasks-per-core=1/g" \
                    -e "s/#PBS -l nodes=$NODES/#SBATCH --nodes=$NODES/g" \
                    -e "s/#PBS -l ppn=$PPN/#SBATCH --cpus-per-task=$PPN/g" \
                    -e 's/#PBS -l walltime/#SBATCH --time/g' \
                    -e 's/#PBS -l mem/#SBATCH --mem/g' \
                    -e 's/#PBS -M/#SBATCH --mail-user/g' \
                    -e 's/#PBS -m ea/#SBATCH --mail-type END/g' \
                    -e 's/#PBS -o/#SBATCH --output/g' \
                    -e 's/#PBS -e/#SBATCH --error/g' \
                    -e '/ module/ i\export HDF5_USE_FILE_LOCKING="FALSE"\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK' \
                    -e '/#PBS -j oe/d' \
            -e '\|^#!/bin/bash$|d' \
            $TORQUE_SCRIPT)
# Substitute job name in SLURM template
SLURM_SCRIPT="${SLURM_SCRIPT//job_updated_name/job_updated_$JOB_NAME}"

# Output the transformed SLURM script to a new file
echo "$SLURM_TEMPLATE" | sed "s|# Your commands here|$(echo "$SLURM_SCRIPT" | sed ':a;N;$!ba;s/\n/\\n/g')|g" > "slurm_$TORQUE_SCRIPT"


