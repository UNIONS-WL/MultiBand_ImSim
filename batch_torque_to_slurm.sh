#!/bin/bash

# Check if there are any Torque scripts
if ! compgen -G "job_*" >/dev/null; then
    echo "No Torque scripts found."
    exit 1
fi

# Iterate over Torque scripts starting with "job_"
for TORQUE_SCRIPT in job_*; do
    # Define SLURM submission script template
    SLURM_TEMPLATE="#!/bin/bash
#SBATCH --job-name=job_updated_name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00

# Your commands here"

    # Extract job name from Torque script
    JOB_NAME=$(awk '/#PBS -N/{print $3}' $TORQUE_SCRIPT)

    # Replace PBS syntax with SLURM syntax
    SLURM_SCRIPT=$(sed -e 's/#PBS -N/#SBATCH --job-name/g' \
                      -e 's/#PBS -l nodes/#SBATCH --nodes/g' \
                      -e 's/#PBS -l ppn/#SBATCH --ntasks-per-node/g' \
                      -e 's/#PBS -l walltime/#SBATCH --time/g' \
                      -e 's/#PBS -l mem/#SBATCH --mem/g' \
                      -e 's/#PBS -M/#SBATCH --mail-user/g' \
                      -e 's/#PBS -m/#SBATCH --mail-type/g' \
                      -e 's/#PBS -o/#SBATCH --output/g' \
                      -e 's/#PBS -e/#SBATCH --error/g' \
                      $TORQUE_SCRIPT)

    # Substitute job name in SLURM template
    SLURM_SCRIPT="${SLURM_SCRIPT//job_updated_name/job_updated_$JOB_NAME}"

    # Output the transformed SLURM script to a new file
    echo "$SLURM_TEMPLATE" | sed "s|# Your commands here|$(echo "$SLURM_SCRIPT" | sed ':a;N;$!ba;s/\n/\\n/g')|g" > "job_updated_$TORQUE_SCRIPT"
done

