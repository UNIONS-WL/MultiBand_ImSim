#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --output=/scratch/cdaley/imsim/logs/run2/%x-%j.out
# Run 2 of the UNIONS twin on Fir (config_grid2.yaml), one stage per batch job:
#
#   sbatch -c 48 --mem=192G -t 6:00:00 -J im2-draw0 unions/run2.sbatch draw 0 4
#   sbatch -c 64 --mem=320G -t 12:00:00 -J im2-coadd unions/run2.sbatch coadd
#   sbatch -c 32 --mem=128G -t 3:00:00 -J im2-check unions/run2.sbatch check 253.246
#
# draw I N: draws every N-th exposure from I (shard 0 writes the truth first;
#   the others wait for it), after waiting for the staged campaign to be complete.
# coadd: coadds every tile in every branch and writes the PSF dictionary.
# check TILE: known-answer checks on one tile.
# Submit from the repository root. A failure is recorded in
# /scratch/cdaley/imsim/export-grid2/FAILED.
set -euo pipefail
trap 'echo "$(date) ${SLURM_JOB_NAME} ${SLURM_JOB_ID} failed (stage $*)" >> /scratch/cdaley/imsim/export-grid2/FAILED' ERR
mkdir -p /scratch/cdaley/imsim/export-grid2
module load apptainer
CFG=unions/config_grid2.yaml
SIF=/project/rrg-kyi/cdaley/unions/containers/sp_validation-develop.sif
CAMP=$(awk '/^campaign_dir:/{print $2}' $CFG)
ROOT=$(awk '/^output_root:/{print $2}' $CFG)
K=$SLURM_CPUS_PER_TASK
py() { apptainer exec --bind /scratch,/project $SIF env OMP_NUM_THREADS=1 python "$@"; }
echo "start $(date) $*"

case $1 in
  draw)
    for e in $(cat $CAMP/exposures.txt); do
      for k in "" .weight .flag; do
        until [ -s $CAMP/exposures/${e}p$k.fits.fz ]; do sleep 60; done
      done
    done
    if [ "$2" != 0 ]; then
      until [ -s $ROOT/grid_2_shared/truth.fits ]; do sleep 30; done
    fi
    py -m unions.twin $CFG -j $K --shard $2 $3 ;;
  coadd)
    py -m unions.twin $CFG -j $K ;;
  check)
    until [ -s $CAMP/tiles/CFIS.$2.r.weight.fits.fz ]; do sleep 60; done
    py -m unions.checks $CFG --tile $2 -j $K ;;
esac
echo "done $(date)"
