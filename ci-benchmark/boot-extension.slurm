#!/bin/bash
# ===========================================================================
# Bootstrap CI Extension — SLURM Array Job for LRZ CoolMUC-4
# ===========================================================================
#
# One array task = one DGP cell; the R script loops over its reps sequentially.
#
# SUBMISSION:
#   # Study 1: 12 DGP cells  (array indices 1..12)
#   sbatch --clusters=serial --export=STUDY=study1 \
#          --array=1-12 boot-extension.slurm
#
#   # Study 2: 18 DGP cells  (array indices 1..18)
#   sbatch --clusters=serial --export=STUDY=study2 \
#          --array=1-18 boot-extension.slurm
#
# SIZING RATIONALE:
#   A single pffr_coefboot() call at B=499 on a Study 1 cell (n=400, small
#   grid) takes ~10–20 min on 8 cores. Study 2 cells have 3 nygrid levels,
#   and the finest (nygrid=120, n=80) takes ~15 min per call. Total wall
#   time per array task (all reps): 150 reps × ~15 min = ~37h for study1
#   at n=400, easily fits in 48h. Study 2 (50 reps × 3 grids) is similar.
#   8 cores balance bootstrap speed vs. serial queue core limits.
#   --mem: 8 cores × ~3GB/core = ~24GB; request 32GB for headroom.
#
# STUDY env var:
#   Set via --export=STUDY=study1 or --export=STUDY=study2 at sbatch time.
#   Falls back to "study1" if unset.
# ===========================================================================

#SBATCH -J boot_ci
#SBATCH -o ./logs/%x.%A_%a.out
#SBATCH -e ./logs/%x.%A_%a.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/refund
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000mb
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE

# Load modules
module load slurm_setup
module load r/4.3.3-gcc13-mkl

# User library (must match the active R version's path)
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# Parallelism: pass allocated CPUs to the R script via NCPUS
export NCPUS="${SLURM_CPUS_PER_TASK:-1}"

# Study to run: set via --export=STUDY=study1 or --export=STUDY=study2
STUDY="${STUDY:-study1}"

echo "=== boot-extension.slurm ==="
echo "  STUDY              : ${STUDY}"
echo "  SLURM_ARRAY_JOB_ID : ${SLURM_ARRAY_JOB_ID}"
echo "  SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "  NCPUS              : ${NCPUS}"
echo "  HOST               : $(hostname)"
echo "  DATE               : $(date)"
echo "============================="

# SLURM_ARRAY_TASK_ID is the dgp_id; the R script maps it to a DGP cell
Rscript ci-benchmark/sim-study-boot-extension.R \
  "${STUDY}" full "${SLURM_ARRAY_TASK_ID}"

echo "Done at $(date)"
