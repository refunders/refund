#!/bin/bash
#SBATCH -J study1-cl2
#SBATCH -o ./logs/%x.%j.out
#SBATCH -e ./logs/%x.%j.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/refund
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24000mb
#SBATCH --time=18:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE

# Load modules
module load slurm_setup
module load r/4.3.3-gcc13-mkl

# Run CL2 extension (150 reps, 8 parallel workers)
Rscript ci-benchmark/sim-study-nongaussian-cl2-extension.R full
