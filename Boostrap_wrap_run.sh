#!/bin/bash
#SBATCH -J BOOTstrap
#SBATCH -o BOOTstrap_.%j.out
#SBATCH -e BOOTstrap_.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00

set -euo pipefail

# Work in the folder that contains your script and results/
cd /Genomics/ayroleslab2/benja/NMR_analysis/My_First_Paper/

# --- choose ONE environment setup ---

# (A) Conda/mamba env with R + ggplot2 + ggrepel + qqman
source ~/miniforge3/etc/profile.d/conda.sh || source ~/mambaforge/etc/profile.d/conda.sh
conda activate rplots  # <-- swap to your actual R env name

# (or) (B) Princeton modules, if you prefer
# module purge
# module load R/4.3.2

# Avoid oversubscription
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

# Run
Rscript --vanilla Repeat_Amanda_bootstrap.R
