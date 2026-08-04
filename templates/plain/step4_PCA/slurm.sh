#!/bin/bash
#SBATCH --job-name=step4_PCA
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --mem=__MEM__
#SBATCH --time=__TIME__
#SBATCH --partition=__PARTITION__
#SBATCH --output=__LOG_DIR__/slurm_%j.log
#SBATCH --error=__LOG_DIR__/slurm_%j.err
#SBATCH --mail-user=__MAIL_USER__
#SBATCH --mail-type=BEGIN,END,FAIL

source ~/.bashrc
conda activate __CONDA_ENV__
cd __PROJECT_ROOT__

echo "[$(date)] Starting step4_PCA/process.R (no Harmony)"
Rscript step4_PCA/scripts/process.R
echo "[$(date)] Starting step4_PCA/plot.R"
Rscript step4_PCA/scripts/plot.R
echo "[$(date)] Done — inspect figs/fig1_elbow.png, confirm config.clustering.n_dims, then proceed to step5"
