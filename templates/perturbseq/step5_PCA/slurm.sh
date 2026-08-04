#!/bin/bash
#SBATCH --job-name=step5_PCA
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

echo "[$(date)] Starting step5_PCA/process.R (no Harmony)"
Rscript step5_PCA/scripts/process.R
echo "[$(date)] Starting step5_PCA/plot.R"
Rscript step5_PCA/scripts/plot.R
echo "[$(date)] Done — inspect figs/fig1_elbow.png, confirm config.clustering.n_dims, then submit step6_batch_effect/scripts/slurm.sh"
