#!/bin/bash
#SBATCH --job-name=step8_seed
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

echo "[$(date)] Starting step8_seed/process.R"
Rscript step8_seed/scripts/process.R
echo "[$(date)] Starting step8_seed/plot.R"
Rscript step8_seed/scripts/plot.R
echo "[$(date)] Done — inspect figs/, choose config.final_clustering.winning_seed, then cp to step9_cluster_final/obj8_cluster_final.qs2"
