#!/bin/bash
#SBATCH --job-name=step7_res
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

echo "[$(date)] Starting step7_res/process.R"
Rscript step7_res/scripts/process.R
echo "[$(date)] Starting step7_res/plot.R"
Rscript step7_res/scripts/plot.R
echo "[$(date)] Done — inspect figs/, choose config.final_clustering.final_resolution, then proceed to step8_seed"
