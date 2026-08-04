#!/bin/bash
#SBATCH --job-name=step4_filter
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

echo "[$(date)] Starting step4_filter/process.R"
Rscript step4_filter/scripts/process.R
echo "[$(date)] Starting step4_filter/plot.R"
Rscript step4_filter/scripts/plot.R
echo "[$(date)] Done — review filter_summary.csv, then proceed to step5"
