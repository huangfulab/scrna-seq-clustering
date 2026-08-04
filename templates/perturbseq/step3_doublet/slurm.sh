#!/bin/bash
#SBATCH --job-name=step3_doublet
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

echo "[$(date)] Starting step3_doublet/process.R"
Rscript step3_doublet/scripts/process.R
echo "[$(date)] Starting step3_doublet/plot.R"
Rscript step3_doublet/scripts/plot.R
echo "[$(date)] Done — review doublet_summary.csv, then proceed to step4"
