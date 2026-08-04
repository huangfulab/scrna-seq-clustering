#!/bin/bash
#SBATCH --job-name=step2_doublet
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

echo "[$(date)] Starting step2_doublet/process.R"
Rscript step2_doublet/scripts/process.R
echo "[$(date)] Starting step2_doublet/plot.R"
Rscript step2_doublet/scripts/plot.R
echo "[$(date)] Done — review fig1_doublet_summary.png / tbl1_doublet_summary.csv, then proceed to step3"
