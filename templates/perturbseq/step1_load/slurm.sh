#!/bin/bash
#SBATCH --job-name=step1_load
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

echo "[$(date)] Starting step1_load/process.R"
Rscript step1_load/scripts/process.R
echo "[$(date)] Starting step1_load/plot.R"
Rscript step1_load/scripts/plot.R
echo "[$(date)] Done — review fig1-3, guide_origin_summary.csv, set config.qc.umi_threshold in config.yaml"
