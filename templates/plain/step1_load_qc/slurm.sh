#!/bin/bash
#SBATCH --job-name=step1_load_qc
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

echo "[$(date)] Starting step1_load_qc/process.R"
Rscript step1_load_qc/scripts/process.R
echo "[$(date)] Starting step1_load_qc/plot.R"
Rscript step1_load_qc/scripts/plot.R
echo "[$(date)] Done — review fig1_qc_violin.png, set config.qc.{min_ncount,min_nfeat,max_pct_mt} in config.yaml"
