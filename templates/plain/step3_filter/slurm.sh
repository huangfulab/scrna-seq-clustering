#!/bin/bash
#SBATCH --job-name=step3_filter
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

echo "[$(date)] Starting step3_filter/process.R"
Rscript step3_filter/scripts/process.R
echo "[$(date)] Starting step3_filter/plot.R"
Rscript step3_filter/scripts/plot.R
echo "[$(date)] Done — review fig1_upset.png / fig2_qc_filter.png / tbl1_filter_summary.csv, then proceed to step4"
