#!/bin/bash
#SBATCH --job-name=step2_assign
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

echo "[$(date)] Starting step2_assign/process.R (config.qc.umi_threshold set in config.yaml)"
Rscript step2_assign/scripts/process.R
echo "[$(date)] Starting step2_assign/plot.R"
Rscript step2_assign/scripts/plot.R
echo "[$(date)] Done — review fig1 cutoffs, then proceed to step3"
