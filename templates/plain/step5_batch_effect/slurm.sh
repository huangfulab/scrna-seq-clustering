#!/bin/bash
#SBATCH --job-name=step5_batch_effect
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

echo "[$(date)] Starting step5_batch_effect/process.R"
Rscript step5_batch_effect/scripts/process.R
echo "[$(date)] Starting step5_batch_effect/plot.R"
Rscript step5_batch_effect/scripts/plot.R
echo "[$(date)] Done — inspect figs/, choose resolution + bad cluster, then fill step6_rm_badcl placeholders"
