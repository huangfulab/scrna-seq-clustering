#!/bin/bash
#SBATCH --job-name=step7_rm_badcl
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

echo "[$(date)] Starting step7_rm_badcl/process.R"
Rscript step7_rm_badcl/scripts/process.R
echo "[$(date)] Starting step7_rm_badcl/plot.R"
Rscript step7_rm_badcl/scripts/plot.R
echo "[$(date)] Done — inspect figs/, then submit step8_res/scripts/slurm.sh"
