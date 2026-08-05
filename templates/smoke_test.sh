#!/bin/bash
#SBATCH --job-name=smoke_test
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --partition=__PARTITION__
#SBATCH --output=__LOG_DIR__/slurm_%j.log
#SBATCH --error=__LOG_DIR__/slurm_%j.err
#SBATCH --mail-user=__MAIL_USER__
#SBATCH --mail-type=BEGIN,END,FAIL

source ~/.bashrc
conda activate __CONDA_ENV__
cd __PROJECT_ROOT__

echo "[$(date)] SLURM smoke test starting"
echo "Partition: __PARTITION__"
echo "Project root (pwd): $(pwd)"
echo "Conda env: $(conda info --envs | grep '\*' || echo 'unknown')"
Rscript -e 'cat("R OK:", R.version.string, "\n"); suppressPackageStartupMessages(library(Seurat)); cat("Seurat OK:", as.character(packageVersion("Seurat")), "\n"); suppressPackageStartupMessages(library(here)); invisible(yaml::read_yaml("config.yaml")); cat("config.yaml readable OK\n")'
echo "[$(date)] SLURM smoke test done — if this line and the R checks above all printed with no errors, and you got the expected BEGIN/END emails, SLURM settings are good"
