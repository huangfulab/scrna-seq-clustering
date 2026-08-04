suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(here)
  library(RColorBrewer)
  library(future)
  library(future.apply)
  library(Matrix)
  library(glue)
  library(patchwork)
  library(yaml)
})

options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = Inf)

pal <- brewer.pal(9, "Set1")

theme_set(theme_bw(base_family = "Arial"))
theme_update(
  plot.title       = element_text(hjust = 0.5),
  plot.subtitle    = element_text(hjust = 0.5),
  plot.caption     = element_text(hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# ---- Config ----------
# Every step script's `source(here::here("init.R"))` gets `config` for free —
# do not re-read config.yaml in any process.R/plot.R.
config <- yaml::read_yaml(here("config.yaml"))

n_cores    <- config$n_cores
lane_names <- unlist(config$lanes)

# NOTE: the "plain" profile has no CRISPR/guide-capture assay at all, so
# unlike templates/perturbseq/init.R there is no nc_set/oeg_set/oeg_gene_map
# reading here — config.yaml's paths.nc_gRNA_csv/oeg_gRNA_csv keys are
# perturbseq-only and are never read by any plain-profile step script.

# ---- Marker gene panel (config.markers.*) ----------
# Shared across step5/step7/step8/step9/step10 plot.R dotplots
# (gene -> marker-category strip labels). Replace config.yaml's markers.*
# block for a new tissue — see the heavily-commented config.yaml.template.
markers       <- unlist(config$markers$panel)
cell_type_map <- unlist(config$markers$cell_type_map)
ct_levels     <- unlist(config$markers$ct_levels)
