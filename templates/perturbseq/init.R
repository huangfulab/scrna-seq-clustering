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

stopifnot(
  file.exists(here(config$paths$nc_gRNA_csv)),
  file.exists(here(config$paths$oeg_gRNA_csv))
)
nc_set  <- read_csv(here(config$paths$nc_gRNA_csv),  show_col_types = FALSE)$name
oeg_set <- read_csv(here(config$paths$oeg_gRNA_csv), show_col_types = FALSE)$name

oeg_gene_map <- setNames(sub("-[0-9]+$", "", oeg_set), oeg_set)

n_nc_total  <- length(nc_set)
n_oeg_total <- length(oeg_set)
n_all_total <- n_nc_total + n_oeg_total

# ---- Marker gene panel (config.markers.*) ----------
# Shared across step6/step8/step9/step10/step11 plot.R dotplots
# (gene -> marker-category strip labels). Replace config.yaml's markers.*
# block for a new tissue — see the heavily-commented config.yaml.template.
markers       <- unlist(config$markers$panel)
cell_type_map <- unlist(config$markers$cell_type_map)
ct_levels     <- unlist(config$markers$ct_levels)
