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
# block for a new tissue — see markers.csv in this run's root.
# `markers` (dotplot row order) comes from cell_type_map's key order, not a
# separate `panel` list — yaml::read_yaml() preserves YAML mapping key order,
# so this is never out of sync with cell_type_map by construction.
cell_type_map <- unlist(config$markers$cell_type_map)
markers       <- names(cell_type_map)
ct_levels     <- unlist(config$markers$ct_levels)

# ---- Panel-grid sizing for wrap_plots() figures faceted by lane ----------
# Used by step1_load/step3_doublet's per-lane plot lists. Reproduces the
# original validated run's exact sizing at n=8 lanes (nrow=2, width=20,
# height=8) while scaling correctly for any other lane count, instead of a
# fixed nrow=2/width=20/height=8 baked in for that one specific run's lane
# count — a run with more or fewer lanes was getting squeezed or wasting
# space (and for nrow specifically, mislaid out) under the old fixed values.
lane_grid_dims <- function(n, max_nrow = 2) {
  nrow <- min(max_nrow, n)
  ncol <- ceiling(n / nrow)
  list(nrow = nrow, width = 5 * ncol, height = 4 * nrow)
}
