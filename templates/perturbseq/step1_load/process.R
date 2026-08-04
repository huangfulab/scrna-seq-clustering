source(here::here("init.R"))

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))
plan(multicore, workers = n_workers)

# ---- Paths ------------------------------------------------------------------
h5_paths  <- setNames(
  here(config$paths$h5_dir,
       map_chr(lane_names, function(lane) glue(config$paths$h5_pattern))),
  lane_names
)
obj_paths <- setNames(
  here("step1_load", "obj0_load", glue("{lane_names}.qs2")),
  lane_names
)

dir.create(here("step1_load", "obj0_load"), recursive = TRUE, showWarnings = FALSE)

# ---- Functions --------------------------------------------------------------
classify_cells <- function(thresh_mat, nc_set, oeg_set) {
  has_nc  <- colSums(
    thresh_mat[rownames(thresh_mat) %in% nc_set,  , drop = FALSE]
  ) > 0
  has_oeg <- colSums(
    thresh_mat[rownames(thresh_mat) %in% oeg_set, , drop = FALSE]
  ) > 0
  case_when(
    has_nc  & has_oeg  ~ "Mixed",
    has_nc  & !has_oeg ~ "NC",
    !has_nc & has_oeg  ~ "OEG",
    TRUE               ~ "Unassigned"
  )
}

# ---- Per-lane loading (parallel, fork) --------------------------------------
cat("Loading", length(lane_names), "lanes...\n")

lane_results <- future_lapply(lane_names, function(lane) {

  raw <- Read10X_h5(h5_paths[[lane]])
  obj <- CreateSeuratObject(
    counts = raw[["Gene Expression"]],
    project = lane, min.cells = 3, min.features = 200
  )
  crispr_counts <- raw[["CRISPR Guide Capture"]][, colnames(obj), drop = FALSE]
  obj[["CRISPR"]] <- CreateAssayObject(counts = crispr_counts)
  crispr_mat <- GetAssayData(obj, assay = "CRISPR", layer = "counts")

  # Precompute threshold matrices once; reuse for metadata + recovery
  thresh_mats <- map(1:10, ~ crispr_mat >= .x)

  if (!any(rownames(thresh_mats[[1]]) %in% nc_set))
    warning("No NC guides matched in thresh_mats")
  if (!any(rownames(thresh_mats[[1]]) %in% oeg_set))
    warning("No OEG guides matched in thresh_mats")

  origin_cols <- set_names(map(thresh_mats, classify_cells, nc_set, oeg_set),
                           glue("guide_origin_UMI{1:10}"))
  nguide_cols <- set_names(map(thresh_mats, ~ as.integer(colSums(.x))),
                           glue("n_guides_UMI{1:10}"))
  meta_df           <- as.data.frame(c(origin_cols, nguide_cols))
  rownames(meta_df) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = meta_df)
  qs_save(obj, file = obj_paths[[lane]])
  cat("[", lane, "]", ncol(obj), "cells\n")

  list(n_cells = ncol(obj))

}, future.seed = TRUE)

names(lane_results) <- lane_names
plan(sequential)
cat("\nCells per lane:\n")
walk(lane_names, function(l) cat(" ", l, ":", lane_results[[l]]$n_cells, "\n"))

cat("\n=== STEP 1 COMPLETE ===\n")
cat("Run step1_load/plot.R to generate QC figures.\n")
cat("Set config$qc$umi_threshold in config.yaml before running step2.\n")
