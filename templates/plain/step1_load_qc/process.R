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
  here("step1_load_qc", "obj0_load_qc", glue("{lane_names}.qs2")),
  lane_names
)

dir.create(here("step1_load_qc", "obj0_load_qc"), recursive = TRUE, showWarnings = FALSE)

# ---- Per-lane loading + QC metrics (parallel, fork) --------------------------
# NOTE: the "plain" profile has no CRISPR/guide-capture assay at all, so this
# step merges the original pipeline's step1_load (h5 loading, Gene Expression
# only — no CRISPR assay) with step2_assign's pct.mt/pct.rb computation into a
# single step. There is no guide-assignment stage to separate them around.
cat("Loading", length(lane_names), "lanes...\n")

lane_results <- future_lapply(lane_names, function(lane) {

  raw <- Read10X_h5(h5_paths[[lane]])
  # Read10X_h5 returns a named list only for multi-modal h5 files (e.g. GEX +
  # Antibody Capture); a single-modality (pure Gene Expression) h5 returns a
  # bare dgCMatrix directly. Handle both so this template works either way.
  counts <- if (is.list(raw)) raw[["Gene Expression"]] else raw
  obj <- CreateSeuratObject(
    counts = counts,
    project = lane, min.cells = 3, min.features = 200
  )
  obj[["pct.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["pct.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

  qs_save(obj, file = obj_paths[[lane]])
  cat("[", lane, "]", ncol(obj), "cells\n")

  list(n_cells = ncol(obj))

}, future.seed = TRUE)

names(lane_results) <- lane_names
plan(sequential)
cat("\nCells per lane:\n")
walk(lane_names, function(l) cat(" ", l, ":", lane_results[[l]]$n_cells, "\n"))

cat("\n=== STEP 1 COMPLETE ===\n")
cat("Run step1_load_qc/plot.R to generate QC violin figures.\n")
