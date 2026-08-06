source(here::here("init.R"))

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))
plan(multicore, workers = n_workers)

# ---- User-set parameter (config.qc.umi_threshold) ----------
stopifnot(!is.null(config$qc$umi_threshold))

dir.create(here("step2_assign", "obj1_assign"), recursive = TRUE, showWarnings = FALSE)

# ---- Functions --------------------------------------------------------------

# Returns per-cell guide names (semicolon-sep) at given UMI threshold.
# Cells with no detected guides get NA.
extract_guides <- function(crispr_mat, threshold) {
  thresh_mat <- crispr_mat >= threshold
  nz         <- summary(thresh_mat)
  nz         <- nz[nz$x, ]                          # keep TRUE entries only
  if (nrow(nz) == 0) {
    return(setNames(rep(NA_character_, ncol(crispr_mat)), colnames(crispr_mat)))
  }
  per_cell <- tapply(
    rownames(crispr_mat)[nz$i], colnames(crispr_mat)[nz$j],
    function(g) paste(sort(g), collapse = ";")
  )
  result <- setNames(rep(NA_character_, ncol(crispr_mat)), colnames(crispr_mat))
  result[names(per_cell)] <- as.character(per_cell)
  result
}

# ---- Per-lane processing ---------------------------------------
cat("Processing", length(lane_names), "lanes — UMI threshold =", config$qc$umi_threshold, "\n")

lane_results <- future_lapply(lane_names, function(lane) {
  cat("Lane:", lane, "\n")
  obj <- qs_read(here("step1_load", "obj0_load", glue("{lane}.qs2")))
  crispr_mat <- GetAssayData(obj, assay = "CRISPR", layer = "counts")

  guides_vec <- extract_guides(crispr_mat, config$qc$umi_threshold)
  target_genes_vec <- map_chr(guides_vec, function(g) {
    if (is.na(g)) return(NA_character_)
    gnames <- strsplit(g, ";")[[1]]
    genes  <- map_chr(gnames, function(gn) {
      if (gn %in% names(oeg_gene_map)) oeg_gene_map[[gn]] else "NegativeControl"
    })
    paste(unique(genes), collapse = ";")
  })

  obj[["guide_origin"]] <- obj@meta.data[[glue("guide_origin_UMI{config$qc$umi_threshold}")]]
  obj[["MOI"]]          <- obj@meta.data[[glue("n_guides_UMI{config$qc$umi_threshold}")]]
  obj[["guides"]]       <- unname(guides_vec[colnames(obj)])
  obj[["target_genes"]] <- unname(target_genes_vec[colnames(obj)])
  obj[["pct.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["pct.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

  # Drop all threshold sweep columns
  drop_cols <- intersect(
    colnames(obj@meta.data),
    c(glue("guide_origin_UMI{1:10}"), glue("n_guides_UMI{1:10}"))
  )
  obj@meta.data <- obj@meta.data[
    , !colnames(obj@meta.data) %in% drop_cols, drop = FALSE
  ]

  qs_save(obj, file = here("step2_assign", "obj1_assign", glue("{lane}.qs2")))
  cat(" saved:", ncol(obj), "cells\n")

  list(n_cells = ncol(obj))

}, future.seed = TRUE)

names(lane_results) <- lane_names
plan(sequential)
cat("\nCells per lane:\n")
walk(lane_names, function(l) cat(" ", l, ":", lane_results[[l]]$n_cells, "\n"))

cat("\n=== STEP 2 COMPLETE ===\n")
cat("Run step2_assign/plot.R to generate QC violin plots.\n")
cat("Proceed to step3_doublet/process.R\n")
