# conda run -n <conda.env_name> Rscript step8_res/scripts/process.R
source(here::here("init.R"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))

stopifnot(!is.null(config$clustering$n_dims))
RESOLUTIONS <- unlist(config$clustering$resolutions)
vars_to_regress <- unlist(config$pca$vars_to_regress)

# Single-lane shortcut: with only one lane there is no batch effect to check
# (nothing to compare one lane against), so step6_batch_effect/step7_rm_badcl
# are skippable as a pair for single-lane runs. Fall back to step5_PCA's
# object directly in that case — everything below recomputes HVG/scale/PCA/
# neighbors/UMAP/clusters from scratch regardless of which object we start
# from, so feeding it step5's un-cleaned object instead of step7's
# cluster-cleaned one is safe.
if (length(config$lanes) == 1) {
  cat("Single-lane run (config$lanes has 1 entry) — loading obj4_PCA directly,\n")
  cat("skipping step6_batch_effect/step7_rm_badcl...\n")
  obj <- qs_read(here("step5_PCA", "obj4_PCA.qs2"), nthreads = n_cores)
} else {
  cat("Loading obj6_rm_badcl...\n")
  obj <- qs_read(here("step7_rm_badcl", "obj6_rm_badcl.qs2"), nthreads = n_cores)
}
cat("Loaded:", ncol(obj), "cells\n")
cat("Phase distribution:\n"); print(table(obj$Phase))

cat(glue("\n[1] FindVariableFeatures (nfeatures = {config$pca$hvg_nfeatures})...\n"))
obj <- FindVariableFeatures(obj, nfeatures = config$pca$hvg_nfeatures, verbose = FALSE)
cat("Top 10 HVGs:", paste(head(VariableFeatures(obj), 10), collapse = ", "), "\n")

cat(glue("\n[2] ScaleData (vars.to.regress: {paste(vars_to_regress, collapse = ', ')})...\n"))
plan("multicore", workers = n_cores)
obj <- ScaleData(obj, features = VariableFeatures(obj),
                 vars.to.regress = vars_to_regress, verbose = TRUE)
plan("sequential")

cat(glue("\n[3] RunPCA (npcs = {config$pca$npcs}, seed = {config$pca$seed})...\n"))
set.seed(config$pca$seed)
obj <- RunPCA(obj, npcs = config$pca$npcs, verbose = FALSE)
cat("Top genes PC1-2:\n")
print(obj[["pca"]], dims = 1:2, nfeatures = 5)

cat(glue("\n[4] FindNeighbors (dims=1:{config$clustering$n_dims}, k={config$clustering$k_param})...\n"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:config$clustering$n_dims,
                     k.param = config$clustering$k_param, verbose = FALSE)

cat(glue("\n[5] RunUMAP (dims=1:{config$clustering$n_dims}, seed={config$clustering$umap_seed})...\n"))
obj <- RunUMAP(obj, reduction = "pca", dims = 1:config$clustering$n_dims,
               seed.use = config$clustering$umap_seed, umap.method = "uwot", verbose = FALSE)

cat("\n[6] FindClusters resolutions:", paste(RESOLUTIONS, collapse = ", "), "\n")
plan("multicore", workers = n_cores)
obj <- FindClusters(obj, resolution = RESOLUTIONS, algorithm = config$clustering$algorithm, verbose = FALSE)
plan("sequential")
walk(RESOLUTIONS, function(res) {
  col <- glue("RNA_snn_res.{res}")
  cat(glue("  res={res}: {nlevels(obj[[col, drop=TRUE]])} clusters\n"))
})

cat("\nSaving obj7_cluster2...\n")
obj7_path <- here("step8_res", "obj7_cluster2.qs2")
qs_save(obj, file = obj7_path, nthreads = n_cores)
cat(glue("Saved: step8_res/obj7_cluster2.qs2 ({round(file.size(obj7_path)/1e6, 1)} MB)\n"))

cat("\n=== STEP 8 PROCESS.R COMPLETE ===\n")
cat("Run step8_res/scripts/plot.R to generate figures and tables.\n")
