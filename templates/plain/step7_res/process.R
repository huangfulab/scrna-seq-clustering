# conda run -n <conda.env_name> Rscript step7_res/scripts/process.R
source(here::here("init.R"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))

RESOLUTIONS <- unlist(config$clustering$resolutions)
N_DIMS      <- config$clustering$n_dims
K_PARAM     <- config$clustering$k_param

# Single-lane shortcut: with only one lane there is no cross-lane batch
# effect to check, so step5_batch_effect/step6_rm_badcl are skipped entirely
# and this step reads straight from step4_PCA's output instead of
# step6_rm_badcl's. Multi-lane runs read the bad-cluster-removed object.
if (length(config$lanes) == 1) {
  cat("Single-lane run — loading obj4_PCA directly (step5/step6 skipped)...\n")
  obj <- qs_read(here("step4_PCA", "obj4_PCA.qs2"), nthreads = n_cores)
} else {
  cat("Loading obj6_rm_badcl...\n")
  obj <- qs_read(here("step6_rm_badcl", "obj6_rm_badcl.qs2"), nthreads = n_cores)
}
cat("Loaded:", ncol(obj), "cells\n")
cat("Phase distribution:\n"); print(table(obj$Phase))

cat("\n[1] FindVariableFeatures (nfeatures =", config$pca$hvg_nfeatures, ")...\n")
obj <- FindVariableFeatures(obj, nfeatures = config$pca$hvg_nfeatures, verbose = FALSE)
cat("Top 10 HVGs:", paste(head(VariableFeatures(obj), 10), collapse = ", "), "\n")

cat("\n[2] ScaleData (vars.to.regress:", paste(unlist(config$pca$vars_to_regress), collapse = ", "), ")...\n")
plan("multicore", workers = n_cores)
obj <- ScaleData(obj, features = VariableFeatures(obj),
                 vars.to.regress = unlist(config$pca$vars_to_regress), verbose = TRUE)
plan("sequential")

cat("\n[3] RunPCA (npcs =", config$pca$npcs, ", seed =", config$pca$seed, ")...\n")
set.seed(config$pca$seed)
obj <- RunPCA(obj, npcs = config$pca$npcs, verbose = FALSE)
cat("Top genes PC1-2:\n")
print(obj[["pca"]], dims = 1:2, nfeatures = 5)

cat(glue("\n[4] FindNeighbors (dims=1:{N_DIMS}, k={K_PARAM})...\n"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:N_DIMS,
                     k.param = K_PARAM, verbose = FALSE)

cat(glue("\n[5] RunUMAP (dims=1:{N_DIMS}, seed={config$clustering$umap_seed})...\n"))
obj <- RunUMAP(obj, reduction = "pca", dims = 1:N_DIMS,
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
obj7_path <- here("step7_res", "obj7_cluster2.qs2")
qs_save(obj, file = obj7_path, nthreads = n_cores)
cat(glue("Saved: step7_res/obj7_cluster2.qs2 ({round(file.size(obj7_path)/1e6, 1)} MB)\n"))

cat("\n=== STEP 7 PROCESS.R COMPLETE ===\n")
cat("Run step7_res/plot.R to generate figures and tables.\n")
