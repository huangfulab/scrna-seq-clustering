# conda run -n <conda.env_name> Rscript step9_seed/scripts/process.R
source(here::here("init.R"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))

stopifnot(!is.null(config$final_clustering$final_resolution))
SEEDS      <- unlist(config$final_clustering$seeds)
RESOLUTION <- config$final_clustering$final_resolution

# No-bad-cluster shortcut: if step6_batch_effect's checkpoint found nothing
# worth removing (config$bad_cluster_removal$no_bad_cluster == TRUE),
# step7_rm_badcl/step8_res are skipped entirely and this step reads straight
# from step6_batch_effect's object instead of step8_res's — recomputing HVG/
# scale/PCA/neighbors/UMAP/clusters on the exact same cells step6_batch_effect
# already clustered would only reproduce that same result.
if (isTRUE(config$bad_cluster_removal$no_bad_cluster)) {
  cat("No bad cluster was found (config$bad_cluster_removal$no_bad_cluster == TRUE) —\n")
  cat("skipping step7_rm_badcl/step8_res, loading obj5_UMAP directly...\n")
  obj <- qs_read(here("step6_batch_effect", "obj5_UMAP.qs2"), nthreads = n_cores)
} else {
  cat("Loading obj7_cluster2...\n")
  obj <- qs_read(here("step8_res", "obj7_cluster2.qs2"), nthreads = n_cores)
}
cat("Loaded:", ncol(obj), "cells\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("Graphs:", paste(Graphs(obj), collapse = ", "), "\n")

# Strip prior UMAP/neighbor/cluster products ----------
cat("\nStripping prior umap reduction, RNA_nn/RNA_snn graphs, RNA_snn_res.*/seurat_clusters columns...\n")
obj[["umap"]]   <- NULL
obj[["RNA_nn"]] <- NULL
obj[["RNA_snn"]] <- NULL
res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
obj@meta.data[, c(res_cols, "seurat_clusters")] <- NULL
cat("Kept reductions:", paste(Reductions(obj), collapse = ", "), "\n")

# Seed sweep: FindNeighbors + RunUMAP + FindClusters per seed, save inside each worker ----------
obj_dir <- here("step9_seed", "obj")
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

cat(glue("\nSweeping {length(SEEDS)} seeds (dims=1:{config$clustering$n_dims}, k={config$clustering$k_param}, res={RESOLUTION})...\n"))
plan("multicore", workers = n_cores)
results <- future_lapply(SEEDS, function(s) {
  tmp <- FindNeighbors(obj, reduction = "pca", dims = 1:config$clustering$n_dims,
                       k.param = config$clustering$k_param, verbose = FALSE)
  tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:config$clustering$n_dims,
                 seed.use = s, umap.method = "uwot", verbose = FALSE)
  tmp <- FindClusters(tmp, resolution = RESOLUTION,
                      random.seed = s, algorithm = config$clustering$algorithm, verbose = FALSE)
  Idents(tmp) <- tmp$seurat_clusters

  out_path <- file.path(obj_dir, glue("seed{s}.qs2"))
  qs_save(tmp, file = out_path, nthreads = 1)

  list(seed = s, n_clusters = nlevels(Idents(tmp)), path = out_path)
}, future.seed = FALSE)
plan("sequential")

# Report ----------
walk(results, function(res) {
  cat(glue("  seed={res$seed}: {res$n_clusters} clusters -> {res$path} ({round(file.size(res$path)/1e9, 2)} GB)\n"))
})

cat("\n=== STEP 9 COMPLETE ===\n")
cat("Run step9_seed/scripts/plot.R to generate figures.\n")
