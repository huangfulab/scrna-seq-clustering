# conda run -n <conda.env_name> Rscript step5_batch_effect/scripts/process.R
source(here::here("init.R"))

N_DIMS      <- config$clustering$n_dims
K_PARAM     <- config$clustering$k_param
RESOLUTIONS <- unlist(config$clustering$resolutions)

# Load ----------
stopifnot(file.exists(here("step4_PCA", "obj4_PCA.qs2")))
cat("Loading obj4_PCA...\n")
obj <- qs_read(here("step4_PCA", "obj4_PCA.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj), "cells\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")

# FindNeighbors ----------
cat(glue("\n[FindNeighbors] reduction=pca, dims=1:{N_DIMS}, k.param={K_PARAM}\n"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:N_DIMS,
                     k.param = K_PARAM, verbose = FALSE)

# RunUMAP ----------
cat(glue("\n[RunUMAP] reduction=pca, dims=1:{N_DIMS}, seed.use={config$clustering$umap_seed}\n"))
obj <- RunUMAP(obj, reduction = "pca", dims = 1:N_DIMS,
               seed.use = config$clustering$umap_seed, umap.method = "uwot", verbose = FALSE)

# FindClusters ----------
cat("\n[FindClusters] Louvain, resolutions:", paste(RESOLUTIONS, collapse = ", "), "\n")
plan("multicore", workers = n_cores)
obj <- FindClusters(obj, resolution = RESOLUTIONS, algorithm = config$clustering$algorithm, verbose = FALSE)
plan("sequential")
walk(RESOLUTIONS, function(res) {
  col <- glue("RNA_snn_res.{res}")
  cat(glue("  res={res}: {nlevels(obj[[col, drop=TRUE]])} clusters\n"))
})

# Save ----------
cat("\nSaving obj5_UMAP...\n")
obj5_path <- here("step5_batch_effect", "obj5_UMAP.qs2")
qs_save(obj, file = obj5_path, nthreads = n_cores)
cat(glue("Saved: step5_batch_effect/obj5_UMAP.qs2 ({round(file.size(obj5_path)/1e6, 1)} MB)\n"))

cat("\n=== STEP 5 COMPLETE ===\n")
cat("Run step5_batch_effect/plot.R to generate UMAP figures.\n")
