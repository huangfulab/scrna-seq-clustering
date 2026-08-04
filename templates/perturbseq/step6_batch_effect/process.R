# conda run -n <conda.env_name> Rscript step6_batch_effect/scripts/process.R
source(here::here("init.R"))

stopifnot(!is.null(config$clustering$n_dims))
RESOLUTIONS <- unlist(config$clustering$resolutions)

# Load ----------
stopifnot(file.exists(here("step5_PCA", "obj4_PCA.qs2")))
cat("Loading obj4_PCA...\n")
obj <- qs_read(here("step5_PCA", "obj4_PCA.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj), "cells\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")

# FindNeighbors ----------
cat(glue("\n[FindNeighbors] reduction=pca, dims=1:{config$clustering$n_dims}, k.param={config$clustering$k_param}\n"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:config$clustering$n_dims,
                     k.param = config$clustering$k_param, verbose = FALSE)

# RunUMAP ----------
cat(glue("\n[RunUMAP] reduction=pca, dims=1:{config$clustering$n_dims}, seed.use={config$clustering$umap_seed}\n"))
obj <- RunUMAP(obj, reduction = "pca", dims = 1:config$clustering$n_dims,
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
obj5_path <- here("step6_batch_effect", "obj5_UMAP.qs2")
qs_save(obj, file = obj5_path, nthreads = n_cores)
cat(glue("Saved: step6_batch_effect/obj5_UMAP.qs2 ({round(file.size(obj5_path)/1e6, 1)} MB)\n"))

cat("\n=== STEP 6 COMPLETE ===\n")
cat("Run step6_batch_effect/scripts/plot.R to generate UMAP figures.\n")
