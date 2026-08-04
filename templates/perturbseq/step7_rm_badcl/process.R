# conda run -n <conda.env_name> Rscript step7_rm_badcl/scripts/process.R
source(here::here("init.R"))

# Filled in config.yaml after inspecting step6_batch_effect/figs — the lowest resolution
# that clearly separates a low-quality/outlier cluster, and its cluster ID.
stopifnot(!is.null(config$bad_cluster_removal$resolution), !is.null(config$bad_cluster_removal$cluster_rm))
RESOLUTION <- config$bad_cluster_removal$resolution
CLUSTER_RM <- as.character(config$bad_cluster_removal$cluster_rm)

cat("Loading obj5_UMAP...\n")
obj5 <- qs_read(here("step6_batch_effect", "obj5_UMAP.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj5), "cells\n")

res_col <- glue("RNA_snn_res.{RESOLUTION}")
keep    <- obj5@meta.data[[res_col]] != CLUSTER_RM
cat("Removing cluster", CLUSTER_RM, "at res", RESOLUTION, ":", sum(!keep), "cells removed\n")

obj6 <- obj5[, keep]
cat("Remaining:", ncol(obj6), "cells\n")

qs_save(obj6, here("step7_rm_badcl", "obj6_rm_badcl.qs2"), nthreads = n_cores)
cat("Saved: step7_rm_badcl/obj6_rm_badcl.qs2\n")
cat("\n=== STEP 7_RM_BADCL PROCESS.R COMPLETE ===\n")
