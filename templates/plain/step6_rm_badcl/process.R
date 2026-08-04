# conda run -n <conda.env_name> Rscript step6_rm_badcl/scripts/process.R
source(here::here("init.R"))

# Filled in by the user after inspecting step5_batch_effect/figs — see
# config.bad_cluster_removal.{resolution,cluster_rm}. This step is only ever
# scaffolded when a bad cluster WAS found. If nothing looked bad at any swept
# resolution, don't scaffold this step at all — set
# config.bad_cluster_removal.no_bad_cluster: true instead and skip straight to
# step8_seed (see SKILL.md section 5's no-bad-cluster shortcut).
RESOLUTION <- config$bad_cluster_removal$resolution
CLUSTER_RM <- config$bad_cluster_removal$cluster_rm
stopifnot(!is.null(RESOLUTION), !is.null(CLUSTER_RM))

cat("Loading obj5_UMAP...\n")
obj5 <- qs_read(here("step5_batch_effect", "obj5_UMAP.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj5), "cells\n")

res_col <- glue("RNA_snn_res.{RESOLUTION}")
keep    <- obj5@meta.data[[res_col]] != CLUSTER_RM
cat("Removing cluster", CLUSTER_RM, "at res", RESOLUTION, ":", sum(!keep), "cells removed\n")

obj6 <- obj5[, keep]
cat("Remaining:", ncol(obj6), "cells\n")

qs_save(obj6, here("step6_rm_badcl", "obj6_rm_badcl.qs2"), nthreads = n_cores)
cat("Saved: step6_rm_badcl/obj6_rm_badcl.qs2\n")
cat("\n=== STEP 6_RM_BADCL PROCESS.R COMPLETE ===\n")
