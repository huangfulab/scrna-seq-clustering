source(here::here("init.R"))

# ---- Directories ----------
# NOTE: this dir is created (matching the original pipeline's step4_filter)
# but the merged object below is saved directly under step3_filter/, not
# inside it — preserved as-is for fidelity to the original.
dir.create(here("step3_filter", "obj"), recursive = TRUE, showWarnings = FALSE)

# ---- Per-lane filter ----------
# NOTE: the plain profile has no guide capture, so the original pipeline's
# MIN_MOI/MAX_MOI filter and f_assign (guide-origin restriction to NC/OEG
# cells only) are entirely absent here — only nCount/nFeature/pct.mt/doublet
# filters apply.
cat("Processing", length(lane_names), "lanes...\n")

obj_list <- map(lane_names, function(lane) {
  cat("\n--- Lane:", lane, "---\n")
  obj   <- qs_read(here("step2_doublet", "obj1_doublet", glue("{lane}.qs2")), nthreads = n_cores)
  n_pre <- ncol(obj)

  f_ncount  <- obj$nCount_RNA   <= config$qc$min_ncount
  f_nfeat   <- obj$nFeature_RNA <= config$qc$min_nfeat
  f_mt      <- obj$pct.mt       >= config$qc$max_pct_mt
  f_doublet <- is.na(obj$scDblFinder.class) | obj$scDblFinder.class != "singlet"

  keep   <- !f_ncount & !f_nfeat & !f_mt & !f_doublet
  n_post <- sum(keep)
  cat(glue("  Pre: {n_pre} | Post: {n_post} ({round(n_post/n_pre*100,1)}%)\n"))

  obj_f <- obj[, keep]
  obj_f <- RenameCells(obj_f, add.cell.id = lane)
  obj_f
})
names(obj_list) <- lane_names

# ---- Merge ----------
cat("\n=== Merging", length(obj_list), "lanes ===\n")
obj3_filter <- merge(obj_list[[1]], y = obj_list[-1], merge.data = TRUE)
cat("Merged:", ncol(obj3_filter), "cells x", nrow(obj3_filter), "genes\n")

qs_save(obj3_filter, file = here("step3_filter", "obj3_filter.qs2"), nthreads = n_cores)
cat("Saved: step3_filter/obj3_filter.qs2\n")

cat("\n=== STEP 3 COMPLETE ===\n")
cat("Run step3_filter/plot.R to generate figures and tables.\n")
