source(here::here("init.R"))

# ---- Directories ----------
dir.create(here("step4_filter", "obj"), recursive = TRUE, showWarnings = FALSE)

# ---- Per-lane filter ----------
cat("Processing", length(lane_names), "lanes...\n")

obj_list <- map(lane_names, function(lane) {
  cat("\n--- Lane:", lane, "---\n")
  obj   <- qs_read(here("step3_doublet", "obj2_doublet", glue("{lane}.qs2")), nthreads = n_cores)
  n_pre <- ncol(obj)

  f_ncount  <- obj$nCount_RNA   <= config$qc$min_ncount
  f_nfeat   <- obj$nFeature_RNA <= config$qc$min_nfeat
  f_mt      <- obj$pct.mt       >= config$qc$max_pct_mt
  f_doublet <- is.na(obj$scDblFinder.class) | obj$scDblFinder.class != "singlet"
  f_moi     <- obj$MOI < config$qc$min_moi | obj$MOI > config$qc$max_moi
  f_assign  <- !(obj$guide_origin %in% c("NC", "OEG"))

  keep   <- !f_ncount & !f_nfeat & !f_mt & !f_doublet & !f_moi & !f_assign
  n_post <- sum(keep)
  cat(glue("  Pre: {n_pre} | Post: {n_post} ({round(n_post/n_pre*100,1)}%)\n"))

  obj_f <- obj[, keep]
  obj_f <- RenameCells(obj_f, add.cell.id = lane)
  cat("  guide_origin:\n"); print(table(obj_f$guide_origin))
  obj_f
})
names(obj_list) <- lane_names

# ---- Merge ----------
cat("\n=== Merging", length(obj_list), "lanes ===\n")
obj3_filter <- merge(obj_list[[1]], y = obj_list[-1], merge.data = TRUE)
cat("Merged:", ncol(obj3_filter), "cells x", nrow(obj3_filter), "genes\n")
cat("guide_origin:\n"); print(table(obj3_filter$guide_origin))

qs_save(obj3_filter, file = here("step4_filter", "obj3_filter.qs2"), nthreads = n_cores)
cat("Saved: step4_filter/obj3_filter.qs2\n")

cat("\n=== STEP 4 COMPLETE ===\n")
cat("Run step4_filter/scripts/plot.R to generate figures and tables.\n")

cat("\n---\n")
cat("STOP: confirm this step finished and its outputs look reasonable before telling Claude to\n")
cat("write step5. There is no new value to decide here — the QC/MOI cutoffs were already locked\n")
cat("in at the step2/step3 checkpoints. Check the printed 'Pre: ... | Post: ...' cell counts per\n")
cat("lane above look sane (not near-zero, not near-100%), then run plot.R and skim figs/fig4_upset.png\n")
cat("and figs/fig5_qc_filter.png for anything alarming.\n")
