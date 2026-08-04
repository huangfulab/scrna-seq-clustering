source(here::here("init.R"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))

# Load ----------
stopifnot(file.exists(here("step4_filter", "obj3_filter.qs2")))
cat("Loading obj3_filter...\n")
obj <- qs_read(here("step4_filter", "obj3_filter.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj), "cells x", nrow(obj), "genes\n")
cat("orig.ident:\n"); print(table(obj$orig.ident))

# NormalizeData ----------
cat("\n[1] NormalizeData...\n")
plan("multicore", workers = n_cores)
options(future.globals.maxSize = Inf)
obj <- NormalizeData(obj, verbose = FALSE)
plan("sequential")

# FindVariableFeatures ----------
cat(glue("\n[2] FindVariableFeatures (nfeatures = {config$pca$hvg_nfeatures})...\n"))
obj <- FindVariableFeatures(obj, nfeatures = config$pca$hvg_nfeatures, verbose = FALSE)
cat("Top 10 HVGs:", paste(head(VariableFeatures(obj), 10), collapse = ", "), "\n")

# CellCycleScoring ----------
cat("\n[3] CellCycleScoring...\n")
s_genes   <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes
obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes,
                         set.ident = FALSE)
cat("Phase distribution:\n"); print(table(obj$Phase))

# ScaleData ----------
vars_to_regress <- unlist(config$pca$vars_to_regress)
cat(glue("\n[4] ScaleData (vars.to.regress: {paste(vars_to_regress, collapse = ', ')})...\n"))
plan("multicore", workers = n_cores)
obj <- ScaleData(obj,
                 features        = VariableFeatures(obj),
                 vars.to.regress = vars_to_regress,
                 verbose         = TRUE)
plan("sequential")

# RunPCA ----------
cat(glue("\n[5] RunPCA (npcs = {config$pca$npcs})...\n"))
set.seed(config$pca$seed)
obj <- RunPCA(obj, npcs = config$pca$npcs, verbose = FALSE)
cat("Top genes PC1-2:\n")
print(obj[["pca"]], dims = 1:2, nfeatures = 5)

# Save ----------
cat("\nSaving obj4_PCA...\n")
obj4_path <- here("step5_PCA", "obj4_PCA.qs2")
qs_save(obj, file = obj4_path, nthreads = n_cores)
cat(glue("Saved: step5_PCA/obj4_PCA.qs2 ({round(file.size(obj4_path)/1e6, 1)} MB)\n"))

# Verification ----------
cat("\n=== Verification ===\n")
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("Layers:", paste(Layers(obj), collapse = ", "), "\n")
cat("\n=== STEP 5 COMPLETE — run step5_PCA/scripts/plot.R, then inspect elbow plot ===\n")
