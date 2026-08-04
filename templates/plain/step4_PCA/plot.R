# conda run -n <conda.env_name> Rscript step4_PCA/scripts/plot.R
source(here::here("init.R"))

# Directories ----------
figs_dir   <- here("step4_PCA", "figs")
tables_dir <- here("step4_PCA", "tables")
dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Load ----------
cat("Loading obj4_PCA...\n")
obj <- qs_read(here("step4_PCA", "obj4_PCA.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj), "cells x", nrow(obj), "genes\n")

# tbl1: variance explained ----------
pca_stdev   <- obj[["pca"]]@stdev
pca_var     <- pca_stdev^2
pca_var_pct <- pca_var / sum(pca_var) * 100
var_df <- data.frame(
  PC          = seq_along(pca_stdev),
  stdev       = pca_stdev,
  var_pct     = round(pca_var_pct, 3),
  cum_var_pct = round(cumsum(pca_var_pct), 3)
)
write_csv(var_df, file.path(tables_dir, "tbl1_variance_explained_by_pc.csv"))
cat("Saved: tables/tbl1_variance_explained_by_pc.csv\n")
print(head(var_df, 15))

# fig1: elbow plot ----------
p <- ElbowPlot(obj, ndims = config$pca$npcs)
ggsave(file.path(figs_dir, "fig1_elbow.png"), p,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig1_elbow.png\n")

cat("\n=== STEP 4 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: don't run step5 yet.\n")
cat("Open figs/fig1_elbow.png. Each point is one principal component (PC); the y-axis is how much\n")
cat("of the cell-to-cell variation that PC explains. Look for the 'elbow' — the point where adding\n")
cat("more PCs stops adding much new information and the curve flattens out. Also check\n")
cat("tables/tbl1_variance_explained_by_pc.csv (cum_var_pct column) for the exact numbers.\n")
cat("Tell Claude how many PCs to use downstream (config.clustering.n_dims — currently",
    config$clustering$n_dims, "). Note this is a DIFFERENT number from config.pca.npcs (",
    config$pca$npcs, "),\n")
cat("which is just how many PCs were computed here (generous, fixed) — n_dims is how many of\n")
cat("them get USED for neighbor-finding/clustering starting in step5. Once you confirm and update\n")
cat("config.yaml, tell Claude to write step5_batch_effect (or, if this is a single-lane run,\n")
cat("step7_res directly — step5/step6 exist only to check and remove cross-lane batch effects).\n")
