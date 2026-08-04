# conda run -n <conda.env_name> Rscript step5_PCA/scripts/plot.R
source(here::here("init.R"))

# Directories ----------
figs_dir   <- here("step5_PCA", "figs")
tables_dir <- here("step5_PCA", "tables")
dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Load ----------
cat("Loading obj4_PCA...\n")
obj <- qs_read(here("step5_PCA", "obj4_PCA.qs2"), nthreads = n_cores)
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

cat("\n=== STEP 5 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: don't run step6 yet.\n")
cat("Open figs/fig1_elbow.png. Each point is one principal component (PC); the y-axis is how much\n")
cat("variance that PC explains. Look for the 'elbow' — the PC number where the curve flattens out\n")
cat("and adding more PCs stops adding much new signal (mostly noise beyond that point). Report that\n")
cat("PC number back to Claude as config.clustering.n_dims (current guess in config.yaml:",
    config$clustering$n_dims, "). This is different from config.pca.npcs (", config$pca$npcs,
    "), which is\n")
cat("just how many PCs were computed — n_dims is how many of those to actually use for clustering.\n")
cat("Once you confirm, step6_batch_effect will be written using that value.\n")
