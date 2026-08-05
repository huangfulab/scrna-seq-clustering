# conda run -n <conda.env_name> Rscript step6_rm_badcl/scripts/plot.R
source(here::here("init.R"))

# Must match the values used in step6_rm_badcl/scripts/process.R
RESOLUTION <- config$bad_cluster_removal$resolution
CLUSTER_RM <- config$bad_cluster_removal$cluster_rm
stopifnot(!is.null(RESOLUTION), !is.null(CLUSTER_RM))

figs_dir <- here("step6_rm_badcl", "figs")
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

# Load ----------
cat("Loading obj5_UMAP (pre-removal)...\n")
obj5 <- qs_read(here("step5_batch_effect", "obj5_UMAP.qs2"), nthreads = n_cores)
cat("Loading obj6_rm_badcl (post-removal)...\n")
obj6 <- qs_read(here("step6_rm_badcl", "obj6_rm_badcl.qs2"), nthreads = n_cores)
cat("Pre:", ncol(obj5), "| Post:", ncol(obj6), "cells\n")

# NOTE: the original pipeline's fig1 (MOI pre/post removal), fig2 (cells per
# OEG target), and fig3 (coverage barplot) are all guide/MOI-specific and are
# dropped entirely here — no guide capture in the plain profile. Only the
# original fig5 (QC violins post-removal) survives, renumbered to fig1 since
# it is now the only figure this step produces.

# fig1 (was fig5): QC violins post-removal ----------
cat("[fig1] QC violins post-removal...\n")
ylim_max_list <- list(nCount_RNA = 30000, nFeature_RNA = 10000, pct.mt = 20)
qc_meta <- obj6@meta.data %>%
  select(nCount_RNA, nFeature_RNA, pct.mt, pct.rb) %>%
  mutate(lane = factor(obj6$orig.ident, levels = lane_names))

make_violin <- function(df, yvar, ylab, ylim_max = NA) {
  p <- ggplot(df, aes(x = lane, y = .data[[yvar]], fill = lane)) +
    geom_violin(alpha = 0.8, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", coef = 0) +
    scale_y_continuous(labels = scales::comma) +
    guides(fill = "none") +
    labs(title = ylab, x = NULL, y = yvar)
  if (!is.na(ylim_max)) p <- p + coord_cartesian(ylim = c(0, ylim_max))
  p
}

p <- wrap_plots(
  make_violin(qc_meta, "nCount_RNA",   "UMI count",        ylim_max_list$nCount_RNA),
  make_violin(qc_meta, "nFeature_RNA", "Gene count",        ylim_max_list$nFeature_RNA),
  make_violin(qc_meta, "pct.mt",       "% mt",              ylim_max_list$pct.mt),
  make_violin(qc_meta, "pct.rb",       "% ribosomal genes"),
  nrow = 2, axis_titles = "collect"
) + plot_annotation(title = glue("QC metrics per lane (post cluster-{CLUSTER_RM} removal)"))
# Fixed 2x2 metric grid — lanes sit on the x-axis WITHIN each panel, not as
# separate panels. Width = ncol panels * (fixed per-panel margin for axis/
# labels + a per-lane slot width * n lanes) — the frame content scales
# purely proportionally with lane count; only the margin term is fixed.
# Reproduces the original validated run's width=14 exactly at n=8 lanes.
qc_violin_margin_per_panel <- 1
qc_violin_width_per_lane   <- 0.75
qc_violin_ncol             <- 2
qc_violin_width <- qc_violin_ncol * (qc_violin_margin_per_panel + qc_violin_width_per_lane * length(lane_names))
ggsave(file.path(figs_dir, "fig1_qc_removal.png"), p, width = qc_violin_width, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig1_qc_removal.png\n")

cat("\n=== STEP 6_RM_BADCL PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: confirm this step finished and its outputs look reasonable before telling Claude to\n")
cat("write step7_res.\n")
cat("Open figs/fig1_qc_removal.png — the same 4-panel QC violin as earlier steps, now showing\n")
cat("only cells that remain after removing cluster", CLUSTER_RM, "(at resolution", RESOLUTION,
    "). Confirm it still\n")
cat("looks reasonable (no new problem visible) and that the cell count drop (see the process.R\n")
cat("console output: Pre vs Post) matches what you expected. No new config value is needed —\n")
cat("just confirm, then tell Claude to write step7_res.\n")
