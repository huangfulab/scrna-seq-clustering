# conda run -n <conda.env_name> Rscript step2_assign/scripts/plot.R
source(here::here("init.R"))

figs_dir <- here("step2_assign", "figs")
tables_dir <- here("step2_assign", "tables")
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

cutoff_list <- list(
  nCount_RNA   = config$qc$min_ncount,
  nFeature_RNA = config$qc$min_nfeat,
  pct.mt       = config$qc$max_pct_mt
)
ylim_max_list <- list(
  nCount_RNA   = 30000,
  nFeature_RNA = 10000,
  pct.mt       = 20
)

# ---- Derive plot data from per-lane objects ----------------------------------
cat("Reading per-lane objects...\n")
qc_meta <- map(lane_names, function(l) {
  obj <- qs_read(here("step2_assign", "obj1_assign", glue("{l}.qs2")), nthreads = n_cores)
  cat(" ", l, "\n")
  cbind(lane = l, obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "pct.mt", "pct.rb")])
}) %>%
  list_rbind() %>%
  mutate(lane = factor(lane, levels = lane_names))

# fill = lane → ggplot2 default hue palette = Seurat VlnPlot default
# coord_cartesian zooms without clipping data (preserves violin shape)
make_violin <- function(df, yvar, ylab, ylim_max = NA, cutoff = NULL) {
  p <- ggplot(df, aes(x = lane, y = .data[[yvar]], fill = lane)) +
    geom_violin(alpha = 0.8, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", coef = 0) +
    scale_y_continuous(labels = scales::comma) +
    guides(fill = "none") +
    labs(title = ylab, x = NULL, y = yvar)
  if (!is.na(ylim_max)) p <- p + coord_cartesian(ylim = c(0, ylim_max))
  if (!is.null(cutoff)) {
    p <- p + geom_hline(
      yintercept = cutoff, linetype = "dashed",
      color = "red", linewidth = 0.7
    )
  }
  p
}

p_ncount <- make_violin(qc_meta, "nCount_RNA", glue("UMI count ({cutoff_list$nCount_RNA})"),
  ylim_max = ylim_max_list$nCount_RNA,
  cutoff = cutoff_list$nCount_RNA
)
p_nfeature <- make_violin(qc_meta, "nFeature_RNA", glue("Gene count ({cutoff_list$nFeature_RNA})"),
  ylim_max = ylim_max_list$nFeature_RNA,
  cutoff = cutoff_list$nFeature_RNA
)
p_pctmt <- make_violin(qc_meta, "pct.mt", glue("% mitochondrial genes ({cutoff_list$pct.mt})"),
  ylim_max = ylim_max_list$pct.mt,
  cutoff = cutoff_list$pct.mt
)
p_pctrb <- make_violin(qc_meta, "pct.rb", "% ribosomal genes",
  cutoff = NULL
)

p <- wrap_plots(p_ncount, p_nfeature, p_pctmt, p_pctrb,
  nrow = 2,
  axis_titles = "collect"
) + plot_annotation(title = "QC metrics per lane")

# Fixed 2x2 metric grid — lanes sit on the x-axis WITHIN each panel, not as
# separate panels. Width = ncol panels * (fixed per-panel margin for axis/
# labels + a per-lane slot width * n lanes) — the frame content scales
# purely proportionally with lane count; only the margin term is fixed.
# Reproduces the original validated run's width=14 exactly at n=8 lanes.
qc_violin_margin_per_panel <- 1
qc_violin_width_per_lane   <- 0.75
qc_violin_ncol             <- 2
qc_violin_width <- qc_violin_ncol * (qc_violin_margin_per_panel + qc_violin_width_per_lane * length(lane_names))
ggsave(file.path(figs_dir, "fig1_qc_violin.png"), p, width = qc_violin_width, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig1_qc_violin.png\n")

# ---- Summary table ----------------------------------------------------------
qc_summary <- qc_meta %>%
  group_by(lane) %>%
  summarise(
    n_cells = n(),
    nCount_median = round(median(nCount_RNA, na.rm = TRUE), 1),
    nCount_cutoff = cutoff_list$nCount_RNA,
    nCount_pass_n = sum(nCount_RNA >= cutoff_list$nCount_RNA, na.rm = TRUE),
    nCount_pass_pct = round(100 * mean(nCount_RNA >= cutoff_list$nCount_RNA, na.rm = TRUE), 1),
    nFeature_median = round(median(nFeature_RNA, na.rm = TRUE), 1),
    nFeature_cutoff = cutoff_list$nFeature_RNA,
    nFeature_pass_n = sum(nFeature_RNA >= cutoff_list$nFeature_RNA, na.rm = TRUE),
    nFeature_pass_pct = round(100 * mean(nFeature_RNA >= cutoff_list$nFeature_RNA, na.rm = TRUE), 1),
    pct_mt_median = round(median(pct.mt, na.rm = TRUE), 1),
    pct_mt_cutoff = cutoff_list$pct.mt,
    pct_mt_pass_n = sum(pct.mt <= cutoff_list$pct.mt, na.rm = TRUE),
    pct_mt_pass_pct = round(100 * mean(pct.mt <= cutoff_list$pct.mt, na.rm = TRUE), 1),
    pct_rb_median = round(median(pct.rb, na.rm = TRUE), 1),
    .groups = "drop"
  )

write_csv(qc_summary, file.path(tables_dir, "tbl1_qc_summary.csv"))
cat("Saved: tables/tbl1_qc_summary.csv\n")

cat("\n=== STEP 2 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: don't run step3 yet.\n")
cat("Open figs/fig1_qc_violin.png. It shows, for each sample, how many RNA molecules (UMI count),\n")
cat("how many genes (Gene count), and what % of RNA is mitochondrial (% mt) each cell has.\n")
cat("The dashed red lines show the CURRENT cutoff guesses in config.yaml (qc.min_ncount =",
    config$qc$min_ncount, ", qc.min_nfeat =", config$qc$min_nfeat,
    ", qc.max_pct_mt =", config$qc$max_pct_mt, "). Decide if they look right (cells below/right\n")
cat("of a bad cutoff are usually dying or empty droplets) and tell Claude the final min UMI count /\n")
cat("min gene count / max %mt to use. Once you confirm and update config.yaml, re-run this plot.R —\n")
cat("fig1 will regenerate with the confirmed lines — before step3 is written.\n")
