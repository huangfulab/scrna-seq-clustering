# conda run -n <conda.env_name> Rscript step3_filter/scripts/plot.R
source(here::here("init.R"))

library(ComplexUpset)

figs_dir   <- here("step3_filter", "figs")
tables_dir <- here("step3_filter", "tables")
dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# NOTE: renumbering vs. the original pipeline's step4_filter — that step had
# 5 figs (fig1 MOI pre/post, fig2 cells per target, fig3 coverage barplot,
# fig4 UpSet, fig5 QC violins post-filter). fig1/fig2/fig3 all depend on
# guide/MOI concepts that don't exist in the plain profile and are dropped
# entirely. The original fig4 (UpSet of disqualification reasons) becomes
# fig1 here, and fig5 (QC violins post-filter) becomes fig2.

# ---- Read step2 objs (fig1 UpSet + filter summary) ----------
cat("Reading step2 per-lane objects...\n")

lane_data <- map(lane_names, function(lane) {
  cat(" ", lane, "\n")
  obj   <- qs_read(here("step2_doublet", "obj1_doublet", glue("{lane}.qs2")), nthreads = n_cores)
  n_pre <- ncol(obj)

  f_ncount  <- obj$nCount_RNA   <= config$qc$min_ncount
  f_nfeat   <- obj$nFeature_RNA <= config$qc$min_nfeat
  f_mt      <- obj$pct.mt       >= config$qc$max_pct_mt
  f_doublet <- is.na(obj$scDblFinder.class) | obj$scDblFinder.class != "singlet"
  keep      <- !f_ncount & !f_nfeat & !f_mt & !f_doublet

  # fig1 (was fig4): disqualification flags (cells failing >= 1 filter)
  # NOTE: the original pipeline's disq_df also had MOI_low/MOI_high/Mixed
  # columns (guide-specific) — dropped here, only QC + doublet flags remain.
  disq_df <- data.frame(
    "nCount"  = as.integer(f_ncount),
    "nFeat"   = as.integer(f_nfeat),
    "pct.mt"  = as.integer(f_mt),
    "Doublet" = as.integer(!is.na(obj$scDblFinder.class) & obj$scDblFinder.class == "doublet"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  disq_df <- disq_df[rowSums(disq_df) > 0, ]

  # filter summary row
  summary_row <- data.frame(
    lane          = lane,
    pre_filter_n  = n_pre,
    post_filter_n = sum(keep),
    pct_retained  = round(sum(keep) / n_pre * 100, 1),
    stringsAsFactors = FALSE
  )

  list(disq_df = disq_df, summary_row = summary_row)
})
names(lane_data) <- lane_names

# ---- Read obj3_filter (fig2 QC violins) ----------
cat("Reading obj3_filter...\n")
obj3_filter <- qs_read(here("step3_filter", "obj3_filter.qs2"), nthreads = n_cores)

# ---- tbl1: filter summary ----------
filter_tbl <- bind_rows(map(lane_data, "summary_row"))
rownames(filter_tbl) <- NULL
write_csv(filter_tbl, file.path(tables_dir, "tbl1_filter_summary.csv"))
cat("Saved: tables/tbl1_filter_summary.csv\n")
print(filter_tbl)

# fig1: UpSet (disqualification) ----------
cat("[fig1] UpSet...\n")
fig1_df <- bind_rows(map(lane_data, "disq_df"))
rownames(fig1_df) <- NULL
fig1_df <- fig1_df %>%
  mutate(across(everything(), as.logical))
n_disq   <- nrow(fig1_df)
set_cols <- names(fig1_df)

bw_no_grid <- theme_bw(base_family = "Arial") + theme(panel.grid = element_blank())

p <- upset(
  fig1_df,
  intersect              = set_cols,
  min_size               = 100,
  sort_sets              = "descending",
  sort_intersections_by  = "cardinality",
  themes = upset_modify_themes(list(
    "Intersection size" = bw_no_grid + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ),
    intersections_matrix = bw_no_grid + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ),
    overall_sizes = bw_no_grid
  )),
  base_annotations = list(
    "Intersection size" = intersection_size(counts = FALSE) +
      stat_summary(
        aes(label = after_stat(scales::comma(y))),
        geom = "text", fun = "sum",
        angle = 45, hjust = 0, vjust = -0.5, size = 3
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  ),
  set_sizes = (
    upset_set_size(geom = geom_bar(width = 0.6)) +
      scale_y_continuous(
        labels = function(x) x / 10000,
        name   = "Set size (10,000)",
        trans  = "reverse"
      ) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.title.y = element_blank())
  ),
  height_ratio = 0.4,
  width_ratio  = 0.25
) +
  plot_annotation(
    title = glue(
      "Disqualified cells by filters\n",
      "(n = {scales::comma(n_disq)} cells, {length(lane_names)} lanes combined, intersection min size = 100 cells)"
    ),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 13))
  )
ggsave(file.path(figs_dir, "fig1_upset.png"), p, width = 9, height = 6, dpi = 300, bg = "white")
cat("Saved: figs/fig1_upset.png\n")

# fig2: QC violins post-filter ----------
cat("[fig2] QC violins post-filter...\n")
cutoff_list   <- list(nCount_RNA = config$qc$min_ncount, nFeature_RNA = config$qc$min_nfeat, pct.mt = config$qc$max_pct_mt)
ylim_max_list <- list(nCount_RNA = 30000, nFeature_RNA = 10000, pct.mt = 20)

qc_meta <- obj3_filter@meta.data %>%
  select(nCount_RNA, nFeature_RNA, pct.mt, pct.rb) %>%
  mutate(lane = factor(obj3_filter$orig.ident, levels = lane_names))

make_violin <- function(df, yvar, ylab, ylim_max = NA, cutoff = NULL) {
  p <- ggplot(df, aes(x = lane, y = .data[[yvar]], fill = lane)) +
    geom_violin(alpha = 0.8, scale = "width") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", coef = 0) +
    scale_y_continuous(labels = scales::comma) +
    guides(fill = "none") +
    labs(title = ylab, x = NULL, y = yvar)
  if (!is.na(ylim_max)) p <- p + coord_cartesian(ylim = c(0, ylim_max))
  if (!is.null(cutoff)) p <- p + geom_hline(yintercept = cutoff, linetype = "dashed",
                                             color = "red", linewidth = 0.7)
  p
}

p <- wrap_plots(
  make_violin(qc_meta, "nCount_RNA",   glue("UMI count ({cutoff_list$nCount_RNA})"),   ylim_max_list$nCount_RNA,   cutoff_list$nCount_RNA),
  make_violin(qc_meta, "nFeature_RNA", glue("Gene count ({cutoff_list$nFeature_RNA})"), ylim_max_list$nFeature_RNA, cutoff_list$nFeature_RNA),
  make_violin(qc_meta, "pct.mt",       glue("% mt ({cutoff_list$pct.mt})"),             ylim_max_list$pct.mt,       cutoff_list$pct.mt),
  make_violin(qc_meta, "pct.rb",       "% ribosomal genes"),
  nrow = 2, axis_titles = "collect"
) + plot_annotation(title = "QC metrics per lane (post-filter)")

# Fixed 2x2 metric grid — lanes sit on the x-axis WITHIN each panel, not as
# separate panels. Width = ncol panels * (fixed per-panel margin for axis/
# labels + a per-lane slot width * n lanes) — the frame content scales
# purely proportionally with lane count; only the margin term is fixed.
# Reproduces the original validated run's width=14 exactly at n=8 lanes.
qc_violin_margin_per_panel <- 1
qc_violin_width_per_lane   <- 0.75
qc_violin_ncol             <- 2
qc_violin_width <- qc_violin_ncol * (qc_violin_margin_per_panel + qc_violin_width_per_lane * length(lane_names))
ggsave(file.path(figs_dir, "fig2_qc_filter.png"), p, width = qc_violin_width, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig2_qc_filter.png\n")

# tbl2: post-filter QC summary ----------
qc_summary <- qc_meta %>%
  group_by(lane) %>%
  summarise(
    n_cells          = n(),
    nCount_median    = round(median(nCount_RNA,   na.rm = TRUE), 0),
    nFeature_median  = round(median(nFeature_RNA, na.rm = TRUE), 0),
    pct_mt_median    = round(median(pct.mt,       na.rm = TRUE), 2),
    pct_rb_median    = round(median(pct.rb,       na.rm = TRUE), 2),
    .groups = "drop"
  )
write_csv(qc_summary, file.path(tables_dir, "tbl2_qc_summary.csv"))
cat("Saved: tables/tbl2_qc_summary.csv\n")

cat("\n=== STEP 3 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: confirm this step finished and its outputs look reasonable before telling Claude to\n")
cat("write step4_PCA.\n")
cat("Open figs/fig1_upset.png — it shows how many cells were dropped for each REASON (low UMI,\n")
cat("low gene count, high %mt, called a doublet) and how those reasons overlap. Open\n")
cat("figs/fig2_qc_filter.png — the same 4-panel QC violin as step1, but now only showing cells\n")
cat("that passed filtering (should look much cleaner/tighter than the step1 version). Check\n")
cat("tables/tbl1_filter_summary.csv for how many cells were retained per lane (pct_retained).\n")
cat("No new config value is needed here — just confirm the retained cell counts look reasonable\n")
cat("(not near-zero, not near-100%) before telling Claude to write step4_PCA.\n")
