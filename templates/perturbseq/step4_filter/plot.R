# conda run -n <conda.env_name> Rscript step4_filter/scripts/plot.R
source(here::here("init.R"))

library(ComplexUpset)

# ---- Parameters (cosmetic display caps — not user-tunable QC cutoffs) ----------
MOI_CAP       <- 40L
COVERAGE_THRS <- c(10L, 25L, 50L, 100L)

figs_dir   <- here("step4_filter", "figs")
tables_dir <- here("step4_filter", "tables")
dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read step3 objs (figs 1, 4 + filter summary) ----------
cat("Reading step3 per-lane objects...\n")

lane_data <- map(lane_names, function(lane) {
  cat(" ", lane, "\n")
  obj   <- qs_read(here("step3_doublet", "obj2_doublet", glue("{lane}.qs2")), nthreads = n_cores)
  n_pre <- ncol(obj)

  f_ncount  <- obj$nCount_RNA   <= config$qc$min_ncount
  f_nfeat   <- obj$nFeature_RNA <= config$qc$min_nfeat
  f_mt      <- obj$pct.mt       >= config$qc$max_pct_mt
  f_doublet <- is.na(obj$scDblFinder.class) | obj$scDblFinder.class != "singlet"
  f_moi     <- obj$MOI < config$qc$min_moi | obj$MOI > config$qc$max_moi
  f_assign  <- !(obj$guide_origin %in% c("NC", "OEG"))
  keep      <- !f_ncount & !f_nfeat & !f_mt & !f_doublet & !f_moi & !f_assign

  # fig1: MOI distribution pre/post
  pre_df  <- as.data.frame(table(MOI_capped = pmin(obj$MOI, MOI_CAP))) %>% mutate(filter = "Pre-filter")
  post_df <- as.data.frame(table(MOI_capped = pmin(obj$MOI[keep], MOI_CAP))) %>% mutate(filter = "Post-filter")
  moi_df  <- bind_rows(pre_df, post_df)

  # fig4: disqualification flags (cells failing >= 1 filter)
  disq_df <- data.frame(
    "nCount"     = as.integer(f_ncount),
    "nFeat"      = as.integer(f_nfeat),
    "pct.mt"     = as.integer(f_mt),
    "Doublet"    = as.integer(!is.na(obj$scDblFinder.class) & obj$scDblFinder.class == "doublet"),
    "MOI_low"    = as.integer(obj$MOI < config$qc$min_moi),
    "MOI_high"   = as.integer(obj$MOI > config$qc$max_moi),
    "Mixed"      = as.integer(obj$guide_origin %in% "Mixed"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  disq_df <- disq_df[rowSums(disq_df) > 0, ]

  # filter summary row
  meta_post <- obj@meta.data[keep, ]
  summary_row <- data.frame(
    lane          = lane,
    pre_filter_n  = n_pre,
    post_filter_n = sum(keep),
    pct_retained  = round(sum(keep) / n_pre * 100, 1),
    OEG_n         = sum(meta_post$guide_origin == "OEG"),
    NC_n          = sum(meta_post$guide_origin == "NC"),
    stringsAsFactors = FALSE
  )

  list(moi_df = moi_df, disq_df = disq_df, summary_row = summary_row,
       moi_raw_pre = obj$MOI, moi_raw_post = obj$MOI[keep])
})
names(lane_data) <- lane_names

# ---- Read obj3_filter (figs 2, 3, 5) ----------
cat("Reading obj3_filter...\n")
obj3_filter <- qs_read(here("step4_filter", "obj3_filter.qs2"), nthreads = n_cores)

# ---- tbl1: filter summary ----------
filter_tbl <- bind_rows(map(lane_data, "summary_row"))
rownames(filter_tbl) <- NULL
write_csv(filter_tbl, file.path(tables_dir, "tbl1_filter_summary.csv"))
cat("Saved: tables/tbl1_filter_summary.csv\n")
print(filter_tbl)

# fig1: MOI pre/post ----------
cat("[fig1] MOI pre/post...\n")
fig1_df <- bind_rows(map(lane_data, "moi_df")) %>%
  mutate(
    MOI_capped = factor(as.character(MOI_capped), levels = as.character(0:MOI_CAP)),
    filter     = factor(filter, levels = c("Pre-filter", "Post-filter"))
  ) %>%
  group_by(MOI_capped, filter) %>%
  summarise(Freq = sum(Freq), .groups = "drop")

med_pre  <- median(unlist(map(lane_data, "moi_raw_pre")))
med_post <- median(unlist(map(lane_data, "moi_raw_post")))

p <- ggplot(fig1_df, aes(x = MOI_capped, y = Freq, fill = filter)) +
  geom_col(position = "identity", alpha = 1, linewidth = 0) +
  scale_fill_manual(
    values = c("Pre-filter" = pal[2], "Post-filter" = pal[5]),
    labels = c("Pre-filter"  = glue("Pre-filter\n(med = {med_pre})"),
               "Post-filter" = glue("Post-filter\n(med = {med_post})")),
    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(breaks = as.character(c(seq(0, 35, 5), MOI_CAP)),
                   labels = c(seq(0, 35, 5), paste0("≥40"))) +
  labs(title = "MOI distribution pre/post filter", x = "MOI", y = "N cells") +
  theme(legend.position = "top")

ggsave(file.path(figs_dir, "fig1_MOI_prepost.png"), p, width = 6, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig1_MOI_prepost.png\n")

# fig2: cells per target ----------
cat("[fig2] Cells per target...\n")
GUIDE_CAP <- 500L
GENE_CAP  <- 3000L

oeg_meta   <- obj3_filter@meta.data[obj3_filter$guide_origin == "OEG", ]
genes_vec  <- unlist(strsplit(na.omit(oeg_meta$target_genes), ";"))
gene_tbl   <- table(genes_vec)
gene_df    <- data.frame(type = "Per gene",  target = names(gene_tbl),
                         n_cells = as.integer(gene_tbl), stringsAsFactors = FALSE)
guides_vec <- unlist(strsplit(na.omit(oeg_meta$guides), ";"))
guides_oeg <- guides_vec[guides_vec %in% oeg_set]
guide_tbl  <- table(guides_oeg)
guide_df   <- data.frame(type = "Per guide", target = names(guide_tbl),
                         n_cells = as.integer(guide_tbl), stringsAsFactors = FALSE)
fig2_df    <- bind_rows(gene_df, guide_df) %>%
  mutate(type = factor(type, levels = c("Per guide", "Per gene")))

med_guide <- median(guide_df$n_cells)
med_gene  <- median(gene_df$n_cells)

p_guide <- ggplot(filter(fig2_df, type == "Per guide") %>% mutate(n_cells = pmin(n_cells, GUIDE_CAP)),
                  aes(x = n_cells)) +
  geom_histogram(bins = 100, fill = pal[2], linewidth = 0) +
  geom_vline(xintercept = pmin(med_guide, GUIDE_CAP), linetype = "dashed", color = "black", linewidth = 0.7) +
  annotate("text", x = pmin(med_guide, GUIDE_CAP), y = Inf,
           label = glue("Median = {scales::comma(med_guide)}"),
           hjust = -0.1, vjust = 2.5, size = 3.5) +
  scale_x_continuous(breaks = seq(0, GUIDE_CAP, 100),
                     labels = c(seq(0, GUIDE_CAP - 100, 100), paste0("≥ ", scales::comma(GUIDE_CAP)))) +
  labs(title = "Per guide", x = "N cells", y = "N guides")

p_gene <- ggplot(filter(fig2_df, type == "Per gene") %>% mutate(n_cells = pmin(n_cells, GENE_CAP)),
                 aes(x = n_cells)) +
  geom_histogram(bins = 100, fill = pal[5], linewidth = 0) +
  geom_vline(xintercept = pmin(med_gene, GENE_CAP), linetype = "dashed", color = "black", linewidth = 0.7) +
  annotate("text", x = pmin(med_gene, GENE_CAP), y = Inf,
           label = glue("Median = {scales::comma(med_gene)}"),
           hjust = -0.1, vjust = 2.5, size = 3.5) +
  scale_x_continuous(breaks = seq(0, GENE_CAP, 1000),
                     labels = c("0", "1,000", "2,000", paste0("≥ ", scales::comma(GENE_CAP)))) +
  labs(title = "Per gene", x = "N cells", y = "N genes")

p <- wrap_plots(p_guide, p_gene) +
  plot_annotation(title = "Cells per OEG target (post-filter)")
ggsave(file.path(figs_dir, "fig2_cells_per_target.png"), p, width = 10, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig2_cells_per_target.png\n")

# fig3: coverage barplot ----------
cat("[fig3] Coverage barplot...\n")
n_gene_denom  <- length(unique(oeg_gene_map))
n_guide_denom <- length(oeg_set)
fig3_df <- map(COVERAGE_THRS, function(thr) {
  bind_rows(
    data.frame(threshold = thr, type = "Per gene",
               n_meet = sum(gene_df$n_cells  >= thr), denom = n_gene_denom,
               stringsAsFactors = FALSE),
    data.frame(threshold = thr, type = "Per guide",
               n_meet = sum(guide_df$n_cells >= thr), denom = n_guide_denom,
               stringsAsFactors = FALSE)
  )
}) %>%
  list_rbind() %>%
  mutate(
    pct       = round(n_meet / denom * 100, 1),
    label     = paste0(n_meet, "\n(", pct, "%)"),
    threshold = factor(threshold, levels = COVERAGE_THRS),
    type      = factor(type, levels = c("Per guide", "Per gene"))
  )

p <- ggplot(fig3_df, aes(x = threshold, y = pct, fill = type, label = label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(position = position_dodge(width = 0.7), vjust = -0.2, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = c("Per gene" = pal[5], "Per guide" = pal[2]), name = NULL) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  labs(title = "Coverage: % guides/genes with >= N cells",
       x = "Minimum N cells threshold", y = "% of total") +
  theme(legend.position = "top")
ggsave(file.path(figs_dir, "fig3_coverage_barplot.png"), p, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig3_coverage_barplot.png\n")

# fig4: UpSet (disqualification) ----------
cat("[fig4] UpSet...\n")
fig4_df <- bind_rows(map(lane_data, "disq_df"))
rownames(fig4_df) <- NULL
fig4_df <- fig4_df %>%
  rename("MOI = 0"              = "MOI_low",
         !!glue("MOI > {config$qc$max_moi}") := "MOI_high") %>%
  mutate(across(everything(), as.logical))
n_disq    <- nrow(fig4_df)
set_cols  <- names(fig4_df)

bw_no_grid <- theme_bw(base_family = "Arial") + theme(panel.grid = element_blank())

p <- upset(
  fig4_df,
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
ggsave(file.path(figs_dir, "fig4_upset.png"), p, width = 9, height = 6, dpi = 300, bg = "white")
cat("Saved: figs/fig4_upset.png\n")

# fig5: QC violins post-filter ----------
cat("[fig5] QC violins post-filter...\n")
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

ggsave(file.path(figs_dir, "fig5_qc_filter.png"), p, width = 14, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig5_qc_filter.png\n")

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

cat("\n=== STEP 4 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: confirm this step finished and its outputs look reasonable before telling Claude to\n")
cat("write step5. No new value is needed — skim figs/fig4_upset.png (which filters overlap most)\n")
cat("and figs/fig5_qc_filter.png (post-filter QC looks clean) as a sanity check.\n")
