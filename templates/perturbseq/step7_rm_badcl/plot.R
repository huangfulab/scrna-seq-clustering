# conda run -n <conda.env_name> Rscript step7_rm_badcl/scripts/plot.R
source(here::here("init.R"))

# Must match config.bad_cluster_removal.* used in step7_rm_badcl/scripts/process.R
stopifnot(!is.null(config$bad_cluster_removal$resolution), !is.null(config$bad_cluster_removal$cluster_rm))
RESOLUTION <- config$bad_cluster_removal$resolution
CLUSTER_RM <- as.character(config$bad_cluster_removal$cluster_rm)

# Cosmetic display caps / thresholds — not user-tunable QC cutoffs.
MOI_CAP       <- 40L
GUIDE_CAP     <- 500L
GENE_CAP      <- 3000L
COVERAGE_THRS <- c(10L, 25L, 50L, 100L)

figs_dir <- here("step7_rm_badcl", "figs")
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

# Load ----------
cat("Loading obj5_UMAP (pre-removal)...\n")
obj5 <- qs_read(here("step6_batch_effect", "obj5_UMAP.qs2"), nthreads = n_cores)
cat("Loading obj6_rm_badcl (post-removal)...\n")
obj6 <- qs_read(here("step7_rm_badcl", "obj6_rm_badcl.qs2"), nthreads = n_cores)
cat("Pre:", ncol(obj5), "| Post:", ncol(obj6), "cells\n")

# fig1: MOI pre/post removal ----------
cat("[fig1] MOI pre/post removal...\n")
med_pre  <- median(obj5$MOI)
med_post <- median(obj6$MOI)

fig1_df <- bind_rows(
  as.data.frame(table(MOI_capped = pmin(obj5$MOI, MOI_CAP))) %>% mutate(filter = "Pre-removal"),
  as.data.frame(table(MOI_capped = pmin(obj6$MOI, MOI_CAP))) %>% mutate(filter = "Post-removal")
) %>%
  mutate(
    MOI_capped = factor(as.character(MOI_capped), levels = as.character(0:MOI_CAP)),
    filter     = factor(filter, levels = c("Pre-removal", "Post-removal"))
  ) %>%
  group_by(MOI_capped, filter) %>%
  summarise(Freq = sum(Freq), .groups = "drop")

p <- ggplot(fig1_df, aes(x = MOI_capped, y = Freq, fill = filter)) +
  geom_col(position = "identity", alpha = 1, linewidth = 0) +
  scale_fill_manual(
    values = c("Pre-removal" = pal[2], "Post-removal" = pal[5]),
    labels = c("Pre-removal"  = glue("Pre-removal\n(med = {med_pre})"),
               "Post-removal" = glue("Post-removal\n(med = {med_post})")),
    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_discrete(breaks = as.character(c(seq(0, 35, 5), MOI_CAP)),
                   labels = c(seq(0, 35, 5), "≥40")) +
  labs(title = glue("MOI distribution pre/post cluster-{CLUSTER_RM} removal"), x = "MOI", y = "N cells") +
  theme(legend.position = "top")
ggsave(file.path(figs_dir, "fig1_MOI_prepost.png"), p, width = 6, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig1_MOI_prepost.png\n")

# fig2: cells per target ----------
cat("[fig2] Cells per target...\n")
oeg_meta   <- obj6@meta.data[obj6$guide_origin == "OEG", ]
genes_vec  <- unlist(strsplit(na.omit(oeg_meta$target_genes), ";"))
gene_tbl   <- table(genes_vec)
gene_df    <- data.frame(target = names(gene_tbl), n_cells = as.integer(gene_tbl),
                         stringsAsFactors = FALSE)
guides_vec <- unlist(strsplit(na.omit(oeg_meta$guides), ";"))
guides_oeg <- guides_vec[guides_vec %in% oeg_set]
guide_tbl  <- table(guides_oeg)
guide_df   <- data.frame(target = names(guide_tbl), n_cells = as.integer(guide_tbl),
                         stringsAsFactors = FALSE)

med_guide <- median(guide_df$n_cells)
med_gene  <- median(gene_df$n_cells)

p_guide <- ggplot(guide_df %>% mutate(n_cells = pmin(n_cells, GUIDE_CAP)), aes(x = n_cells)) +
  geom_histogram(bins = 100, fill = pal[2], linewidth = 0) +
  geom_vline(xintercept = pmin(med_guide, GUIDE_CAP), linetype = "dashed",
             color = "black", linewidth = 0.7) +
  annotate("text", x = pmin(med_guide, GUIDE_CAP), y = Inf,
           label = glue("Median = {scales::comma(med_guide)}"),
           hjust = -0.1, vjust = 2.5, size = 3.5) +
  scale_x_continuous(breaks = seq(0, GUIDE_CAP, 100),
                     labels = c(seq(0, GUIDE_CAP - 100, 100), "≥ 500")) +
  labs(title = "Per guide", x = "N cells", y = "N guides")

p_gene <- ggplot(gene_df %>% mutate(n_cells = pmin(n_cells, GENE_CAP)), aes(x = n_cells)) +
  geom_histogram(bins = 100, fill = pal[5], linewidth = 0) +
  geom_vline(xintercept = pmin(med_gene, GENE_CAP), linetype = "dashed",
             color = "black", linewidth = 0.7) +
  annotate("text", x = pmin(med_gene, GENE_CAP), y = Inf,
           label = glue("Median = {scales::comma(med_gene)}"),
           hjust = -0.1, vjust = 2.5, size = 3.5) +
  scale_x_continuous(breaks = seq(0, GENE_CAP, 1000),
                     labels = c("0", "1,000", "2,000", "≥ 3,000")) +
  labs(title = "Per gene", x = "N cells", y = "N genes")

p <- wrap_plots(p_guide, p_gene) +
  plot_annotation(title = glue("Cells per OEG target (post cluster-{CLUSTER_RM} removal)"))
ggsave(file.path(figs_dir, "fig2_cells_per_target.png"), p, width = 10, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig2_cells_per_target.png\n")

# fig3: coverage barplot ----------
cat("[fig3] Coverage barplot...\n")
n_gene_denom  <- length(unique(oeg_gene_map))
n_guide_denom <- length(oeg_set)
fig3_df <- map(COVERAGE_THRS, function(thr) {
  bind_rows(
    data.frame(threshold = thr, type = "Gene",
               n_meet = sum(gene_df$n_cells >= thr), denom = n_gene_denom,
               stringsAsFactors = FALSE),
    data.frame(threshold = thr, type = "Guide",
               n_meet = sum(guide_df$n_cells >= thr), denom = n_guide_denom,
               stringsAsFactors = FALSE)
  )
}) %>%
  list_rbind() %>%
  mutate(
    pct       = round(n_meet / denom * 100, 1),
    label     = paste0(n_meet, "\n(", pct, "%)"),
    threshold = factor(threshold, levels = COVERAGE_THRS),
    type      = factor(type, levels = c("Guide", "Gene"))
  )

p <- ggplot(fig3_df, aes(x = threshold, y = pct, fill = type, label = label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(position = position_dodge(width = 0.7), vjust = -0.2, size = 3, lineheight = 0.9) +
  scale_fill_manual(values = c("Gene" = pal[5], "Guide" = pal[2]), name = NULL) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  labs(title = "Coverage: % guides/genes with >= N cells",
       x = "Minimum N cells threshold", y = "% of total") +
  theme(legend.position = "top")
ggsave(file.path(figs_dir, "fig3_coverage_barplot.png"), p, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: figs/fig3_coverage_barplot.png\n")

# fig5: QC violins post-removal ----------
cat("[fig5] QC violins post-removal...\n")
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
ggsave(file.path(figs_dir, "fig5_qc_removal.png"), p, width = qc_violin_width, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig5_qc_removal.png\n")

cat("\n=== STEP 7_RM_BADCL PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: confirm this step finished and its outputs look reasonable before telling Claude to\n")
cat("write step8. No new value is needed — skim figs/fig5_qc_removal.png to confirm the removed\n")
cat("cluster is really gone and the remaining QC distributions still look reasonable.\n")
