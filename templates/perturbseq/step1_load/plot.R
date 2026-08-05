# conda run -n <conda.env_name> Rscript step1_load/scripts/plot.R
source(here::here("init.R"))

figs_dir <- here("step1_load", "figs")
tables_dir <- here("step1_load", "tables")
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

assignment_colors <- c(NC = pal[3], OEG = pal[2], Mixed = pal[1], Unassigned = "grey70")

# ---- Derive plot data — one lane at a time to keep peak memory low ----------
cat("Deriving plot data from per-lane objects...\n")
lane_data <- setNames(
  map(lane_names, function(l) {
    obj <- qs_read(here("step1_load", "obj0_load", glue("{l}.qs2")), nthreads = n_cores)
    crispr_mat <- GetAssayData(obj, assay = "CRISPR", layer = "counts")
    thresh_mats <- map(1:10, ~ crispr_mat >= .x)

    # fig1: non-zero UMI counts 1-10 per guide-cell pair
    crispr_nz <- summary(crispr_mat)
    umi_values <- as.integer(crispr_nz$x)
    plot1_df <- data.frame(
      umi = factor(umi_values[umi_values >= 1 & umi_values <= 10], levels = 1:10)
    ) %>% count(umi, .drop = FALSE)

    # fig2: guide origin counts per threshold (already in metadata)
    plot2_df <- map(1:10, function(t) {
      tbl <- table(factor(
        obj@meta.data[[glue("guide_origin_UMI{t}")]],
        levels = c("NC", "OEG", "Mixed", "Unassigned")
      ))
      data.frame(
        threshold = t, category = names(tbl), n_cells = as.integer(tbl),
        stringsAsFactors = FALSE
      )
    }) %>% list_rbind()

    # fig3: assigned guide recovery per threshold
    plot3_df <- map(1:10, function(t) {
      m <- thresh_mats[[t]]
      assigned <- rownames(m)[rowSums(m) > 0]
      det_nc <- assigned[assigned %in% nc_set]
      det_oeg <- assigned[assigned %in% oeg_set]
      data.frame(
        threshold = t,
        category = c("Total", "NC", "OEG"),
        n_assigned = c(length(det_nc) + length(det_oeg), length(det_nc), length(det_oeg)),
        n_total = c(n_all_total, n_nc_total, n_oeg_total),
        stringsAsFactors = FALSE
      )
    }) %>%
      list_rbind() %>%
      mutate(recovery = n_assigned / n_total * 100)

    cat(" ", l, "\n")
    list(plot1 = plot1_df, plot2 = plot2_df, plot3 = plot3_df)
  }),
  lane_names
)

lane_grid <- lane_grid_dims(length(lane_names))

# fig1: UMI distribution ----------
cat("[fig1] UMI distribution...\n")
p_list <- map(lane_names, function(l) {
  df <- lane_data[[l]]$plot1
  ggplot(df, aes(x = umi, y = n)) +
    geom_col(fill = pal[2]) +
    scale_y_log10(labels = scales::comma) +
    labs(title = l, x = "UMI per guide per cell", y = "Frequency")
})
p <- wrap_plots(p_list, nrow = lane_grid$nrow, axis_titles = "collect") +
  plot_annotation(title = "UMI distribution per lane (zoom range: 1-10)")
ggsave(file.path(figs_dir, "fig1_UMI_distribution.png"), p, width = lane_grid$width, height = lane_grid$height, dpi = 300, bg = "white")
cat("Saved: figs/fig1_UMI_distribution.png\n")

# fig2: Guide origin by threshold ----------
cat("[fig2] Guide origin...\n")
p_list <- map(lane_names, function(l) {
  df <- lane_data[[l]]$plot2 %>%
    mutate(
      category  = factor(category, levels = c("NC", "OEG", "Mixed", "Unassigned")),
      threshold = factor(threshold, levels = 1:10)
    )
  ggplot(df, aes(x = threshold, y = n_cells, fill = category)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = assignment_colors, name = "Guide origin") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = l, x = "UMI threshold", y = "Proportion")
})
p <- wrap_plots(p_list, nrow = lane_grid$nrow, axis_titles = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Guide origin by UMI threshold") &
  theme(legend.position = "bottom")
ggsave(file.path(figs_dir, "fig2_guide_origin.png"), p, width = lane_grid$width, height = lane_grid$height, dpi = 300, bg = "white")
cat("Saved: figs/fig2_guide_origin.png\n")

# fig3: gRNA recovery ----------
cat("[fig3] gRNA recovery...\n")
p_list <- map(lane_names, function(l) {
  df <- lane_data[[l]]$plot3 %>%
    mutate(
      threshold = factor(threshold, levels = 1:10),
      category = factor(category, levels = c("Total", "NC", "OEG")),
      unassigned = 100 - recovery
    )
  ggplot(df, aes(x = threshold, y = unassigned, color = category, group = category)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(
      values = c(Total = pal[1], NC = pal[2], OEG = pal[3]),
      name = "Guide set"
    ) +
    scale_y_continuous(
      limits = c(0, NA), expand = expansion(mult = c(0, 0.05)),
      labels = function(x) glue("{x}%")
    ) +
    labs(title = l, x = "UMI threshold", y = "% guides unassigned")
})
p <- wrap_plots(p_list, nrow = lane_grid$nrow, axis_titles = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "gRNA recovery by UMI threshold") &
  theme(legend.position = "bottom")
ggsave(file.path(figs_dir, "fig3_gRNA_recovery.png"), p, width = lane_grid$width, height = lane_grid$height, dpi = 300, bg = "white")
cat("Saved: figs/fig3_gRNA_recovery.png\n")

# tbl1: guide origin summary ----------
guide_summary <- map(lane_names, function(l) {
  df <- lane_data[[l]]$plot2
  total_n <- sum(df$n_cells[df$threshold == 1L])
  df %>%
    pivot_wider(
      names_from = category, values_from = n_cells,
      names_glue = "{category}_n_cells"
    ) %>%
    mutate(
      lane = l,
      total_cells = total_n,
      NC_pct = round(NC_n_cells / total_n * 100, 1),
      OEG_pct = round(OEG_n_cells / total_n * 100, 1),
      Mixed_pct = round(Mixed_n_cells / total_n * 100, 1),
      Unassigned_pct = round(Unassigned_n_cells / total_n * 100, 1),
      Mixed_Unassigned_n_cells = Mixed_n_cells + Unassigned_n_cells,
      Mixed_Unassigned_pct = round(Mixed_Unassigned_n_cells / total_n * 100, 1)
    ) %>%
    select(
      lane, threshold, total_cells,
      NC_n_cells, OEG_n_cells, Mixed_n_cells, Unassigned_n_cells,
      NC_pct, OEG_pct, Mixed_pct, Unassigned_pct,
      Mixed_Unassigned_n_cells, Mixed_Unassigned_pct
    )
}) %>% list_rbind()
write_csv(guide_summary, file.path(tables_dir, "tbl1_guide_origin_summary.csv"))
cat("Saved: tables/tbl1_guide_origin_summary.csv\n")

cat("\n=== STEP 1 PLOT.R COMPLETE ===\n")
cat("Review fig2/fig3 to confirm config$qc$umi_threshold before running step2.\n")

cat("\n---\n")
cat("STOP: don't run step2 yet.\n")
cat("Open figs/fig2_guide_origin.png and figs/fig3_gRNA_recovery.png. fig2 shows, for each lane,\n")
cat("what fraction of cells are classified NC / OEG / Mixed / Unassigned as the minimum UMI count\n")
cat("needed to call a guide 'detected' is swept from 1 to 10. fig3 shows what % of your NC/OEG\n")
cat("guides fail to be recovered (detected in at least one cell) at each threshold. Pick the\n")
cat("smallest UMI threshold where the guide_origin proportions in fig2 have stabilized (Mixed/\n")
cat("Unassigned stop shrinking as you raise the threshold further) while guide recovery in fig3\n")
cat("is still high. Current guess in config.yaml (qc.umi_threshold):", config$qc$umi_threshold, "\n")
cat("Tell Claude the confirmed integer threshold so step2_assign can be written with it.\n")
