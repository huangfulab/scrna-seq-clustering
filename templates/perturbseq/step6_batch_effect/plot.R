# conda run -n <conda.env_name> Rscript step6_batch_effect/scripts/plot.R
source(here::here("init.R"))

stopifnot(!is.null(config$clustering$n_dims))
RESOLUTIONS <- unlist(config$clustering$resolutions)

# markers / cell_type_map / ct_levels come from init.R

# ---- Helpers ----------
composition_bar <- function(obj, group_var, cluster_col) {
  clvls <- levels(obj@meta.data[[cluster_col]])
  FetchData(obj, vars = c(cluster_col, group_var)) %>%
    rename_with(~ "cluster", .cols = 1) %>%
    mutate(cluster = factor(cluster, levels = clvls)) %>%
    group_by(cluster, .data[[group_var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(cluster) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    ggplot(aes(x = cluster, y = pct, fill = .data[[group_var]])) +
    geom_bar(stat = "identity") +
    labs(x = "Cluster", y = "% cells", fill = group_var) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- Subdir creation ----------
tbl_dirs <- setNames(
  map_chr("tbl1_cluster_summary",
          ~ { d <- here("step6_batch_effect", "tables", .x); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }),
  "summary"
)

fig_dirs <- setNames(
  map_chr(
    c("fig1_umap", "fig2_dotplot", "fig3_violin_qc", "fig4_lane_comp"),
    ~ { d <- here("step6_batch_effect", "figs", .x); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  ),
  c("umap", "dotplot", "violin", "lane")
)

# ---- Load ----------
cat("Loading obj5_UMAP...\n")
obj <- qs_read(here("step6_batch_effect", "obj5_UMAP.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj), "cells\n")

# fig5: UMAP by cell-cycle phase (embedding is resolution-independent, so one plot suffices)
p <- DimPlot(obj, reduction = "umap", group.by = "Phase",
             pt.size = 1, raster = TRUE, shuffle = TRUE) +
  labs(title = "UMAP — cell-cycle phase") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  coord_cartesian(clip = "on")
ggsave(here("step6_batch_effect", "figs", "fig5_umap_phase.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")
cat("Saved: figs/fig5_umap_phase.png\n")

# ---- Per-resolution loop ----------
walk(RESOLUTIONS, function(res) {
  col        <- glue("RNA_snn_res.{res}")
  clvls      <- levels(obj@meta.data[[col]])
  n_clusters <- length(clvls)
  bar_width  <- max(8, n_clusters * 0.6 + 3)
  fname      <- glue("res{res}.png")
  cat(glue("\n[res={res}] {n_clusters} clusters\n"))

  # fig1: UMAP cluster
  p <- DimPlot(obj, reduction = "umap", group.by = col,
               label = TRUE, label.size = 4, pt.size = 1,
               raster = TRUE, shuffle = TRUE) +
    labs(title = glue("UMAP — res={res}, n_dims={config$clustering$n_dims}, k={config$clustering$k_param}")) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    coord_cartesian(clip = "on")
  ggsave(file.path(fig_dirs["umap"], fname), p, width = 8, height = 7, dpi = 300, bg = "white")

  # fig2: dotplot markers
  markers_present <- markers[markers %in% rownames(obj)]
  p_base <- DotPlot(obj, features = markers_present, group.by = col, scale = TRUE, dot.scale = 6)
  dot_df <- p_base$data %>%
    mutate(
      features.plot = factor(features.plot, levels = markers_present),
      strip_labels  = factor(cell_type_map[as.character(features.plot)], levels = ct_levels)
    )
  p <- ggplot(dot_df, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_color_distiller(palette = "RdYlBu") +
    scale_size(range = c(0, 6)) +
    facet_grid(strip_labels ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Cluster", y = NULL, color = "Avg exp\n(scaled)", size = "% expressed") +
    theme(strip.text.y = element_text(angle = 0))
  ggsave(file.path(fig_dirs["dotplot"], fname), p,
         width  = max(8, n_clusters * 0.55 + 3),
         height = max(10, length(ct_levels) * 0.38 + 4),
         dpi = 300, bg = "white")

  # fig3: violin QC
  qc_vars  <- c("nCount_RNA", "nFeature_RNA", "pct.mt", "pct.rb", "MOI")
  vln_long <- FetchData(obj, vars = c(col, qc_vars)) %>%
    rename_with(~ "cluster", .cols = 1) %>%
    mutate(cluster = factor(cluster, levels = clvls)) %>%
    pivot_longer(-cluster, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = qc_vars))
  p <- ggplot(vln_long, aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(scale = "width", show.legend = FALSE) +
    facet_wrap(~ metric, ncol = 3, nrow = 2, scales = "free") +
    labs(x = "Cluster", y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dirs["violin"], fname), p,
         width = max(18, bar_width), height = 10, dpi = 300, bg = "white")

  # fig4: lane composition
  p <- composition_bar(obj, "orig.ident", col)
  ggsave(file.path(fig_dirs["lane"], fname), p, width = bar_width, height = 5, dpi = 300, bg = "white")

  # tbl1: cluster summary
  tbl <- FetchData(obj, vars = c(col, "guide_origin", "nCount_RNA", "nFeature_RNA",
                                 "pct.mt", "MOI", "Phase")) %>%
    rename_with(~ "cluster", .cols = 1) %>%
    mutate(cluster = factor(cluster, levels = clvls)) %>%
    group_by(cluster) %>%
    summarise(
      n_cells    = n(),
      pct_OEG    = round(mean(guide_origin == "OEG") * 100, 1),
      pct_NC     = round(mean(guide_origin == "NC")  * 100, 1),
      med_nCount = round(median(nCount_RNA),   0),
      med_nFeat  = round(median(nFeature_RNA), 0),
      med_pct_mt = round(median(pct.mt),       2),
      med_MOI    = round(median(MOI),          1),
      pct_G1     = round(mean(Phase == "G1")  * 100, 1),
      pct_S      = round(mean(Phase == "S")   * 100, 1),
      pct_G2M    = round(mean(Phase == "G2M") * 100, 1),
      .groups = "drop"
    ) %>%
    mutate(pct_cells = round(n_cells / sum(n_cells) * 100, 2)) %>%
    relocate(pct_cells, .after = n_cells)
  write_csv(tbl, file.path(tbl_dirs["summary"], glue("res{res}.csv")))

  cat(glue("  4 figs + 1 tbl saved\n"))
})

cat("\n=== STEP 6 PLOT.R COMPLETE ===\n")
png_files <- list.files(here("step6_batch_effect", "figs"), pattern = "\\.png$", recursive = TRUE)
cat(glue("{length(png_files)} PNGs total\n"))

cat("\n---\n")
cat("STOP: don't run step7 yet.\n")
cat("figs/fig5_umap_phase.png colors the UMAP by cell-cycle phase (G1/S/G2M) — useful as a\n")
cat("sanity check if pca.vars_to_regress includes S.Score/G2M.Score: phases should look mixed\n")
cat("within each cluster rather than forming their own sub-blobs.\n")
cat("For each resolution in figs/fig1_umap/, figs/fig2_dotplot/, figs/fig3_violin_qc/ and\n")
cat("figs/fig4_lane_comp/ (filenames res0.1.png ... res0.6.png), look for a cluster that stands out\n")
cat("as low quality or technical rather than biological — e.g. very low nCount/nFeature in\n")
cat("fig3_violin_qc, near-zero marker expression across the board in fig2_dotplot, or a cluster\n")
cat("that is almost entirely one lane in fig4_lane_comp. Pick the LOWEST resolution where such a\n")
cat("cluster is clearly and stably isolated, and note its cluster ID (the label shown on the UMAP\n")
cat("in fig1_umap, e.g. \"3\"). Tell Claude both values so config.yaml can be updated:\n")
cat("  bad_cluster_removal.resolution  (e.g. 0.1)\n")
cat("  bad_cluster_removal.cluster_rm  (e.g. \"3\")\n")
cat("If nothing looks clearly bad, say so — it is valid to skip removal; ask Claude how to proceed.\n")
