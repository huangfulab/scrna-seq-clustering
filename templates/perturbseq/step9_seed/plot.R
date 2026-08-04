# conda run -n <conda.env_name> Rscript step9_seed/scripts/plot.R
source(here::here("init.R"))

stopifnot(!is.null(config$final_clustering$final_resolution))
SEEDS      <- unlist(config$final_clustering$seeds)
RESOLUTION <- config$final_clustering$final_resolution

# markers / cell_type_map / ct_levels come from init.R

# ---- Subdir creation ----------
fig_dirs <- setNames(
  map_chr(
    c("fig1_umap", "fig2_dotplot"),
    ~ { d <- here("step9_seed", "figs", .x); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  ),
  c("umap", "dotplot")
)

# ---- Per-seed loop (each seed's obj loaded individually) ----------
walk(SEEDS, function(s) {
  cat(glue("\nLoading obj/seed{s}.qs2...\n"))
  obj <- qs_read(here("step9_seed", "obj", glue("seed{s}.qs2")), nthreads = n_cores)
  clvls      <- levels(obj$seurat_clusters)
  n_clusters <- length(clvls)
  fname      <- glue("seed{s}.png")
  cat(glue("[seed={s}] {n_clusters} clusters\n"))

  # fig1: UMAP cluster
  p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
               label = TRUE, label.size = 4, pt.size = 1,
               raster = TRUE, shuffle = TRUE) +
    labs(title = glue("UMAP — seed={s}, res={RESOLUTION}, n_dims={config$clustering$n_dims}, k={config$clustering$k_param}")) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    coord_cartesian(clip = "on")
  ggsave(file.path(fig_dirs["umap"], fname), p, width = 8, height = 7, dpi = 300, bg = "white")

  # fig2: dotplot markers
  markers_present <- markers[markers %in% rownames(obj)]
  p_base <- DotPlot(obj, features = markers_present, group.by = "seurat_clusters", scale = TRUE, dot.scale = 6)
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

  cat(glue("  2 figs saved\n"))
})

cat("\n=== STEP 9 PLOT.R COMPLETE ===\n")
png_files <- list.files(here("step9_seed", "figs"), pattern = "\\.png$", recursive = TRUE)
cat(glue("{length(png_files)} PNGs total\n"))

cat("\n---\n")
cat("STOP: don't run step10 yet.\n")
cat("Compare figs/umap/seed*.png and figs/dotplot/seed*.png across all seeds in\n")
cat("config.final_clustering.seeds. The cluster COUNT and marker signatures should look mostly\n")
cat("consistent across seeds if the clustering is robust — pick the seed whose UMAP/dotplot look\n")
cat("most representative (typical cluster count, clean marker separation, no seed-specific\n")
cat("fragmentation). Report that seed number back to Claude as config.final_clustering.winning_seed.\n")
cat("Once set, step10_cluster_final has no process.R to run — you (or Claude on your behalf) will\n")
cat("manually copy the winning seed's object forward:\n")
cat("  cp step9_seed/obj/seed<winning_seed>.qs2 step10_cluster_final/obj8_cluster_final.qs2\n")
cat("before running step10_cluster_final/scripts/plot.R.\n")
