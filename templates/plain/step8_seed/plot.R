# conda run -n <conda.env_name> Rscript step8_seed/scripts/plot.R
source(here::here("init.R"))

SEEDS      <- unlist(config$final_clustering$seeds)
N_DIMS     <- config$clustering$n_dims
K_PARAM    <- config$clustering$k_param
RESOLUTION <- config$final_clustering$final_resolution
stopifnot(!is.null(RESOLUTION))

# markers / cell_type_map / ct_levels come from init.R

# ---- Subdir creation ----------
fig_dirs <- setNames(
  map_chr(
    c("fig1_umap", "fig2_dotplot"),
    ~ { d <- here("step8_seed", "figs", .x); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  ),
  c("umap", "dotplot")
)

# ---- Per-seed loop (each seed's obj loaded individually) ----------
walk(SEEDS, function(s) {
  cat(glue("\nLoading obj/seed{s}.qs2...\n"))
  obj <- qs_read(here("step8_seed", "obj", glue("seed{s}.qs2")), nthreads = n_cores)
  clvls      <- levels(obj$seurat_clusters)
  n_clusters <- length(clvls)
  fname      <- glue("seed{s}.png")
  cat(glue("[seed={s}] {n_clusters} clusters\n"))

  # fig1: UMAP cluster
  p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
               label = TRUE, label.size = 4, pt.size = 1,
               raster = TRUE, shuffle = TRUE) +
    labs(title = glue("UMAP — seed={s}, res={RESOLUTION}, n_dims={N_DIMS}, k={K_PARAM}")) +
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

cat("\n=== STEP 8 PLOT.R COMPLETE ===\n")
png_files <- list.files(here("step8_seed", "figs"), pattern = "\\.png$", recursive = TRUE)
cat(glue("{length(png_files)} PNGs total\n"))

cat("\n---\n")
cat("STOP: don't run step9 yet.\n")
cat("This step reran UMAP + clustering at the SAME resolution (", RESOLUTION, ") but with",
    length(SEEDS), "different\n")
cat("random seeds, to check how sensitive the result is to the random seed alone. Open\n")
cat("figs/fig1_umap/seedX.png for each seed and compare: the cluster SHAPES/COUNT should look\n")
cat("essentially the same across seeds (UMAP layout/orientation can differ — that's expected and\n")
cat("fine — but the number and biological identity of clusters, visible in figs/fig2_dotplot/\n")
cat("seedX.png, should be stable). Pick the seed that gives the most typical/representative,\n")
cat("cleanly-separated result and tell Claude which one (config.final_clustering.winning_seed).\n")
cat("Once you update config.yaml, tell Claude to write step9_cluster_final — the winning seed's\n")
cat("object (step8_seed/obj/seed<N>.qs2) will need to be copied to\n")
cat("step9_cluster_final/obj8_cluster_final.qs2 before that step's plot.R can run.\n")
