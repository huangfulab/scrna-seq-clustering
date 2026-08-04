# conda run -n <conda.env_name> Rscript step9_cluster_final/scripts/plot.R
#
# No process.R for this step — obj8_cluster_final.qs2 arrives via a manual:
#   cp step8_seed/obj/seed{WINNING_SEED}.qs2 step9_cluster_final/obj8_cluster_final.qs2
# Do NOT write a process.R that reassigns obj$seurat_clusters from another object:
# in Seurat v5 that only sets metadata and does NOT update Idents(), which silently
# breaks every cluster-grouped plot below. The winning seed's file already has
# Idents(tmp) <- tmp$seurat_clusters baked in from step8's worker — the cp preserves it.
source(here::here("init.R"))

RESOLUTION <- config$final_clustering$final_resolution  # for the title label only
out_dir    <- here("step9_cluster_final", "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# markers / cell_type_map / ct_levels come from init.R

# NOTE: renumbering vs. the original pipeline's step10_cluster_final — that
# step had 9 figs. fig4 (guide_origin UMAP, OEG subsampled to NC count) and
# fig8 (guide_origin composition) are both guide-specific and dropped
# entirely — no guide capture in the plain profile. Everything after fig3 is
# shifted down by one (fig5->4, fig6->5), and everything after the dropped
# fig8 is shifted down by two (fig7->6, fig9->7).

# =============================================================================
# LOAD
# =============================================================================
cat("Loading obj8_cluster_final...\n")
obj <- qs_read(here("step9_cluster_final", "obj8_cluster_final.qs2"), nthreads = n_cores)
n_clusters <- nlevels(Idents(obj))
cat("Loaded:", ncol(obj), "cells,", n_clusters, "clusters\n")


# Extract colors matching DimPlot output for a given grouping variable
extract_cols <- function(obj, group.by) {
  vals <- obj@meta.data[[group.by]]
  lvls <- if (is.factor(vals)) levels(vals) else sort(unique(as.character(vals)))
  # When variable is not a factor, DimPlot will order groups alphabetically, so we sort levels accordingly
  p    <- DimPlot(obj, group.by = group.by)
  bld  <- ggplot_build(p)
  df   <- unique(bld$data[[1]][, c("colour", "group")])
  df   <- df[order(df$group), ]
  stopifnot(nrow(df) == length(lvls))
  setNames(df$colour, lvls)
}

lane_cols  <- extract_cols(obj, "orig.ident")
phase_cols <- extract_cols(obj, "Phase")

# =============================================================================
# fig1: UMAP clusters
# =============================================================================
cat("[fig1] UMAP clusters...\n")
p <- DimPlot(obj, label = TRUE, label.size = 4, pt.size = 1.2,
             raster = TRUE, shuffle = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  coord_cartesian(clip = "on") +
  labs(title = glue("Clusters (res={RESOLUTION})")) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.title   = element_text(hjust = 0.5)
  )
ggsave(file.path(out_dir, "fig1_umap_clusters.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")

# =============================================================================
# fig2: UMAP lane
# =============================================================================
cat("[fig2] UMAP lane...\n")
p <- DimPlot(obj, group.by = "orig.ident", cols = lane_cols, pt.size = 1.2,
             raster = TRUE, shuffle = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  coord_cartesian(clip = "on") +
  labs(title = "Lane") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
ggsave(file.path(out_dir, "fig2_umap_lane.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")

# =============================================================================
# fig3: UMAP phase
# =============================================================================
cat("[fig3] UMAP phase...\n")
p <- DimPlot(obj, group.by = "Phase", cols = phase_cols, pt.size = 1.2,
             raster = TRUE, shuffle = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  coord_cartesian(clip = "on") +
  labs(title = "Cell Cycle Phase") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
ggsave(file.path(out_dir, "fig3_umap_phase.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")

# =============================================================================
# fig4 (was fig5): Dotplot markers (faceted by cell type)
# =============================================================================
cat("[fig4] Dotplot markers...\n")
markers_present <- markers[markers %in% rownames(obj)]
p_base <- DotPlot(obj, features = markers_present, group.by = "seurat_clusters",
                  scale = TRUE, dot.scale = 6)
dot_df <- p_base$data
dot_df$features.plot <- factor(dot_df$features.plot, levels = markers_present)
dot_df$strip_labels  <- factor(cell_type_map[as.character(dot_df$features.plot)],
                               levels = ct_levels)

p <- ggplot(dot_df, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_distiller(palette = "RdYlBu") +
  scale_size(range = c(0, 6)) +
  facet_grid(strip_labels ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Cluster", y = NULL,
       color = "Avg exp\n(scaled)", size = "% expressed") +
  theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(out_dir, "fig4_dotplot_markers.png"), p,
       width  = max(8, n_clusters * 0.55 + 3),
       height = max(10, length(ct_levels) * 0.38 + 4),
       dpi = 300, bg = "white")

# =============================================================================
# fig5 (was fig6): Violin QC — 4 metrics (MOI dropped, no guide capture)
# =============================================================================
cat("[fig5] Violin QC...\n")
qc_vars <- c("nCount_RNA", "nFeature_RNA", "pct.mt", "pct.rb")
vln_df <- FetchData(obj, vars = c(qc_vars, "seurat_clusters")) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = levels(obj$seurat_clusters)))
vln_long <- vln_df %>%
  pivot_longer(-seurat_clusters, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = qc_vars))

p <- ggplot(vln_long, aes(x = seurat_clusters, y = value, fill = seurat_clusters)) +
  geom_violin(scale = "width", show.legend = FALSE) +
  facet_wrap(~ metric, ncol = 3, nrow = 2, scales = "free") +
  labs(x = "Cluster", y = NULL)
ggsave(file.path(out_dir, "fig5_violin_qc.png"), p,
       width = 18, height = 10, dpi = 300, bg = "white")

# Helper: stacked bar composition (% cells per cluster by a grouping variable)
composition_bar <- function(obj, group_var, fill_cols, legend_name) {
  df <- FetchData(obj, vars = c("seurat_clusters", group_var)) %>%
    group_by(seurat_clusters, .data[[group_var]]) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(pct = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  ggplot(df, aes(x = seurat_clusters, y = pct, fill = .data[[group_var]])) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fill_cols, name = legend_name) +
    labs(x = "Cluster", y = "% cells")
}

bar_width <- max(8, n_clusters * 0.6 + 3)

# =============================================================================
# fig6 (was fig7): Lane composition
# =============================================================================
cat("[fig6] Lane composition...\n")
p <- composition_bar(obj, "orig.ident", lane_cols, "Lane")
ggsave(file.path(out_dir, "fig6_lane_composition.png"), p,
       width = bar_width, height = 5, dpi = 300, bg = "white")

# =============================================================================
# fig7 (was fig9): Phase composition
# =============================================================================
cat("[fig7] Phase composition...\n")
p <- composition_bar(obj, "Phase", phase_cols, "Phase")
ggsave(file.path(out_dir, "fig7_phase_composition.png"), p,
       width = bar_width, height = 5, dpi = 300, bg = "white")

# =============================================================================
# VERIFICATION
# =============================================================================
cat("\n=== PLOT.R COMPLETE ===\n")
png_files <- sort(list.files(out_dir, pattern = "^fig.*\\.png$"))
cat(glue("{length(png_files)} PNG(s) saved:\n"))
cat(paste0("  ", png_files, collapse = "\n"), "\n")

cat("\n---\n")
cat("STOP: don't run step10 yet.\n")
cat("This is the FINAL cluster set for the whole pipeline. Open fig1_umap_clusters.png (the\n")
cat("cluster layout) alongside fig4_dotplot_markers.png (which marker genes are high in each\n")
cat("cluster). For every cluster number, decide what cell type it represents based on which\n")
cat("marker rows light up (high avg.exp.scaled + high % expressed) for that cluster's column.\n")
cat("Also sanity-check fig2_umap_lane.png / fig6_lane_composition.png (clusters shouldn't be\n")
cat("driven by one lane) and fig3_umap_phase.png / fig7_phase_composition.png / fig5_violin_qc.png\n")
cat("(cell-cycle phase and QC metrics should look reasonable per cluster, not flag a hidden\n")
cat("low-quality cluster). Tell Claude the cluster number -> cell type mapping\n")
cat("(config.cell_type.cluster_celltype_map) and, if you want a specific display order for cell\n")
cat("types in later plots, the order (config.cell_type.cluster_ct_order). Once you update\n")
cat("config.yaml, tell Claude to write step10_cell_type.\n")
