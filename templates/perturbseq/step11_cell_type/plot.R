# conda run -n <conda.env_name> Rscript step11_cell_type/scripts/plot.R
source(here::here("init.R"))

out_dir <- here("step11_cell_type", "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# markers / cell_type_map / ct_levels come from init.R

# =============================================================================
# LOAD
# =============================================================================
cat("Loading obj9_cell_type...\n")
obj <- qs_read(here("step11_cell_type", "obj9_cell_type.qs2"), nthreads = n_cores)
n_celltypes <- nlevels(Idents(obj))
cat("Loaded:", ncol(obj), "cells,", n_celltypes, "cell types\n")


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

celltype_cols     <- extract_cols(obj, "CellType")
lane_cols         <- extract_cols(obj, "orig.ident")
phase_cols        <- extract_cols(obj, "Phase")
guide_origin_cols <- extract_cols(obj, "guide_origin")

# =============================================================================
# fig1: UMAP cell types
# =============================================================================
cat("[fig1] UMAP cell types...\n")
p <- DimPlot(obj, group.by = "CellType", cols = celltype_cols,
             label = TRUE, label.size = 4, pt.size = 1.2,
             raster = TRUE, shuffle = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  coord_cartesian(clip = "on") +
  theme(
    text         = element_text(family = "Arial"),
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank()
  )
ggsave(file.path(out_dir, "fig1_umap_celltype.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")

# =============================================================================
# fig2: Dotplot markers (faceted by cell type, grouped by CellType)
# =============================================================================
cat("[fig2] Dotplot markers...\n")
markers_present <- markers[markers %in% rownames(obj)]
p_base <- DotPlot(obj, features = markers_present, group.by = "CellType",
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
  labs(x = "Cell type", y = NULL,
       color = "Avg exp\n(scaled)", size = "% expressed") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text.y = element_text(angle = 0))
ggsave(file.path(out_dir, "fig2_dotplot_markers.png"), p,
       width  = max(8, n_celltypes * 0.55 + 5),
       height = max(10, length(ct_levels) * 0.38 + 4),
       dpi = 300, bg = "white")

# =============================================================================
# fig3: Violin QC — 5 metrics, 3-col × 2-row, grouped by CellType
# =============================================================================
cat("[fig3] Violin QC...\n")
qc_vars <- c("nCount_RNA", "nFeature_RNA", "pct.mt", "pct.rb", "MOI")
vln_df <- FetchData(obj, vars = c(qc_vars, "CellType")) %>%
  mutate(CellType = factor(CellType, levels = levels(obj$CellType)))
vln_long <- vln_df %>%
  pivot_longer(-CellType, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = qc_vars))

p <- ggplot(vln_long, aes(x = CellType, y = value, fill = CellType)) +
  geom_violin(scale = "width", show.legend = FALSE) +
  facet_wrap(~ metric, ncol = 3, nrow = 2, scales = "free") +
  labs(x = "Cell type", y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "fig3_violin_qc.png"), p,
       width = 18, height = 10, dpi = 300, bg = "white")

# Helper: stacked bar composition (% cells per cell type by a grouping variable)
composition_bar <- function(obj, group_var, fill_cols, legend_name) {
  df <- FetchData(obj, vars = c("CellType", group_var)) %>%
    group_by(CellType, .data[[group_var]]) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(CellType) %>%
    mutate(pct = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  ggplot(df, aes(x = CellType, y = pct, fill = .data[[group_var]])) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fill_cols, name = legend_name) +
    labs(x = "Cell type", y = "% cells") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

bar_width <- max(8, n_celltypes * 0.9 + 3)

# =============================================================================
# fig4: Lane composition
# =============================================================================
cat("[fig4] Lane composition...\n")
p <- composition_bar(obj, "orig.ident", lane_cols, "Lane")
ggsave(file.path(out_dir, "fig4_lane_composition.png"), p,
       width = bar_width, height = 5, dpi = 300, bg = "white")

# =============================================================================
# fig5: Guide origin composition
# =============================================================================
cat("[fig5] Guide origin composition...\n")
p <- composition_bar(obj, "guide_origin", guide_origin_cols, "Guide origin")
ggsave(file.path(out_dir, "fig5_guide_origin_composition.png"), p,
       width = bar_width, height = 5, dpi = 300, bg = "white")

# =============================================================================
# fig6: Phase composition
# =============================================================================
cat("[fig6] Phase composition...\n")
p <- composition_bar(obj, "Phase", phase_cols, "Phase")
ggsave(file.path(out_dir, "fig6_phase_composition.png"), p,
       width = bar_width, height = 5, dpi = 300, bg = "white")

# =============================================================================
# VERIFICATION
# =============================================================================
cat("\n=== PLOT.R COMPLETE ===\n")
png_files <- sort(list.files(out_dir, pattern = "^fig.*\\.png$"))
cat(glue("{length(png_files)} PNG(s) saved:\n"))
cat(paste0("  ", png_files, collapse = "\n"), "\n")

cat("\n---\n")
cat("STOP: pipeline complete.\n")
cat("Open figs/fig1_umap_celltype.png (final cell-type labels on the UMAP) and\n")
cat("figs/fig2_dotplot_markers.png (marker expression per cell type) as a final sanity check that\n")
cat("each cell type's dominant marker category matches its assigned label. Also skim\n")
cat("figs/fig3_violin_qc.png, figs/fig4_lane_composition.png and figs/fig6_phase_composition.png\n")
cat("for anything that looks like a lingering batch or QC artifact. There is nothing further to\n")
cat("confirm back to Claude — this is the end of the clustering pipeline. obj9_cell_type.qs2 is the\n")
cat("final annotated object for downstream analysis (e.g. 02_on-target/-style differential expression).\n")
