# conda run -n <conda.env_name> Rscript step10_cell_type/scripts/process.R
source(here::here("init.R"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = n_cores))

# Filled in by the user after inspecting step9_cluster_final/figs (fig1
# clusters + fig4 dotplot) — see config.cell_type.{cluster_celltype_map,
# cluster_ct_order}. Cluster ID (character) -> CellType label. Must cover
# every cluster level exhaustively. This map is bespoke to this run's
# seed/resolution and is distinct from init.R's shared gene-level
# cell_type_map (marker -> category); do not conflate the two.
celltype_map <- unlist(config$cell_type$cluster_celltype_map)
ct_order     <- unlist(config$cell_type$cluster_ct_order)
stopifnot(length(celltype_map) > 0, length(ct_order) > 0)

# Load ----------
cat("Loading obj8_cluster_final...\n")
obj <- qs_read(here("step9_cluster_final", "obj8_cluster_final.qs2"), nthreads = n_cores)
cat("Loaded:", ncol(obj), "cells,", nlevels(obj$seurat_clusters), "clusters\n")

# Assign CellType ----------
obj$CellType <- factor(unname(celltype_map[as.character(obj$seurat_clusters)]),
                       levels = ct_order)
stopifnot(!anyNA(obj$CellType))
Idents(obj) <- obj$CellType

cat("\nCellType distribution:\n"); print(table(obj$CellType))
cat("\nCluster → CellType crosstab:\n"); print(table(obj$seurat_clusters, obj$CellType))

# Save ----------
cat("\nSaving obj9_cell_type...\n")
out_path <- here("step10_cell_type", "obj9_cell_type.qs2")
qs_save(obj, file = out_path, nthreads = n_cores)
cat(glue("Saved: step10_cell_type/obj9_cell_type.qs2 ({round(file.size(out_path)/1e9, 2)} GB)\n"))

cat("\n=== STEP 10 COMPLETE ===\n")
cat("Run step10_cell_type/plot.R to generate figures.\n")
