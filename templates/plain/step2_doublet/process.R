source(here::here("init.R"))

library(scDblFinder)
library(BiocParallel)

n_total_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_total_cores)) stop("SLURM_CPUS_PER_TASK not set — run via sbatch")
n_lanes_parallel <- 5L
n_cores_per_lane <- n_total_cores %/% n_lanes_parallel   # 8 when total=40

# ---- Directories ----------
dir.create(here("step2_doublet", "obj1_doublet"), recursive = TRUE, showWarnings = FALSE)

# ---- Per-lane doublet detection (5 lanes parallel) ----------
cat("Running scDblFinder on", length(lane_names), "lanes",
    "(", n_lanes_parallel, "parallel,", n_cores_per_lane, "threads/lane)...\n")

plan(multicore, workers = n_lanes_parallel)

lane_results <- future_lapply(lane_names, function(lane) {
  obj <- qs_read(here("step1_load_qc", "obj0_load_qc", glue("{lane}.qs2")), nthreads = n_cores_per_lane)

  obj_sub <- subset(obj, subset = nCount_RNA > config$doublet$min_counts)

  sce <- scDblFinder(
    GetAssayData(obj_sub, assay = "RNA", layer = "counts"),
    BPPARAM   = SerialParam(),
    clusters  = FALSE,
    dims      = config$doublet$dims,
    nfeatures = config$doublet$nfeatures,
    verbose   = FALSE
  )

  # Transfer to full object; cells not in subset get NA
  idx <- match(rownames(obj@meta.data), colnames(sce))
  obj$scDblFinder.score <- sce$scDblFinder.score[idx]
  obj$scDblFinder.class <- as.character(sce$scDblFinder.class[idx])

  qs_save(obj, file = here("step2_doublet", "obj1_doublet", glue("{lane}.qs2")), nthreads = n_cores_per_lane)

  data.frame(
    lane         = lane,
    n_total      = ncol(obj),
    n_eligible   = sum(!is.na(obj$scDblFinder.class)),
    n_singlet    = sum(obj$scDblFinder.class == "singlet", na.rm = TRUE),
    n_doublet    = sum(obj$scDblFinder.class == "doublet", na.rm = TRUE),
    n_ineligible = sum(is.na(obj$scDblFinder.class)),
    stringsAsFactors = FALSE
  )
}, future.seed = TRUE)

plan(sequential)

# ---- Print summary ----------
bind_rows(lane_results) %>%
  mutate(doublet_rate = round(n_doublet / n_eligible * 100, 2)) %>%
  print()

cat("\n=== STEP 2 COMPLETE ===\n")
cat("Run step2_doublet/plot.R to generate figures and tables/tbl1_doublet_summary.csv\n")
