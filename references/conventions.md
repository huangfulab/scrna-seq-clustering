# Conventions

Shared rules that apply across both profiles (`perturbseq` and `plain`). Read this once per new run before scaffolding step1.

## `here()` root-per-run

Each scaffolded run is its own independent `here::here()` root — an empty `.here` file lives at the destination's top level, and every script's first line is `source(here::here("init.R"))`. Do not nest a new run inside another run's root, and do not point a new run's `here()` at another pipeline stage's root. Each pipeline stage should own its own `here()` root rather than sharing one at a project's top level.

If the new run needs to reach a shared resource one level up (e.g. a `gRNA/` directory shared by multiple pipeline stages), use `here("..", "some_dir", ...)` — never a hardcoded absolute path.

## config.yaml

Every step script reads `config <- yaml::read_yaml(here::here("config.yaml"))` once, via `init.R` (so it's available everywhere `source(here::here("init.R"))` runs — no step repeats the read). Full schema, with the exact key path each pipeline step reads it from:

| Key | Meaning | Read by |
|---|---|---|
| `profile` | `"perturbseq"` or `"plain"` — which template set was scaffolded | scaffold.R only |
| `lanes` | Sample/lane IDs, in the order figures should display them | every step |
| `paths.h5_dir`, `paths.h5_pattern` | Where CellRanger `.h5` files live and how they're named (`{lane}` gets substituted) | step1 (load) |
| `paths.nc_gRNA_csv`, `paths.oeg_gRNA_csv` | Negative-control / target guide reference CSVs | perturbseq only, step1-2 |
| `qc.umi_threshold` | Minimum UMI count for a guide to "count" as detected in a cell | perturbseq only, step1→step2 checkpoint |
| `qc.min_ncount`, `qc.min_nfeat`, `qc.max_pct_mt` | Per-cell QC filter thresholds (UMI count, gene count, % mitochondrial) | step2 checkpoint → applied at the filter step |
| `qc.min_moi`, `qc.max_moi` | Multiplicity-of-infection (guides per cell) filter range | perturbseq only, step3 checkpoint → applied at the filter step |
| `doublet.min_counts`, `doublet.dims`, `doublet.nfeatures` | scDblFinder parameters | doublet-detection step |
| `pca.hvg_nfeatures`, `pca.npcs`, `pca.seed`, `pca.vars_to_regress` | Normalization/HVG/PCA parameters. `pca.npcs` is how many PCs get *computed* (generous, fixed) — not the same as `clustering.n_dims` (how many get *used*, the elbow-plot checkpoint) | PCA step |
| `clustering.n_dims`, `clustering.k_param`, `clustering.resolutions`, `clustering.umap_seed`, `clustering.algorithm` | FindNeighbors/RunUMAP/FindClusters parameters | clustering steps |
| `bad_cluster_removal.resolution`, `bad_cluster_removal.cluster_rm` | Which resolution + cluster ID to drop as an outlier | filled at the batch-effect checkpoint |
| `final_clustering.final_resolution`, `final_clustering.seeds`, `final_clustering.winning_seed` | Committed resolution, the seed sweep list, and which seed won | filled at the resolution-sweep and seed-sweep checkpoints |
| `markers.panel`, `markers.cell_type_map`, `markers.ct_levels` | Marker gene panel for dotplots and their category labels — **always study-specific, ask the user, never reuse the cardiac-progenitor example** | clustering/final/cell-type steps |
| `cell_type.cluster_celltype_map`, `cell_type.cluster_ct_order` | Cluster ID → cell-type label, and display order | filled at the final checkpoint, applied in the last step |
| `conda.env_name` | Conda/mamba environment every `slurm.sh` activates | scaffold.R (substituted into slurm.sh) |
| `slurm.mail_user`, `slurm.log_dir`, `slurm.partition`, `slurm.project_root`, `slurm.resources.*` | SLURM header values and per-step resource requests | scaffold.R (substituted into slurm.sh) |

R lists loaded from YAML sequences (`clustering.resolutions`, `final_clustering.seeds`) need `unlist()` to behave as plain numeric vectors — the templates already do this; don't remove it if you're hand-editing a scaffolded script.

## SLURM

Every `slurm.sh` this skill writes must keep the four directives required project-wide (from the global `CLAUDE.md`):

```bash
#SBATCH --output=<slurm.log_dir>/slurm_%j.log
#SBATCH --error=<slurm.log_dir>/slurm_%j.err
#SBATCH --mail-user=<slurm.mail_user>
#SBATCH --mail-type=BEGIN,END,FAIL
```

plus `--cpus-per-task`/`--mem`/`--time`/`--partition` from `slurm.resources.<step>`, a `conda activate <conda.env_name>` line, and `cd <slurm.project_root>`. `scripts/scaffold.R` fills these in from `config.yaml` when it copies a step's `slurm.sh` template — the template itself only has `__PLACEHOLDER__` tokens, since a shell script can't read YAML on its own.

**The skill hands the user the `sbatch` command — it never runs `sbatch` itself.** See `SKILL.md`'s non-negotiable behaviors.

## Conda environment

`templates/environment.yml.template` lists the packages every step needs (Seurat v5, qs2, ggplot2, dplyr, tidyr, purrr, readr, here, RColorBrewer, future, future.apply, Matrix, glue, patchwork, yaml, plus per-step extras: scDblFinder + BiocParallel for doublet detection, ComplexUpset for the disqualification UpSet plot). Give the user the exact `conda create`/`mamba create` command from that file before scaffolding step1 — every generated `slurm.sh` assumes the named environment already exists.

## LOG.md

Once a destination has produced its first real output, suggest starting a `LOG.md` there with a parameters table, an append-only changelog (what changed and why), a SLURM job-history table, and an Output Files / Figures inventory. This skill doesn't write LOG.md automatically — that's the user's (or the assisting Claude session's) lab notebook to maintain interactively as the run actually happens, not something to template in advance.
