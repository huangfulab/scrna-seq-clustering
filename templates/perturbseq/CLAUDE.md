# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this
repository.

This run was scaffolded by the scrna-seq-clustering skill, profile: `perturbseq`.

## Running scripts

```bash
conda activate <config.yaml conda.env_name>
cd <this directory>
Rscript stepN_xxx/scripts/process.R   # data processing (a few steps have no process.R)
Rscript stepN_xxx/scripts/plot.R      # figures & tables
```

Or via SLURM, where a step has a `slurm.sh`:
```bash
sbatch stepN_xxx/scripts/slurm.sh
```

## config.yaml — single source of truth

Every tunable number/path for this pipeline lives in `config.yaml`, not hardcoded in any
step script. Do not hand-edit constants inside `stepN_xxx/scripts/*.R` files — edit
`config.yaml` instead and re-run. `config.yaml` itself has no comments (it's read/written by
Claude Code, not meant for hand-editing) — see `references/conventions.md` in the skill that
scaffolded this run for what each key controls and which step reads it.

## init.R — what it provides

`source(here::here("init.R"))` at the top of every script. Do not repeat any of the
following in a step script:

- Packages, plotting theme, palette (`pal <- brewer.pal(9, "Set1")`)
- `config` — the parsed config.yaml (`config <- yaml::read_yaml(here("config.yaml"))`)
- `n_cores`, `lane_names` (from config.yaml)
- `markers`, `cell_type_map`, `ct_levels` (from config.yaml's `markers.*` block)
- `nc_set`, `oeg_set`, `oeg_gene_map` (perturbseq only — guide-capture CSVs)

## Pipeline design — one step at a time

This run is scaffolded ONE step folder at a time, not all at once (via this skill's
`scripts/scaffold.R`). Most steps end their `plot.R` with a `STOP:` message asking you to
inspect a figure and report a value back before the next step's `config.yaml` fields are
filled in and its folder is scaffolded. Do not skip ahead or hand-write a later step's
scripts yourself — regenerate them from the templates once the checkpoint is confirmed.

## Seurat pipeline constraints

- **Never call `JoinLayers()`** unless you have a specific reason to and understand why —
  the original validated pipeline never needed it.
- **Never copy clusters via `obj$seurat_clusters <- x` expecting `Idents()` to update** —
  `$<-` does NOT update `Idents()` in Seurat v5. Always call `Idents(obj) <- x` explicitly
  if downstream code (e.g. `DimPlot`) relies on `Idents()`.
- **RNA assay rownames are gene symbols**, not Ensembl IDs — match genes on symbol,
  not on any `*_id` column in your guide CSVs.

## Plotting

- Use Seurat plotting functions first (`VlnPlot`, `FeaturePlot`, `DimPlot`, etc.) with the
  Seurat default palette.
- Fall back to ggplot2 + `pal` only when Seurat plotting is not feasible.
- Always pass `bg = "white"` to every `ggsave()` call.

## plot.R conventions

Each step's `stepN_xxx/scripts/plot.R` is the sole entry point for figures and tables.
Do not duplicate any figure or table output in `process.R` if `plot.R` already produces it.

**Output structure**:
- `stepN_xxx/figs/figN_<name>.png` — PNG figures
- `stepN_xxx/tables/tblN_<name>.csv` — CSV tables
