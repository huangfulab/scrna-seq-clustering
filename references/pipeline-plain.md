# Plain scRNA-seq profile (no CRISPR/guide capture) — 10 steps

Same underlying pipeline as `references/pipeline-perturbseq.md`, with every guide-assignment/MOI/UMI-threshold step removed entirely (not branched-in — a plain-profile run never has a `guide_origin`, `MOI`, `guides`, or `target_genes` column at all). Templates live in `templates/plain/stepN_xxx/`. Use this doc for the plain-language checkpoint explanations; `SKILL.md`'s table is the quick-reference version.

| plain step | corresponds to (perturbseq profile) |
|---|---|
| step1_load_qc | step1_load + step2_assign's QC-metric half (no guide logic at all) |
| step2_doublet | step3_doublet |
| step3_filter | step4_filter |
| step4_PCA | step5_PCA |
| step5_batch_effect | step6_batch_effect |
| step6_rm_badcl | step7_rm_badcl |
| step7_res | step8_res |
| step8_seed | step9_seed |
| step9_cluster_final | step10_cluster_final |
| step10_cell_type | step11_cell_type |

## step1_load_qc — load the data, compute QC metrics

Reads each lane's CellRanger `.h5` (Gene Expression only — no CRISPR Guide Capture assay to load), creates a Seurat object per lane, computes `pct.mt` (% mitochondrial reads — a proxy for cell stress/death) and `pct.rb` (% ribosomal reads). There's no UMI-threshold sweep in this profile, since there's no guide capture to threshold.

**Checkpoint (fig1 — QC violin, re-run after confirming):** UMI count, gene count, %mt violin plots per lane, with dashed reference lines at the current `config.qc.*` guesses. Explain: cells to the left of the UMI/gene-count line are probably empty droplets or fragments, not real cells; cells above the %mt line are probably dying/dead. Ask the user to eyeball whether the lines look reasonable, or move them. → `config.qc.min_ncount`, `config.qc.min_nfeat`, `config.qc.max_pct_mt`. **Once confirmed, re-run plot.R so fig1 reflects the real decision before scaffolding step2.**

## step2_doublet — flag doublets

Runs `scDblFinder` per lane (only on cells above `doublet.min_counts` UMIs) to flag likely doublets (two cells captured together as one). Nothing new to ask — still a hard stop so the user can confirm the doublet rate looks sane before moving on.

## step3_filter — apply everything, merge lanes

Purely mechanical: applies the confirmed QC cutoffs and drops non-singlets, then merges all lanes into one Seurat object. Nothing new to ask — still a hard stop so the user can sanity-check the retention rate before moving on.

## step4_PCA — normalize, find variable genes, score cell cycle, PCA

Identical to the perturbseq profile's step5_PCA: normalization → 3000 highly-variable genes → cell-cycle scoring → scale (regressing `%mt`) → PCA (50 components). No batch-correction — this pipeline checks for batch effects visually at the next step instead.

**Checkpoint (fig1 — elbow plot):** each point is one principal component; y-axis is variance explained. Look for where the curve flattens ("the elbow") — use enough PCs to capture real structure without clustering on noise. → `config.clustering.n_dims`.

## Single-lane shortcut

If `config.lanes` has exactly one entry, there's no second lane to compare against, so step5_batch_effect and step6_rm_badcl are skipped as a pair — scaffold step7_res directly after step4_PCA is confirmed. step7_res's `process.R` reads its input object conditionally (from step4_PCA when single-lane, from step6_rm_badcl when multi-lane) so nothing else about it changes. If a bad cluster still shows up in that single clustering pass, there's no scripted removal step for it in this profile — see `SKILL.md`'s "if something doesn't fit the template" section. In this case step7_res becomes the first dotplot-producing step, so the marker-panel confirmation (see step5_batch_effect below) belongs right before scaffolding step7_res instead.

## step5_batch_effect — cluster across a resolution sweep, check for batch effects

**Before scaffolding this step**, confirm the marker panel with the user — this is the first step whose `plot.R` produces a dotplot (`fig2_dotplot`) from `config.markers.panel`/`cell_type_map`/`ct_levels`. Show them the current panel and category grouping and get explicit sign-off (or changes) before proceeding; see `SKILL.md` section 6 for why this is a dedicated checkpoint rather than folded into the general config interview.

`FindNeighbors` → `RunUMAP` → `FindClusters` at 6 resolutions using the confirmed `n_dims`. Produces, per resolution: a UMAP colored by cluster, a marker-gene dotplot, QC violins per cluster, and a lane-composition bar chart.

**Checkpoint (all four fig sets, across resolutions):** low resolution merges similar cell states into fewer, bigger clusters; high resolution splits into more, smaller ones — some of which may be noise rather than real biology. Ask the user to pick the lowest resolution where one cluster clearly stands out as low-quality (globally low gene count, or concentrated in one lane rather than spread across lanes — a technical artifact, not real biology) so it can be removed. → `config.bad_cluster_removal.resolution`, `config.bad_cluster_removal.cluster_rm`.

## step6_rm_badcl — drop the flagged cluster

Purely mechanical: drops every cell in `bad_cluster_removal.cluster_rm` at `bad_cluster_removal.resolution`. Nothing new to ask — still a hard stop so the user can confirm the removal count looks sane.

## step7_res — recluster from scratch on the trimmed data

Full recompute (HVGs → scale → PCA → neighbors → UMAP → clusters) on the post-removal object, across the same 6 resolutions.

**Checkpoint (resolution-sweep figs):** same resolution tradeoff as step5_batch_effect, but now committing to the resolution for the final analysis. Look at whether the marker-gene dotplot's clusters map onto biologically distinct, interpretable groups at each resolution. → `config.final_clustering.final_resolution`.

## step8_seed — check that the clustering is stable across random seeds

Reruns `FindNeighbors → RunUMAP → FindClusters` at the committed resolution across 8 seeds, saving each independently.

**Checkpoint (seed-sweep figs):** if clusters look basically the same shape/count across seeds, the result is trustworthy; if they vary a lot, the resolution or PC count may be poorly chosen. Pick whichever seed gives the cleanest, most stable-looking structure. → `config.final_clustering.winning_seed`.

## step9_cluster_final — no new computation, just the winning seed's figures

No `process.R` — the object arrives via copying the winning seed's file from step8's output (this sidesteps Seurat v5's `Idents()` update gotcha: never reassign `obj$seurat_clusters <- x` expecting `Idents()` to follow — always copy the file that already has `Idents()` baked in). `plot.R` regenerates the figure suite (UMAP by cluster/lane/phase, marker dotplot, QC violins, composition bars — no guide-origin panel in this profile).

**Checkpoint:** show the user the final cluster UMAP and marker dotplot together. The value collected here (which cluster maps to which cell type) feeds step10 — each row block in the dotplot is a marker-category (from `markers.cell_type_map`); columns lighting up strongly for a category's markers suggest that cluster is that cell type. → `config.cell_type.cluster_celltype_map`.

## step10_cell_type — apply the cluster → cell-type map

Applies `cell_type.cluster_celltype_map` (must cover every cluster, exhaustively), regenerates the figure suite grouped by `CellType`. End of pipeline — no new value needed, just confirm the final labeled UMAP looks right.
