# Perturb-seq / CRISPR-screen profile — 11 steps

Ported from a validated run analyzing cardiac-progenitor organoids with an overexpression-gene (OEG) CRISPR screen. Templates live in `templates/perturbseq/stepN_xxx/`. Each step folder has `process.R` (computation, saves a `.qs2` Seurat object) and `plot.R` (all figures/tables — never duplicated in process.R), and (except step10/step11, which run interactively) a `slurm.sh` that runs both in sequence.

Use this doc for the plain-language explanation to give the user at each checkpoint — `SKILL.md`'s checkpoint table is the quick-reference version of this.

## step1_load — load the data, sweep UMI thresholds

Reads each lane's CellRanger `.h5` (Gene Expression + CRISPR Guide Capture assays), creates a Seurat object per lane, and — purely as a diagnostic, not a filter — tries every UMI threshold from 1 to 10 to see how many cells would be classified NC (negative control guide), OEG (target guide), Mixed (both), or Unassigned at each threshold.

**Checkpoint (fig2 — guide origin by threshold):** for each lane, this shows what fraction of cells fall into each category as the UMI threshold moves from 1 to 10. Explain to the user: a threshold that's too low lets noise (ambient RNA, PCR crosstalk) masquerade as a real guide; too high throws away real signal from cells that were sequenced shallowly. Look for a threshold where the categories stop changing much as you move it — that plateau is usually the right choice. → `config.qc.umi_threshold`.

## step2_assign — fix the UMI threshold, compute per-cell guide identity + QC

Applies the confirmed `qc.umi_threshold`, extracts each cell's guide identity (semicolon-joined if multiple), MOI (guides per cell), maps guide → target gene, computes `pct.mt` (% mitochondrial reads — a proxy for cell stress/death) and `pct.rb` (% ribosomal reads).

**Checkpoint (fig1 — QC violin, re-run after confirming):** four violin plots per lane — UMI count, gene count, %mt, %ribosomal — with dashed reference lines at the current `config.qc.*` guesses. Explain: cells to the left of the UMI/gene-count line are probably empty droplets or fragments, not real cells; cells above the %mt line are probably dying/dead (their cytoplasmic RNA degraded, leaving a higher relative share of mitochondrial transcripts). Ask the user to eyeball whether the current lines look reasonable for their data's distribution, or to move them. → `config.qc.min_ncount`, `config.qc.min_nfeat`, `config.qc.max_pct_mt`. **Once confirmed, re-run plot.R so fig1 reflects the real decision, then show the user the regenerated fig1 and get explicit confirmation they're satisfied with it — not just the re-run itself — before scaffolding step3.**

## step3_doublet — flag doublets, decide MOI cutoffs

Runs `scDblFinder` per lane (only on cells above `doublet.min_counts` UMIs) to flag cells that are likely two cells captured together (a "doublet") rather than one.

**Checkpoint (fig3 — MOI vs UMI count, re-run after confirming):** doublets are more common at high MOI (a cell with many guides detected is more likely to actually be two+ real cells stuck together). Explain: this figure shows median UMI count as a function of MOI, with a reference line at the current MOI cutoff. A high MOI paired with an unusually high UMI count is a red flag for doublets even beyond what scDblFinder caught. Ask for the min/max MOI to keep. → `config.qc.min_moi`, `config.qc.max_moi`. **Once confirmed, re-run plot.R so fig3 reflects the real decision, then show the user the regenerated fig3 and get explicit confirmation they're satisfied with it — not just the re-run itself — before scaffolding step4.**

## step4_filter — apply everything, merge lanes

Purely mechanical: applies `qc.min_ncount`/`qc.min_nfeat`/`qc.max_pct_mt` (step2), `qc.min_moi`/`qc.max_moi` (step3), drops non-singlets (step3's doublet call), drops cells that aren't cleanly NC-only or OEG-only, then merges all lanes into one Seurat object. Nothing new to ask — still a hard stop so the user can sanity-check the retention rate (how many cells survived) before moving on.

## step5_PCA — normalize, find variable genes, score cell cycle, PCA

Standard Seurat normalization → 3000 highly-variable genes → cell-cycle phase scoring (using Seurat's built-in human gene lists) → scale (regressing out `%mt`) → PCA computing 50 components. No batch-correction (e.g. Harmony) — this pipeline deliberately clusters on raw PCA and checks for batch effects visually at the next step instead of pre-correcting them away.

**Checkpoint (fig1 — elbow plot):** each point is one principal component; y-axis is how much variance it explains. Explain: after some number of PCs, each additional one adds mostly noise, not real biological signal — you want to use enough PCs to capture the real structure but not so many that you're clustering on noise. Look for where the curve visibly flattens ("the elbow"). → `config.clustering.n_dims`.

## Single-lane shortcut

If `config.lanes` has exactly one entry, there's no second lane to compare against, so step6_batch_effect and step7_rm_badcl are skipped as a pair — scaffold step8_res directly after step5_PCA is confirmed. step8_res's `process.R` reads its input object conditionally (from step5_PCA when single-lane, from step7_rm_badcl when multi-lane) so nothing else about it changes. If a bad cluster still shows up in that single clustering pass, there's no scripted removal step for it in this profile — see `SKILL.md`'s "if something doesn't fit the template" section. In this case step8_res becomes the first dotplot-producing step, so the marker-panel confirmation (see step6_batch_effect below) belongs right before scaffolding step8_res instead.

## step6_batch_effect — cluster across a resolution sweep, check for batch effects

**Before scaffolding this step**, confirm the marker panel with the user — this is the first step whose `plot.R` produces a dotplot (`fig2_dotplot`) from `config.markers.cell_type_map`/`ct_levels`. Show them the current panel and category grouping and get explicit sign-off (or changes) before proceeding; see `SKILL.md` section 6 for why this is a dedicated checkpoint rather than folded into the general config interview.

`FindNeighbors` → `RunUMAP` → `FindClusters` at 6 resolutions (0.1 to 0.6) using the confirmed `n_dims`. Produces, per resolution: a UMAP colored by cluster, a marker-gene dotplot, QC violins per cluster, and a lane-composition bar chart.

**Checkpoint (all four fig sets, across resolutions):** explain resolution as a "how finely to split" knob — low resolution merges similar cell states into fewer, bigger clusters; high resolution splits into more, smaller ones, some of which may just be noise or a QC-driven artifact rather than real biology. Ask the user to pick the lowest resolution where one cluster clearly stands out as low-quality (e.g. globally low gene count, or all in one lane instead of spread across lanes — a sign of a batch/technical artifact rather than real biology) so it can be removed. → `config.bad_cluster_removal.resolution`, `config.bad_cluster_removal.cluster_rm`.

## No-bad-cluster shortcut

If nothing looked bad across the resolution sweep — no cluster stands out as low-quality or lane-driven at any of the 6 resolutions — there's nothing to remove. Skip `step7_rm_badcl` and `step8_res` as a pair rather than scaffolding either as a no-op or a redundant recompute: recomputing HVG/scale/PCA/neighbors/UMAP/clusters on the exact same cells `step6_batch_effect` already clustered would only reproduce that same result. Set `config.bad_cluster_removal.no_bad_cluster: true`, leave `cluster_rm: null`, and copy the chosen `bad_cluster_removal.resolution` into `config.final_clustering.final_resolution` — that resolution IS the final one, since there's no separate resweep to reconsider it from. Scaffold `step9_seed` directly after this; its `process.R` reads `step6_batch_effect`'s `obj5_UMAP.qs2` directly instead of `step8_res`'s `obj7_cluster2.qs2` when this flag is set.

## step7_rm_badcl — drop the flagged cluster

Purely mechanical: drops every cell in `bad_cluster_removal.cluster_rm` at `bad_cluster_removal.resolution`. Nothing new to ask — still a hard stop so the user can confirm the removal count looks sane.

## step8_res — recluster from scratch on the trimmed data

Full recompute (HVGs → scale → PCA → neighbors → UMAP → clusters) on the post-removal object, across the same 6 resolutions. Same fig/table shapes as step6.

**Checkpoint (resolution-sweep figs):** same "how finely to split" explanation as step6, but now the user is picking the resolution to actually commit to for the final analysis, not just diagnosing a bad cluster. Encourage looking at whether the marker-gene dotplot's clusters map onto biologically distinct, interpretable groups at each resolution. → `config.final_clustering.final_resolution`.

## step9_seed — check that the clustering is stable across random seeds

UMAP and Louvain clustering both have a random-seed dependency. This step reruns `FindNeighbors → RunUMAP → FindClusters` at the committed `final_resolution` across 8 different seeds, saving each independently.

**Checkpoint (seed-sweep figs):** explain: if the clusters look basically the same shape/count across seeds, the result is trustworthy; if they vary a lot, the resolution or PC count may be poorly chosen. Ask the user to pick whichever seed gives the cleanest, most interpretable UMAP (not necessarily most clusters — most *stable-looking* structure). → `config.final_clustering.winning_seed`.

## step10_cluster_final — no new computation, just the winning seed's figures

There is no `process.R` here by design — the object arrives via copying the winning seed's file from step9's output into this step's folder: `cp step9_seed/obj/seed<winning_seed>.qs2 step10_cluster_final/obj8_cluster_final.qs2` (`<winning_seed>` = `config.final_clustering.winning_seed`, the value just confirmed at step9's checkpoint). Show the user this exact `cp` command per non-negotiable behavior 1 rather than running it silently — it's a real file operation, not a step script. Do this copy, not a script (copying the object programmatically risks Seurat v5's `Idents()` update gotcha: never reassign `obj$seurat_clusters <- x` expecting `Idents()` to follow — always copy the file that already has `Idents()` baked in). `plot.R` regenerates the full figure suite (UMAP by cluster/lane/phase/guide-origin, marker dotplot, QC violins, composition bars) for the committed seed+resolution.

**Checkpoint:** show the user the final cluster UMAP and the marker-gene dotplot together. Nothing new to compute, but the value collected here (which cluster maps to which cell type, from the dotplot) feeds step11 — explain that each row block in the dotplot is a marker-category (from `markers.cell_type_map`), and columns lighting up strongly for a category's markers suggest that cluster is that cell type. → `config.cell_type.cluster_celltype_map`.

## step11_cell_type — apply the cluster → cell-type map

Applies `cell_type.cluster_celltype_map` (must cover every cluster, exhaustively), regenerates the same figure suite grouped by `CellType` instead of raw cluster ID. End of pipeline — no new value needed, just confirm the final labeled UMAP looks right.
