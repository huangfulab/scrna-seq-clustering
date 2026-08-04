# conda run -n <conda.env_name> Rscript step3_doublet/scripts/plot.R
source(here::here("init.R"))

MOI_DISPLAY_CAP  <- 40L
MOI_FILTER_LINE  <- config$qc$max_moi   # MAX_MOI cutoff applied in step4
doublet_colors   <- c(singlet = "grey70", doublet = pal[2])

figs_dir   <- here("step3_doublet", "figs")
tables_dir <- here("step3_doublet", "tables")
dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Derive plot data from per-lane objects ----------------------------------
cat("Reading per-lane objects...\n")
all_meta <- map(lane_names, function(l) {
  obj <- qs_read(here("step3_doublet", "obj2_doublet", glue("{l}.qs2")), nthreads = n_cores)
  cat(" ", l, "\n")
  obj@meta.data %>%
    select(MOI, nCount_RNA, scDblFinder.class) %>%
    mutate(lane = l)
}) %>%
  list_rbind() %>%
  mutate(lane = factor(lane, levels = lane_names))

# ---- tbl1: doublet summary --------------------------------------------------
doublet_summary <- all_meta %>%
  group_by(lane) %>%
  summarise(
    n_total      = n(),
    n_eligible   = sum(!is.na(scDblFinder.class)),
    n_singlet    = sum(scDblFinder.class == "singlet",  na.rm = TRUE),
    n_doublet    = sum(scDblFinder.class == "doublet",  na.rm = TRUE),
    n_ineligible = sum(is.na(scDblFinder.class)),
    doublet_rate = round(n_doublet / n_eligible * 100, 2),
    .groups = "drop"
  )
write_csv(doublet_summary, file.path(tables_dir, "tbl1_doublet_summary.csv"))
cat("Saved: tables/tbl1_doublet_summary.csv\n")

plot_meta <- all_meta %>%
  filter(!is.na(scDblFinder.class)) %>%
  mutate(
    scDblFinder.class = factor(scDblFinder.class, levels = c("singlet", "doublet")),
    MOI_plot = pmin(MOI, MOI_DISPLAY_CAP)
  )

# fig1: MOI vs doublet count ----------
p_list <- map(lane_names, function(l) {
  df <- filter(plot_meta, lane == l)
  ggplot(df, aes(x = MOI_plot, fill = scDblFinder.class)) +
    geom_bar(linewidth = 0) +
    scale_fill_manual(values = doublet_colors, name = NULL) +
    scale_x_continuous(
      breaks = seq(0, MOI_DISPLAY_CAP, by = 5),
      labels = c(seq(0, MOI_DISPLAY_CAP - 5, by = 5), glue("≥{MOI_DISPLAY_CAP}"))
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = l, x = "MOI", y = "N cells")
})

p <- wrap_plots(p_list, nrow = 2, axis_titles = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Doublet count by MOI per lane") &
  theme(legend.position = "bottom")

ggsave(file.path(figs_dir, "fig1_MOI_doublet_count.png"), p, width = 20, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig1_MOI_doublet_count.png\n")

# fig2: MOI vs doublet fraction ----------
p_list <- map(lane_names, function(l) {
  df <- filter(plot_meta, lane == l, MOI_plot > 0) %>%
    count(MOI_plot, scDblFinder.class, .drop = FALSE) %>%
    group_by(MOI_plot) %>%
    mutate(frac = n / sum(n)) %>%
    filter(scDblFinder.class == "doublet")
  ggplot(df, aes(x = MOI_plot, y = frac)) +
    geom_col(fill = doublet_colors[["doublet"]], width = 0.8) +
    scale_x_continuous(
      breaks = seq(0, MOI_DISPLAY_CAP, by = 5),
      labels = c(seq(0, MOI_DISPLAY_CAP - 5, by = 5), glue("≥{MOI_DISPLAY_CAP}"))
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = l, x = "MOI", y = "Doublet fraction")
})

p <- wrap_plots(p_list, nrow = 2, axis_titles = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Doublet fraction by MOI per lane")

ggsave(file.path(figs_dir, "fig2_MOI_doublet_percent.png"), p, width = 20, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig2_MOI_doublet_percent.png\n")

# fig3: MOI vs nCount_RNA (all vs singlets) ----------
ncount_colors    <- c("All Cells" = "grey60", "Singlets Only" = pal[2])

summarise_ncount <- function(df) {
  df %>%
    group_by(MOI_plot) %>%
    summarise(
      median_nCount = median(nCount_RNA, na.rm = TRUE),
      q25 = quantile(nCount_RNA, 0.25, na.rm = TRUE),
      q75 = quantile(nCount_RNA, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

p_list <- map(lane_names, function(l) {
  df <- bind_rows(
    filter(plot_meta, lane == l) %>%
      summarise_ncount() %>% mutate(group = "All Cells"),
    filter(plot_meta, lane == l, scDblFinder.class == "singlet") %>%
      summarise_ncount() %>% mutate(group = "Singlets Only")
  ) %>%
    mutate(group = factor(group, levels = c("All Cells", "Singlets Only"))) %>%
    filter(MOI_plot < MOI_DISPLAY_CAP)

  ggplot(df, aes(x = MOI_plot, y = median_nCount, color = group, fill = group)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = ncount_colors, name = NULL) +
    scale_fill_manual(values  = ncount_colors, name = NULL) +
    scale_x_continuous(
      breaks = seq(0, MOI_DISPLAY_CAP, by = 5),
      labels = c(seq(0, MOI_DISPLAY_CAP - 5, by = 5), glue("≥{MOI_DISPLAY_CAP}"))
    ) +
    geom_vline(xintercept = MOI_FILTER_LINE, linetype = "dashed", color = pal[1], linewidth = 0.7) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = l, x = "MOI", y = "nCount_RNA")
})

p <- wrap_plots(p_list, nrow = 2, axis_titles = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(title = "MOI vs nCount per lane") &
  theme(legend.position = "bottom")

ggsave(file.path(figs_dir, "fig3_MOI_nCount.png"), p, width = 20, height = 8, dpi = 300, bg = "white")
cat("Saved: figs/fig3_MOI_nCount.png\n")

cat("\n=== STEP 3 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: don't run step4 yet.\n")
cat("Open figs/fig3_MOI_nCount.png. It shows, for each lane, median RNA UMI count as a function\n")
cat("of MOI (number of distinct guides detected per cell) — for all cells and for singlets only.\n")
cat("The dashed vertical line is the CURRENT config.qc.max_moi guess (", MOI_FILTER_LINE,
    "). Cells past\n")
cat("this line usually have implausibly many guides for a real single cell — likely missed\n")
cat("multiplets. Also check figs/fig1_MOI_doublet_count.png / fig2_MOI_doublet_percent.png for\n")
cat("where the doublet rate climbs sharply with MOI. Decide the final min/max MOI to keep and\n")
cat("tell Claude config.qc.min_moi / config.qc.max_moi. Once you update config.yaml, re-run this\n")
cat("plot.R to regenerate fig3 with the confirmed line — before step4 is written.\n")
