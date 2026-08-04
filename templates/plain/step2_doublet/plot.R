# conda run -n <conda.env_name> Rscript step2_doublet/scripts/plot.R
source(here::here("init.R"))

doublet_colors <- c(singlet = "grey70", doublet = pal[2])

figs_dir   <- here("step2_doublet", "figs")
tables_dir <- here("step2_doublet", "tables")
dir.create(figs_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Derive plot data from per-lane objects ----------------------------------
cat("Reading per-lane objects...\n")
all_meta <- map(lane_names, function(l) {
  obj <- qs_read(here("step2_doublet", "obj1_doublet", glue("{l}.qs2")), nthreads = n_cores)
  cat(" ", l, "\n")
  obj@meta.data %>%
    select(scDblFinder.class) %>%
    mutate(lane = l)
}) %>%
  list_rbind() %>%
  mutate(lane = factor(lane, levels = lane_names))

# ---- tbl1: doublet summary --------------------------------------------------
# NOTE: identical logic to the original pipeline's tbl1_doublet_summary.csv —
# nothing here references MOI, so nothing to drop for the plain profile.
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

# fig1: doublet count/rate per lane ----------
# NOTE: the original pipeline's fig1 (MOI vs doublet count), fig2 (MOI vs
# doublet fraction), and fig3 (MOI vs nCount) all plot doublet behavior AS A
# FUNCTION OF MOI (number of distinct guides per cell) — a concept that does
# not exist in the plain profile (no guide capture). Replaced with a single
# simple per-lane summary figure so this step still gives the user something
# to sanity-check: a stacked singlet/doublet count bar (left) and a doublet
# rate % bar (right), one bar per lane.
p_count <- ggplot(all_meta %>% filter(!is.na(scDblFinder.class)),
                   aes(x = lane, fill = factor(scDblFinder.class, levels = c("singlet", "doublet")))) +
  geom_bar() +
  scale_fill_manual(values = doublet_colors, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Singlet / doublet count per lane", x = NULL, y = "N cells")

p_rate <- ggplot(doublet_summary, aes(x = lane, y = doublet_rate)) +
  geom_col(fill = doublet_colors[["doublet"]]) +
  scale_y_continuous(labels = function(x) glue("{x}%"),
                      expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Doublet rate per lane", x = NULL, y = "Doublet rate (%)")

p <- wrap_plots(p_count, p_rate, nrow = 1) +
  plot_annotation(title = "Doublet detection summary (scDblFinder)")

ggsave(file.path(figs_dir, "fig1_doublet_summary.png"), p, width = 12, height = 6, dpi = 300, bg = "white")
cat("Saved: figs/fig1_doublet_summary.png\n")

cat("\n=== STEP 2 PLOT.R COMPLETE ===\n")

cat("\n---\n")
cat("STOP: confirm this step finished and its outputs look reasonable before telling Claude to\n")
cat("write step3.\n")
cat("Open figs/fig1_doublet_summary.png and tables/tbl1_doublet_summary.csv. The left panel shows\n")
cat("how many cells per lane were called singlet (real single cell) vs doublet (likely two cells\n")
cat("captured together) by scDblFinder; the right panel shows the doublet rate (%) per lane —\n")
cat("typically a few percent, roughly proportional to how many cells were loaded per lane. No new\n")
cat("config value is needed here — the doublet calls themselves are used directly in step3_filter.\n")
cat("Just confirm the rates look plausible (not 0% and not wildly high, e.g. >30%), then tell\n")
cat("Claude to write step3_filter.\n")
