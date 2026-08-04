# scripts/scaffold.R
#
# Scaffolds a new scRNA-seq/Perturb-seq clustering pipeline run directory from
# this skill's templates/. This script is NOT executed automatically by the
# skill — it is handed to you (or to Claude acting on your behalf) to run
# interactively, one step at a time, matching the pipeline's "hard stop after
# every step" design: it never scaffolds more than one NEW step folder per
# call, because each step's figures/values must be inspected and confirmed
# before the next step's code makes sense to write.
#
# ---------------------------------------------------------------------------
# Invocation — deliberately TWO separate calls, in this order, never combined
# ---------------------------------------------------------------------------
# Source this file once per R session, then call its two functions directly:
#
#   source("scripts/scaffold.R")
#
#   # 1) One-time initial scaffold: creates target_dir, then copies .here,
#   #    CLAUDE.md, init.R, and config.yaml (from config.yaml.template, with
#   #    its `profile:` field corrected to match the argument below). Copies
#   #    NO step folder yet — see why below.
#   scaffold_init(profile = "perturbseq", target_dir = "/abs/path/to/new/run")
#
#   # <-- Now edit target_dir/config.yaml with the run's real values (lanes,
#   #     h5 paths, gRNA CSVs, marker panel, slurm.*/conda.* contact info)
#   #     BEFORE the next call. This matters: slurm.sh's __PLACEHOLDER__
#   #     tokens are substituted from config.yaml at the moment a step folder
#   #     is copied, and sbatch scripts can't read YAML themselves at submit
#   #     time — if you copy a step before config.yaml has its final values,
#   #     that step's slurm.sh permanently bakes in stale/placeholder text
#   #     and has to be manually fixed or re-copied. Do the edit first.
#
#   # 2) Every turn (including the very first step): copies exactly one more
#   #    step folder into target_dir, substituting slurm.sh from config.yaml
#   #    as it currently stands.
#   scaffold_next_step(target_dir = "/abs/path/to/new/run",
#                       step_folder_name = "step2_assign")
#   # step_folder_name may be omitted — it then auto-picks the next
#   # not-yet-scaffolded step folder, in step-number order (this is how you
#   # scaffold step1 itself: call it once with no steps yet on disk).
#
# Both functions can also be run non-interactively:
#   Rscript scripts/scaffold.R perturbseq /abs/path/to/new/run
#     -> runs scaffold_init("perturbseq", "/abs/path/to/new/run")
#   Rscript scripts/scaffold.R next /abs/path/to/new/run step2_assign
#     -> runs scaffold_next_step("/abs/path/to/new/run", "step2_assign")
#   Rscript scripts/scaffold.R next /abs/path/to/new/run
#     -> runs scaffold_next_step("/abs/path/to/new/run")  (auto-picks next)
#
# Neither function ever calls sbatch/Rscript on your behalf — after copying a
# step's files, each prints the exact command you (or Claude) should run next.
# This mirrors the skill's core rule: it only ever hands you code to run
# yourself, it never executes anything.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(glue)
})

# glue() trims exactly one leading/trailing newline from its template by
# default (`.trim = TRUE`), which silently eats the "\n" in every
# `cat(glue("...\n"))` call below. Every such call in this file relies on
# literal \n placement for output formatting, so disable that trimming here
# rather than special-casing every call site.
glue <- function(..., .envir = parent.frame()) glue::glue(..., .trim = FALSE, .envir = .envir)

# =============================================================================
# Locate this skill's templates/ directory, regardless of how/where this
# script is sourced or Rscript'd from.
# =============================================================================
.find_skill_root <- function() {
  # Manual escape hatch, in case auto-detection below guesses wrong (e.g.
  # this file was copied out of the skill folder before being run).
  env_override <- Sys.getenv("SCRNA_SKILL_ROOT", unset = NA)
  if (!is.na(env_override) && nzchar(env_override)) {
    return(normalizePath(env_override, mustWork = TRUE))
  }
  # Case: `source("scripts/scaffold.R")` — inspect the call stack for a
  # literal `source(...)` call and evaluate its `file` argument. This works
  # regardless of the `keep.source` option (srcref-based detection does not).
  # Checked BEFORE the --file= case below on purpose: if this file was
  # source()'d from inside another Rscript-invoked file (e.g. a wrapper or
  # test harness run via `Rscript wrapper.R`), commandArgs()'s --file= would
  # point at that OUTER wrapper, not at scaffold.R itself — the source() call
  # stack is the more specific signal and must win whenever both are present.
  calls <- sys.calls()
  for (i in seq_along(calls)) {
    call_i <- calls[[i]]
    if (is.call(call_i) && identical(call_i[[1]], quote(source))) {
      matched <- tryCatch(match.call(base::source, call_i), error = function(e) NULL)
      if (!is.null(matched) && !is.null(matched[["file"]])) {
        parent_env <- if (i > 1) sys.frame(i - 1) else globalenv()
        file_val <- tryCatch(eval(matched[["file"]], envir = parent_env), error = function(e) NULL)
        if (is.character(file_val) && length(file_val) == 1 && nzchar(file_val)) {
          return(normalizePath(file.path(dirname(file_val), ".."), mustWork = TRUE))
        }
      }
    }
  }
  # Case: `Rscript scripts/scaffold.R ...` (this file itself is the top-level
  # script — no enclosing source() call was found above).
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    this_file <- sub("^--file=", "", file_arg)
    return(normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE))
  }
  # Last-resort fallback: assume the current working directory IS the skill
  # root (e.g. you `setwd()`'d into the skill folder before sourcing this).
  # If this guess is wrong, set the SCRNA_SKILL_ROOT environment variable
  # instead, or edit SKILL_ROOT below by hand.
  warning(
    "Could not auto-detect the skill root directory (templates/ location). ",
    "Falling back to getwd(). If templates/ is not found from there, set ",
    "Sys.setenv(SCRNA_SKILL_ROOT = \"/path/to/scrna-seq-clustering-pipeline\") ",
    "before sourcing this file."
  )
  normalizePath(getwd(), mustWork = TRUE)
}

SKILL_ROOT    <- .find_skill_root()
TEMPLATES_DIR <- file.path(SKILL_ROOT, "templates")

# =============================================================================
# Internal helpers
# =============================================================================

# List step folders for a profile, in step-number order, by reading
# templates/<profile>/ directly — never hardcode step names, so a new step
# folder added later to the templates just works.
.list_step_folders <- function(profile) {
  d <- file.path(TEMPLATES_DIR, profile)
  if (!dir.exists(d)) stop(glue("No such profile template folder: {d}"))
  entries <- list.dirs(d, full.names = FALSE, recursive = FALSE)
  entries <- entries[grepl("^step[0-9]+", entries)]
  if (length(entries) == 0) stop(glue("No stepN_* folders found under {d}"))
  # Sort by the leading step number, not alphabetically (avoids "step10" <
  # "step2" string-sort bugs).
  step_num <- as.integer(sub("^step([0-9]+).*$", "\\1", entries))
  entries[order(step_num)]
}

.detect_profile <- function(target_dir) {
  config_path <- file.path(target_dir, "config.yaml")
  if (!file.exists(config_path)) {
    stop(glue("{config_path} not found — run scaffold_init() first."))
  }
  config <- yaml::read_yaml(config_path)
  profile <- config$profile
  if (is.null(profile) || !profile %in% c("perturbseq", "plain")) {
    stop("config.yaml's top-level `profile:` key must be 'perturbseq' or 'plain'.")
  }
  profile
}

# Which step folders (from the profile's full ordered list) already exist
# under target_dir/?
.already_scaffolded_steps <- function(target_dir, profile) {
  all_steps <- .list_step_folders(profile)
  all_steps[dir.exists(file.path(target_dir, all_steps))]
}

# Substitute every slurm.* / conda.* placeholder token in a slurm.sh's text,
# using values read from target_dir/config.yaml. sbatch scripts cannot read
# YAML themselves at submit time, so this substitution has to happen now, at
# scaffold time, not later.
.substitute_slurm_placeholders <- function(text, config, step_folder_name) {
  res <- config$slurm$resources[[step_folder_name]]
  if (is.null(res)) {
    warning(glue(
      "No config$slurm$resources entry for '{step_folder_name}' in config.yaml ",
      "— __CPUS__/__MEM__/__TIME__ will be left as literal placeholder text. ",
      "Add an entry under slurm.resources.{step_folder_name} in config.yaml ",
      "and re-scaffold this step (delete the folder first) if that happens."
    ))
    res <- list(cpus = "__CPUS__", mem = "__MEM__", time = "__TIME__")
  }
  repl <- c(
    "__MAIL_USER__"    = as.character(config$slurm$mail_user    %||% "__MAIL_USER__"),
    "__LOG_DIR__"      = as.character(config$slurm$log_dir      %||% "__LOG_DIR__"),
    "__PARTITION__"    = as.character(config$slurm$partition    %||% "__PARTITION__"),
    "__PROJECT_ROOT__" = as.character(config$slurm$project_root %||% "__PROJECT_ROOT__"),
    "__CONDA_ENV__"    = as.character(config$conda$env_name     %||% "__CONDA_ENV__"),
    "__CPUS__"         = as.character(res$cpus %||% "__CPUS__"),
    "__MEM__"          = as.character(res$mem  %||% "__MEM__"),
    "__TIME__"         = as.character(res$time %||% "__TIME__")
  )
  for (token in names(repl)) {
    text <- gsub(token, repl[[token]], text, fixed = TRUE)
  }
  text
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Copy one step folder's flat files (process.R / plot.R / slurm.sh — whichever
# are present; step10/step11-equivalents don't have all three, by design) from
# templates/<profile>/<step_folder_name>/ into
# target_dir/<step_folder_name>/scripts/, substituting slurm.sh placeholders
# on the way. Skips (with a message) any destination file that already
# exists, so re-running scaffold_next_step() on an already-scaffolded step is
# a safe no-op rather than clobbering edits you've made.
.copy_step_folder <- function(profile, target_dir, step_folder_name) {
  src_dir <- file.path(TEMPLATES_DIR, profile, step_folder_name)
  if (!dir.exists(src_dir)) {
    available <- paste(.list_step_folders(profile), collapse = ", ")
    stop(glue(
      "No such template step folder: {src_dir}\n",
      "Available step folders for profile '{profile}': {available}"
    ))
  }
  dest_scripts_dir <- file.path(target_dir, step_folder_name, "scripts")
  dir.create(dest_scripts_dir, recursive = TRUE, showWarnings = FALSE)

  config_path <- file.path(target_dir, "config.yaml")
  config <- if (file.exists(config_path)) yaml::read_yaml(config_path) else NULL

  src_files <- list.files(src_dir, full.names = FALSE)
  for (f in src_files) {
    src_path  <- file.path(src_dir, f)
    dest_path <- file.path(dest_scripts_dir, f)
    if (file.exists(dest_path)) {
      cat(glue("  [skip] {step_folder_name}/scripts/{f} already exists — not overwriting\n"))
      next
    }
    if (identical(tolower(f), "slurm.sh")) {
      if (is.null(config)) {
        stop(glue(
          "{target_dir}/config.yaml not found — cannot substitute slurm.sh ",
          "placeholders. Run scaffold_init() first (it creates config.yaml from ",
          "config.yaml.template)."
        ))
      }
      raw_text <- paste(readLines(src_path, warn = FALSE), collapse = "\n")
      sub_text <- .substitute_slurm_placeholders(raw_text, config, step_folder_name)
      writeLines(sub_text, dest_path)
      Sys.chmod(dest_path, "0755")
    } else {
      file.copy(src_path, dest_path, overwrite = FALSE)
    }
    cat(glue("  [ok]   {step_folder_name}/scripts/{f}\n"))
  }
  invisible(dest_scripts_dir)
}

# Print the exact command(s) the user should run next for a freshly-copied
# step, and remind them the skill will not run it for them.
.print_next_command <- function(target_dir, step_folder_name) {
  scripts_dir <- file.path(target_dir, step_folder_name, "scripts")
  has_slurm   <- file.exists(file.path(scripts_dir, "slurm.sh"))
  has_process <- file.exists(file.path(scripts_dir, "process.R"))
  has_plot    <- file.exists(file.path(scripts_dir, "plot.R"))

  cat("\n---\n")
  cat(glue("Scaffolded {step_folder_name}/scripts/. This skill will NOT run anything for you.\n\n"))
  if (has_slurm) {
    cat("Submit it yourself (from inside the run directory, i.e. slurm.sh's own `cd` target):\n\n")
    cat(glue("  sbatch {step_folder_name}/scripts/slurm.sh\n\n"))
    cat("That runs process.R then plot.R in sequence and emails you on completion (per\n")
    cat("--mail-type=BEGIN,END,FAIL). Check the .log/.err files in slurm.log_dir when it's done.\n")
  } else {
    cat("This step has no slurm.sh — it's meant to be run interactively/manually. Run directly:\n\n")
    if (has_process) cat(glue("  Rscript {step_folder_name}/scripts/process.R\n"))
    if (has_plot)    cat(glue("  Rscript {step_folder_name}/scripts/plot.R\n"))
    cat("\n")
  }
  cat("Once it finishes, inspect its figs/tables per the plot.R script's STOP message, report\n")
  cat("the requested value (if any) back to Claude, and Claude will update config.yaml and call\n")
  cat("scaffold_next_step() for the following step.\n")
}

.claude_md_template <- function(profile) {
  bullets <- c(
    "Every tunable number/path for this pipeline lives in `config.yaml`, not hardcoded in any",
    "step script. Do not hand-edit constants inside `stepN_xxx/scripts/*.R` files — edit",
    "`config.yaml` instead and re-run. See `config.yaml`'s own comments for what each key",
    "controls and which step reads it."
  )
  init_bullets <- c(
    "- Packages, plotting theme, palette (`pal <- brewer.pal(9, \"Set1\")`)",
    "- `config` — the parsed config.yaml (`config <- yaml::read_yaml(here(\"config.yaml\"))`)",
    "- `n_cores`, `lane_names` (from config.yaml)",
    "- `markers`, `cell_type_map`, `ct_levels` (from config.yaml's `markers.*` block)"
  )
  if (identical(profile, "perturbseq")) {
    init_bullets <- c(init_bullets,
      "- `nc_set`, `oeg_set`, `oeg_gene_map` (perturbseq only — guide-capture CSVs)")
  }
  seurat_bullets <- c(
    "- **Never call `JoinLayers()`** unless you have a specific reason to and understand why —",
    "  the original validated pipeline never needed it.",
    "- **Never copy clusters via `obj$seurat_clusters <- x` expecting `Idents()` to update** —",
    "  `$<-` does NOT update `Idents()` in Seurat v5. Always call `Idents(obj) <- x` explicitly",
    "  if downstream code (e.g. `DimPlot`) relies on `Idents()`."
  )
  if (identical(profile, "perturbseq")) {
    seurat_bullets <- c(seurat_bullets,
      "- **RNA assay rownames are gene symbols**, not Ensembl IDs — match genes on symbol,",
      "  not on any `*_id` column in your guide CSVs.")
  }
  c(
    "# CLAUDE.md",
    "",
    "This file provides guidance to Claude Code (claude.ai/code) when working with code in this",
    "repository.",
    "",
    glue("This run was scaffolded by the scrna-seq-clustering-pipeline skill, profile: `{profile}`."),
    "",
    "## Running scripts",
    "",
    "```bash",
    "conda activate <config.yaml conda.env_name>",
    "cd <this directory>",
    "Rscript stepN_xxx/scripts/process.R   # data processing (a few steps have no process.R)",
    "Rscript stepN_xxx/scripts/plot.R      # figures & tables",
    "```",
    "",
    "Or via SLURM, where a step has a `slurm.sh`:",
    "```bash",
    "sbatch stepN_xxx/scripts/slurm.sh",
    "```",
    "",
    "## config.yaml — single source of truth",
    "",
    bullets,
    "",
    "## init.R — what it provides",
    "",
    "`source(here::here(\"init.R\"))` at the top of every script. Do not repeat any of the",
    "following in a step script:",
    "",
    init_bullets,
    "",
    "## Pipeline design — one step at a time",
    "",
    "This run is scaffolded ONE step folder at a time, not all at once (via this skill's",
    "`scripts/scaffold.R`). Most steps end their `plot.R` with a `STOP:` message asking you to",
    "inspect a figure and report a value back before the next step's `config.yaml` fields are",
    "filled in and its folder is scaffolded. Do not skip ahead or hand-write a later step's",
    "scripts yourself — regenerate them from the templates once the checkpoint is confirmed.",
    "",
    "## Seurat pipeline constraints",
    "",
    seurat_bullets,
    "",
    "## Plotting",
    "",
    "- Use Seurat plotting functions first (`VlnPlot`, `FeaturePlot`, `DimPlot`, etc.) with the",
    "  Seurat default palette.",
    "- Fall back to ggplot2 + `pal` only when Seurat plotting is not feasible.",
    "- Always pass `bg = \"white\"` to every `ggsave()` call.",
    "",
    "## plot.R conventions",
    "",
    "Each step's `stepN_xxx/scripts/plot.R` is the sole entry point for figures and tables.",
    "Do not duplicate any figure or table output in `process.R` if `plot.R` already produces it.",
    "",
    "**Output structure**:",
    "- `stepN_xxx/figs/figN_<name>.png` — PNG figures",
    "- `stepN_xxx/tables/tblN_<name>.csv` — CSV tables",
    ""
  )
}

# =============================================================================
# scaffold_init(): one-time initial scaffold of a brand-new run directory.
# Deliberately does NOT copy any step folder — see the file header comment
# for why (slurm.sh placeholder substitution needs config.yaml's FINAL
# values, and this function only creates config.yaml, it doesn't fill it in).
# Call scaffold_next_step() afterwards, once config.yaml has been edited, to
# copy step1 (and then every subsequent step).
# =============================================================================
scaffold_init <- function(profile = c("perturbseq", "plain"), target_dir) {
  profile <- match.arg(profile)
  stopifnot(is.character(target_dir), nzchar(target_dir))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  target_dir <- normalizePath(target_dir, mustWork = TRUE)

  cat(glue("Initializing new '{profile}' run at {target_dir} ...\n\n"))

  # .here (empty file — marks this directory as the here::here() root) -------
  here_path <- file.path(target_dir, ".here")
  if (!file.exists(here_path)) {
    file.create(here_path)
    cat("  [ok]   .here\n")
  } else {
    cat("  [skip] .here already exists\n")
  }

  # CLAUDE.md ----------
  claude_md_path <- file.path(target_dir, "CLAUDE.md")
  if (!file.exists(claude_md_path)) {
    writeLines(.claude_md_template(profile), claude_md_path)
    cat("  [ok]   CLAUDE.md\n")
  } else {
    cat("  [skip] CLAUDE.md already exists — not overwriting\n")
  }

  # init.R ----------
  init_dest <- file.path(target_dir, "init.R")
  if (!file.exists(init_dest)) {
    file.copy(file.path(TEMPLATES_DIR, profile, "init.R"), init_dest)
    cat("  [ok]   init.R\n")
  } else {
    cat("  [skip] init.R already exists — not overwriting\n")
  }

  # config.yaml (from config.yaml.template) ----------
  # There is only one config.yaml.template, shared by both profiles (its own
  # comments explain why), and it ships with `profile: perturbseq` literally
  # in the text. Fix that line at copy time so a "plain" scaffold doesn't
  # silently end up with the wrong profile recorded — scaffold_next_step()
  # reads this exact field to decide which templates/<profile>/ folder to
  # pull the next step from, so getting it wrong would copy the wrong
  # pipeline's steps. Done as a targeted line substitution, not a full
  # yaml::read_yaml()/write_yaml() round-trip, so config.yaml.template's
  # extensive `#` comments (the main documentation a user reads) survive
  # intact — write_yaml() would silently drop every comment.
  config_dest <- file.path(target_dir, "config.yaml")
  if (!file.exists(config_dest)) {
    template_lines <- readLines(file.path(TEMPLATES_DIR, "config.yaml.template"), warn = FALSE)
    template_lines <- sub("^profile:\\s*perturbseq\\s*$", glue("profile: {profile}"), template_lines)
    writeLines(template_lines, config_dest)
    cat(glue("  [ok]   config.yaml (profile: {profile}) — EDIT THIS before scaffolding step1!\n"))
  } else {
    cat("  [skip] config.yaml already exists — not overwriting\n")
  }

  cat("\n=== INIT COMPLETE — no step folder copied yet ===\n")
  cat(glue("Next: open {config_dest} and fill in at minimum\n"))
  cat("paths.h5_dir / paths.h5_pattern / lanes, slurm.project_root/mail_user/log_dir/partition\n")
  cat("(needed to fill in every slurm.sh's placeholders), and the markers.* panel for your tissue\n")
  cat("(the shipped example panel is cardiac-progenitor-specific and will not apply to your data).\n")
  cat("THEN call scaffold_next_step(target_dir) to copy step1 — do this only after editing\n")
  cat("config.yaml, since step1's slurm.sh gets its placeholders substituted from whatever is in\n")
  cat("config.yaml at the moment it's copied.\n")
  invisible(target_dir)
}

# =============================================================================
# scaffold_next_step(): copy exactly one more step folder into an existing run
# =============================================================================
scaffold_next_step <- function(target_dir, step_folder_name = NULL) {
  stopifnot(dir.exists(target_dir))
  target_dir <- normalizePath(target_dir, mustWork = TRUE)
  profile <- .detect_profile(target_dir)

  if (is.null(step_folder_name)) {
    all_steps <- .list_step_folders(profile)
    done      <- .already_scaffolded_steps(target_dir, profile)
    remaining <- setdiff(all_steps, done)
    if (length(remaining) == 0) {
      cat(glue("All {length(all_steps)} step folders for profile '{profile}' are already scaffolded.\n"))
      return(invisible(target_dir))
    }
    step_folder_name <- remaining[1]
    cat(glue("No step_folder_name given — auto-picked next unscaffolded step: {step_folder_name}\n\n"))
  }

  cat(glue("Scaffolding next step ({step_folder_name}) for '{profile}' run at {target_dir} ...\n\n"))
  .copy_step_folder(profile, target_dir, step_folder_name)
  cat("\n=== NEXT STEP SCAFFOLDED ===\n")
  .print_next_command(target_dir, step_folder_name)
  invisible(target_dir)
}

# =============================================================================
# Rscript entry point (only fires when THIS file is the top-level script
# being run via `Rscript scripts/scaffold.R ...`; sourcing it — interactively,
# or from inside another Rscript-invoked file, e.g. a wrapper/test harness —
# just defines the two functions above with no side effects).
#
# Checking for "--file=" alone is not enough: commandArgs(trailingOnly=FALSE)
# reflects how the CURRENT R process was launched, so if some other script is
# run via `Rscript wrapper.R args...` and that wrapper does
# `source("scaffold.R")`, the check below would otherwise be TRUE and this
# block would misinterpret wrapper.R's own arguments as scaffold.R's CLI
# arguments. Guard against that by requiring the "--file=" target to actually
# resolve to THIS file (scaffold.R), using the same SKILL_ROOT already
# resolved above (which itself correctly tells the two invocation styles
# apart) rather than trusting commandArgs() alone.
# =============================================================================
.this_file       <- file.path(SKILL_ROOT, "scripts", "scaffold.R")
.rscript_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.is_rscript_invocation <- length(.rscript_file_arg) == 1 &&
  file.exists(.this_file) &&
  normalizePath(sub("^--file=", "", .rscript_file_arg), mustWork = FALSE) ==
    normalizePath(.this_file, mustWork = FALSE)
if (.is_rscript_invocation) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (length(.args) >= 2 && .args[1] %in% c("perturbseq", "plain")) {
    scaffold_init(profile = .args[1], target_dir = .args[2])
  } else if (length(.args) >= 2 && .args[1] == "next") {
    scaffold_next_step(target_dir = .args[2],
                        step_folder_name = if (length(.args) >= 3) .args[3] else NULL)
  } else if (length(.args) > 0) {
    stop(paste(
      "Usage:\n",
      "  Rscript scaffold.R <perturbseq|plain> <target_dir>\n",
      "  Rscript scaffold.R next <target_dir> [step_folder_name]\n"
    ))
  }
  # 0 args: nothing to do when Rscript'd directly with no args — functions
  # are defined but this is an unusual way to invoke it (normally you'd
  # `source()` this file and call scaffold_init()/scaffold_next_step() yourself).
}
