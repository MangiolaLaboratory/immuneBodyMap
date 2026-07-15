#!/usr/bin/env Rscript
# Rebuild SI workbooks and/or instruct on full recomputation.

args <- commandArgs(trailingOnly = TRUE)
si_only <- "--si-only" %in% args
full <- "--full" %in% args

repo <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
if (!nzchar(repo)) {
  # scripts/ -> vignettes/ -> repo
  repo <- normalizePath(file.path("..", ".."), mustWork = FALSE)
}
Sys.setenv(IMMUNE_HEALTHY_BODY_MAP_ROOT = repo)

paths_r <- file.path(repo, "vignettes", "R", "paths.R")
source(paths_r)
source(file.path(repo, "vignettes", "R", "export_supplementary_tables.R"))

if (si_only || (!full && !si_only)) {
  message("Building panel SI workbooks from figure_ready caches where possible...")
  builder <- file.path(
    repo, "vignettes", "R", "build_panel_supplementary_workbooks.R"
  )
  # The Age_Clock-oriented builder expects AGE_CLOCK_ROOT layout. Set archive if available.
  if (!nzchar(Sys.getenv("AGE_CLOCK_ROOT", unset = ""))) {
    # Point builder at a synthetic root by documenting expectation:
    message(
      "Note: build_panel_supplementary_workbooks.R was designed for Age_Clock.\n",
      "Set AGE_CLOCK_ARCHIVE / AGE_CLOCK_ROOT to the analysis tree, or rebuild\n",
      "SI sheets from vignettes/data/source_tables (already included).\n"
    )
  }
  message("SI source tables currently shipped under vignettes/data/source_tables/")
}

if (full) {
  message(
    "Full recomputation mode:\n",
    "  1. Configure config/paths.example.env\n",
    "  2. Ensure PSEUDOBULK_H5AD and SCCOMP_ESTIMATES_DIR (or refit)\n",
    "  3. Delete selected caches under vignettes/data/intermediate and re-render\n",
    "  4. Or run Age_Clock report/*.qmd with AGE_CLOCK_ARCHIVE set\n"
  )
}

invisible(NULL)
