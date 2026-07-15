#!/usr/bin/env Rscript
# Driver: run the per–cell-type gene-ageing targets pipeline (NEW_CELL_TYPE layout).
#
# Configure (see also vignettes/config/paths.example.env):
#   IMMUNE_HEALTHY_BODY_MAP_ROOT  - repo root
#   NEW_CELL_TYPE_ROOT            - output root for V1_<celltype>/_targets stores
#   PSEUDOBULK_H5AD               - primary input (Zenodo)
#   TAR_PIPELINE_ROOT             - optional override of this directory
#   TARGETS_USE_CREW              - true (HPC Slurm) / false (local)
#   TARGETS_MAKE_NAMES            - optional; default "adjust_age"
#                                  also useful: hypothesis_age_monotonic_and_adjust_tissue
#
# Usage:
#   export NEW_CELL_TYPE_ROOT=/path/to/NEW_CELL_TYPE
#   export PSEUDOBULK_H5AD=/path/to/pseudobulk_se.h5ad
#   Rscript vignettes/targets_pipeline/batch_run.R

library(targets)
library(tidyverse)
library(cli)
library(glue)

cell_type_list <- c(
  "plasma",
  "nk",
  "cd4 naive",
  "cd8 naive",
  "cd4 th1/th17 em",
  "cd4 th17 em",
  "mait",
  "tgd",
  "cd4 fh em",
  "cd4 th2 em",
  "cd8 tem",
  "cd4 tcm",
  "b memory",
  "b naive",
  "macrophage",
  "cd8 tcm",
  "treg",
  "monocytic"
)

resolve_new_cell_type_root <- function() {
  env <- Sys.getenv("NEW_CELL_TYPE_ROOT", unset = "")
  if (!nzchar(env)) {
    stop(
      "Set NEW_CELL_TYPE_ROOT to the directory that will hold ",
      "V1_<celltype>/_targets stores.",
      call. = FALSE
    )
  }
  dir.create(env, recursive = TRUE, showWarnings = FALSE)
  normalizePath(env, winslash = "/", mustWork = TRUE)
}

resolve_script_path <- function() {
  env <- Sys.getenv("TAR_PIPELINE_ROOT", unset = "")
  if (nzchar(env)) {
    p <- file.path(env, "dynamic_tar_script.R")
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  root <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
  if (nzchar(root)) {
    p <- file.path(root, "vignettes", "targets_pipeline", "dynamic_tar_script.R")
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  # Same directory as this script when sourced / Rscript'd via relative path
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    here <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE))
    p <- file.path(here, "dynamic_tar_script.R")
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Could not locate dynamic_tar_script.R. Set TAR_PIPELINE_ROOT or ",
    "IMMUNE_HEALTHY_BODY_MAP_ROOT.",
    call. = FALSE
  )
}

root_path <- resolve_new_cell_type_root()
script_path <- resolve_script_path()
Sys.setenv(TAR_PIPELINE_ROOT = dirname(script_path))

make_names <- Sys.getenv("TARGETS_MAKE_NAMES", unset = "adjust_age")

cli_alert_info("NEW_CELL_TYPE_ROOT: {root_path}")
cli_alert_info("dynamic_tar_script: {script_path}")
cli_alert_info("tar_make names: {make_names}")

done_list <- list()
error_list <- list()

cli_progress_bar("Processing cell types", total = length(cell_type_list))

for (cur_ct in cell_type_list) {
  # expose cell type for dynamic_tar_script.R
  assign("cur_ct", cur_ct, envir = .GlobalEnv)

  setwd(root_path)

  tryCatch(
    {
      cli_alert_info(glue("Working on {cur_ct}..."))

      cur_path <- glue("{root_path}/V1_{make.names(cur_ct)}")
      if (!dir.exists(cur_path)) {
        dir.create(cur_path, recursive = TRUE)
      }
      cli_alert_success(glue("{cur_ct} directory ready: {cur_path}"))

      setwd(cur_path)
      source(script_path)
      tar_option_set(debug = character(0))
      tar_make(
        names = make_names,
        shortcut = TRUE,
        callr_function = NULL
      )

      cli_alert_success(glue("{cur_ct} ALL done!"))
      done_list <- append(done_list, cur_ct)
    },
    error = function(e) {
      cli_alert_danger(glue("Error in {cur_ct}: {e$message}"))
      error_list[[cur_ct]] <<- e$message
    }
  )

  cli_progress_update()
}

cli_progress_done()

cli_alert_info("Done: {length(done_list)} ; errors: {length(error_list)}")
if (length(error_list)) {
  print(error_list)
}
