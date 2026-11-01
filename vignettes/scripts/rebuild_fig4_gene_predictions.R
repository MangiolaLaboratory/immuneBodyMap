#!/usr/bin/env Rscript
# Rebuild Fig4 gene-expression ageing predictions (pred_fixed.rds / pred_total.rds)
# from NEW_CELL_TYPE brms targets stores.
#
# Primary inputs:
#   - NEW_CELL_TYPE_ROOT  : directory of V*_celltype/_targets stores (exemption)
#   - gene_hyp_df         : small gene-hypothesis table (shipped under
#                           vignettes/data/intermediate/fig4/new_gene_hyp_df.rds)
#
# Outputs (under NEW_PRED_DIR, default vignettes/data/intermediate/fig4):
#   - pred_fixed.rds
#   - pred_total.rds
#
# Stomach / missing-tissue extensions remain under AGE_CLOCK_ARCHIVE:
#   TAR_stomach/stomach_pseudobulk_data.R
#   TAR_missing_tissue/missing_tissue_data.R
# See rebuild_fig4_gene_predictions.md in this directory.
#
# Portable adaptation of Age_Clock/new_pred/pred_fixed_total.R

suppressPackageStartupMessages({
  library(tidySummarizedExperiment)
  library(HDF5Array)
  library(magrittr)
  library(targets)
  library(tidyverse)
  library(tidybulk)
  library(brms)
})

repo <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
if (!nzchar(repo)) {
  repo <- normalizePath(file.path("..", ".."), mustWork = FALSE)
}
source(file.path(repo, "vignettes", "R", "paths.R"))

new_cell_type_root <- Sys.getenv("NEW_CELL_TYPE_ROOT", unset = "")
if (!nzchar(new_cell_type_root) || !dir.exists(new_cell_type_root)) {
  stop(
    "Set NEW_CELL_TYPE_ROOT to the NEW_CELL_TYPE directory containing ",
    "V*_celltype/_targets stores (targets-cache exemption).",
    call. = FALSE
  )
}

out_dir <- Sys.getenv(
  "NEW_PRED_DIR",
  unset = path_intermediate("fig4")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_hyp_path <- file.path(path_intermediate("fig4"), "new_gene_hyp_df.rds")
if (!file.exists(gene_hyp_path)) {
  stop("Missing gene hypothesis table: ", gene_hyp_path, call. = FALSE)
}
gene_hyp_df <- readRDS(gene_hyp_path)

get_brms_tar <- function(cell_type, gene_id, base_path = new_cell_type_root) {
  cli::cli_alert_info("Locating cell type `{cell_type}` under {base_path} ...")

  ct_list <- list.dirs(base_path, recursive = FALSE) %>% basename()
  ct_list <- ct_list[ct_list %>% stringr::str_starts("V[0-9]_")]
  ct_name <- stringr::str_remove(ct_list, "V[0-9]_")

  if (!make.names(cell_type) %in% ct_name) {
    stop("Cell type not found: ", cell_type, call. = FALSE)
  }

  ct_list <- ct_list[ct_list %>% stringr::str_detect(make.names(cell_type))]
  version_no <- ct_list %>%
    stringr::str_extract("^V\\d+") %>%
    stringr::str_remove("^V") %>%
    as.integer() %>%
    max(na.rm = TRUE)

  targets_path <- glue::glue(
    "{base_path}/V{version_no}_{make.names(cell_type)}/_targets"
  )
  cli::cli_alert_info("Using version V{version_no} at {targets_path}")

  mapping_table <- readr::read_delim(
    glue::glue("{targets_path}/meta/meta"),
    delim = "|",
    show_col_types = FALSE
  ) %>%
    dplyr::filter(name %>% stringr::str_starts("estimates_chunk")) %>%
    dplyr::filter(type == "branch" & is.na(error)) %>%
    dplyr::mutate(
      gene = warnings %>%
        stringr::str_extract("(?<=Gene:___)ENSG\\d+(?=___)")
    )

  if (!gene_id %in% mapping_table$gene) {
    stop("Gene ID not found in mapping: ", gene_id, call. = FALSE)
  }

  target_name <- mapping_table %>%
    dplyr::filter(gene == gene_id) %>%
    dplyr::pull(name)

  if (length(target_name) != 1) {
    stop("Expected 1 target for gene ", gene_id, ", found ", length(target_name),
         call. = FALSE)
  }

  brms_path <- glue::glue("{targets_path}/objects/{target_name}")
  if (!file.exists(brms_path)) {
    stop("Target file not found: ", brms_path, call. = FALSE)
  }

  qs2::qs_read(brms_path)
}

n_workers <- as.integer(Sys.getenv("FIG4_PRED_WORKERS", unset = "4"))

# --- fixed effects -----------------------------------------------------------
short_list_fixed <- gene_hyp_df %>%
  dplyr::filter(component == "fixed") %>%
  dplyr::group_by(cell_type, .feature) %>%
  dplyr::summarise(
    n = dplyr::n(),
    max_abs_est = max(abs_est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(max_abs_est), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

message("Building pred_fixed.rds (", nrow(short_list_fixed), " gene/cell-type rows)...")
future::plan(future::multisession, workers = n_workers)
pred_fixed <- short_list_fixed %>%
  furrr::future_pmap_dfr(~ {
    get_brms_tar(cell_type = ..1, gene_id = ..2) %>%
      dplyr::mutate(pred_df = purrr::map(brms_fit, ~ {
        newdf <- .x$data %>%
          dplyr::distinct(age_decade, sex, ethnicity_groups) |>
          dplyr::arrange(age_decade) %>%
          dplyr::mutate(
            disease_groups_altered = "Normal",
            tissue_groups = NA,
            assay_groups_altered = NA,
            dataset_id_altered = NA,
            offset = 0
          )
        pred <- .x %>% stats::predict(
          newdata = newdf,
          summary = FALSE,
          re_formula = ~ age_decade * sex + ethnicity_groups,
          allow_new_levels = TRUE
        )
        SummarizedExperiment::SummarizedExperiment(
          assays = list(counts = pred),
          colData = newdf
        ) %>%
          tibble::as_tibble() %>%
          dplyr::rename(draws = .feature) %>%
          dplyr::group_by(age_decade, draws) %>%
          dplyr::reframe(count_avg = mean(counts, na.rm = TRUE))
      })) %>%
      dplyr::select(-brms_fit)
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

saveRDS(pred_fixed, file.path(out_dir, "pred_fixed.rds"))
message("Wrote ", file.path(out_dir, "pred_fixed.rds"))

# --- total (tissue) effects --------------------------------------------------
short_list_total <- gene_hyp_df %>%
  dplyr::filter(component == "total") %>%
  dplyr::group_by(tissue, cell_type, .feature) %>%
  dplyr::summarise(
    n = dplyr::n(),
    max_abs_est = max(abs_est, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(tissue, cell_type) %>%
  dplyr::arrange(dplyr::desc(n), dplyr::desc(max_abs_est), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

message("Building pred_total.rds (", nrow(short_list_total), " rows)...")
pred_total <- short_list_total %>%
  furrr::future_pmap_dfr(~ {
    get_brms_tar(cell_type = ..2, gene_id = ..3) %>%
      dplyr::mutate(pred_df = purrr::map(brms_fit, ~ {
        newdf <- .x$data %>%
          dplyr::distinct(age_decade, tissue_groups, sex, ethnicity_groups) |>
          dplyr::arrange(age_decade) %>%
          dplyr::mutate(
            disease_groups_altered = "Normal",
            assay_groups_altered = NA,
            dataset_id_altered = NA,
            offset = 0
          )
        pred <- .x %>% stats::predict(
          newdata = newdf,
          summary = FALSE,
          re_formula = ~ age_decade * sex + ethnicity_groups +
            (age_decade * sex + ethnicity_groups | tissue_groups),
          allow_new_levels = FALSE
        )
        SummarizedExperiment::SummarizedExperiment(
          assays = list(counts = pred),
          colData = newdf
        ) %>%
          tibble::as_tibble() %>%
          dplyr::rename(draws = .feature) %>%
          dplyr::group_by(age_decade, tissue_groups, draws) %>%
          dplyr::reframe(count_avg = mean(counts, na.rm = TRUE))
      })) %>%
      dplyr::select(-brms_fit)
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

saveRDS(pred_total, file.path(out_dir, "pred_total.rds"))
message("Wrote ", file.path(out_dir, "pred_total.rds"))
message(
  "Done. For stomach / missing-tissue extensions see ",
  "vignettes/scripts/rebuild_fig4_gene_predictions.md"
)
