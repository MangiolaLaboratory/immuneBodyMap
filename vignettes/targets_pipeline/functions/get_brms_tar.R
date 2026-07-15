#' Read one gene's brms fit from a NEW_CELL_TYPE targets store.
#'
#' @param cell_type Cell-type label (matched via \code{make.names}).
#' @param gene_id Ensembl gene id (e.g. \code{"ENSG00000134419"}).
#' @param base_path Root containing \code{V*_celltype/_targets}. Defaults to
#'   \code{NEW_CELL_TYPE_ROOT}.
#' @return A qs-deserialized target object (typically a tibble with \code{brms_fit}).
get_brms_tar <- function(
    cell_type,
    gene_id,
    base_path = Sys.getenv("NEW_CELL_TYPE_ROOT", unset = "")
) {
  if (!nzchar(base_path)) {
    stop(
      "Set base_path or NEW_CELL_TYPE_ROOT to the NEW_CELL_TYPE directory.",
      call. = FALSE
    )
  }
  base_path <- normalizePath(base_path, winslash = "/", mustWork = TRUE)
  # ensure trailing slash semantics for glue paths below
  if (!grepl("/$", base_path)) {
    base_path <- paste0(base_path, "/")
  }

  cli::cli_alert_info("Locating cell type `{cell_type}` under {base_path} ...")

  ct_list <- list.dirs(base_path, recursive = FALSE) %>% basename()
  ct_list <- ct_list[ct_list %>% stringr::str_starts("V[0-9]_")]
  ct_name <- stringr::str_remove(ct_list, "V[0-9]_")

  if (!make.names(cell_type) %in% ct_name) {
    cli::cli_alert_danger("Cell type '{cell_type}' not found.")
    stop("Cell type not found.", call. = FALSE)
  }

  ct_list <- ct_list[ct_list %>% stringr::str_detect(make.names(cell_type))]
  version_no <- ct_list %>%
    stringr::str_extract("^V\\d+") %>%
    stringr::str_remove("^V") %>%
    as.integer() %>%
    max(na.rm = TRUE)

  targets_path <- glue::glue(
    "{base_path}V{version_no}_{make.names(cell_type)}/_targets"
  )

  cli::cli_alert_info("Using version V{version_no} at {targets_path}")

  mapping_table <-
    readr::read_delim(glue::glue("{targets_path}/meta/meta"), delim = "|") %>%
    dplyr::filter(name %>% stringr::str_starts("estimates_chunk")) %>%
    dplyr::filter(type == "branch" & is.na(error)) %>%
    dplyr::mutate(
      gene = warnings %>%
        stringr::str_extract("(?<=Gene:___)ENSG\\d+(?=___)")
    )

  if (!gene_id %in% mapping_table$gene) {
    cli::cli_alert_danger("Gene ID '{gene_id}' not found in mapping.")
    stop("Gene ID not found.", call. = FALSE)
  }

  target_name <- mapping_table %>%
    dplyr::filter(gene == gene_id) %>%
    dplyr::pull(name)

  if (length(target_name) != 1) {
    cli::cli_alert_danger(
      "Expected 1 target for gene '{gene_id}', found {length(target_name)}."
    )
    stop("Ambiguous or missing target.", call. = FALSE)
  }

  brms_path <- glue::glue("{targets_path}/objects/{target_name}")

  if (!file.exists(brms_path)) {
    cli::cli_alert_danger("Target file not found at {brms_path}")
    stop("Target file missing.", call. = FALSE)
  }

  cli::cli_alert_success("Successfully loaded target for gene '{gene_id}'.")

  qs2::qs_read(brms_path)
}
