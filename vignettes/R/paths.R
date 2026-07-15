# Portable path resolution for Nature Aging publication vignettes.
#
# Primary configuration:
#   IMMUNE_HEALTHY_BODY_MAP_DATA  - root for publication data (default: vignettes/data)
#   IMMUNE_HEALTHY_BODY_MAP_ROOT  - package / repo root (default: here::here())
#   CELLNEXUS_METADATA            - optional local parquet override (default: cloud via get_metadata_url)
#   CELLNEXUS_CACHE_DIR           - cache dir for downloaded cellNexus metadata
#   SCCOMP_ESTIMATES_DIR          - directory of fitted sccomp estimate RDS files
#   PSEUDOBULK_H5AD               - harmonised pseudobulk h5ad (Zenodo)
#   AGE_CLOCK_ARCHIVE             - optional mirror of Age_Clock analysis tree (full mode)
#   NEW_PRED_DIR                  - output/input dir for gene-ageing predictions (Fig4 rebuild)
#   NEW_CELL_TYPE_ROOT            - NEW_CELL_TYPE targets root (Fig1 adjust lists / Fig4 brms)
#   TAR_PIPELINE_ROOT             - optional override of vignettes/targets_pipeline

publication_repo_root <- function() {
  env <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  # When this file is sourced, vignettes/R/paths.R → repo root is ../..
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)$ofile
    if (!is.null(f) && grepl("[/\\\\]vignettes[/\\\\]R[/\\\\]paths\\.R$", f)) {
      return(normalizePath(
        file.path(dirname(f), "..", ".."),
        winslash = "/",
        mustWork = FALSE
      ))
    }
  }
  if (requireNamespace("here", quietly = TRUE)) {
    cand <- normalizePath(here::here(), winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(cand, "vignettes"))) {
      return(cand)
    }
  }
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

publication_data_root <- function() {
  env <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_DATA", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  normalizePath(
    file.path(publication_repo_root(), "vignettes", "data"),
    winslash = "/",
    mustWork = FALSE
  )
}

publication_vignette_root <- function() {
  normalizePath(
    file.path(publication_repo_root(), "vignettes"),
    winslash = "/",
    mustWork = FALSE
  )
}

#' Informative stop when a required file is missing.
require_publication_file <- function(
    path,
    description,
    zenodo_hint = FALSE,
    env_override = NULL
) {
  if (file.exists(path) || dir.exists(path)) {
    return(invisible(normalizePath(path, winslash = "/", mustWork = FALSE)))
  }
  msg <- paste0(
    "Missing required publication file:\n",
    "  description : ", description, "\n",
    "  expected at : ", path, "\n"
  )
  if (!is.null(env_override) && nzchar(env_override)) {
    msg <- paste0(
      msg,
      "  override via: Sys.setenv(", env_override, " = \"<path>\") ",
      "or export ", env_override, "=<path>\n"
    )
  }
  if (isTRUE(zenodo_hint)) {
    msg <- paste0(
      msg,
      "  Pseudobulk: https://doi.org/10.5281/zenodo.15798373 ",
      "(pseudobulk_se.h5ad; md5 88c71c0fd1d6ce2fe15eccdd7b36110f).\n",
      "  Other large companions may still use <ZENODO_DOI_PENDING>.\n",
      "  Place files under IMMUNE_HEALTHY_BODY_MAP_DATA/zenodo_release/ ",
      "or set the env override above, then retry.\n"
    )
  }
  stop(msg, call. = FALSE)
}

path_source_tables <- function(...) {
  file.path(publication_data_root(), "source_tables", ...)
}

path_figure_ready <- function(...) {
  file.path(publication_data_root(), "figure_ready", ...)
}

path_intermediate <- function(...) {
  file.path(publication_data_root(), "intermediate", ...)
}

path_processed <- function(...) {
  file.path(publication_data_root(), "processed", ...)
}

#' Publication copy of the NEW_CELL_TYPE targets pipeline (from TAR_SCRIPTS).
path_targets_pipeline <- function(...) {
  env <- Sys.getenv("TAR_PIPELINE_ROOT", unset = "")
  root <- if (nzchar(env)) {
    normalizePath(env, winslash = "/", mustWork = FALSE)
  } else {
    file.path(publication_vignette_root(), "targets_pipeline")
  }
  file.path(root, ...)
}

path_new_cell_type_root <- function() {
  env <- Sys.getenv("NEW_CELL_TYPE_ROOT", unset = "")
  if (!nzchar(env)) {
    stop(
      "Set NEW_CELL_TYPE_ROOT to the NEW_CELL_TYPE directory ",
      "(V*_celltype/_targets stores; targets-cache exemption).",
      call. = FALSE
    )
  }
  normalizePath(env, winslash = "/", mustWork = FALSE)
}

#' Resolve edit_covariates.R without hard-coded HPC paths.
path_edit_covariates <- function() {
  candidates <- c(
    Sys.getenv("EDIT_COVARIATES_PATH", unset = ""),
    path_targets_pipeline("functions", "edit_covariates.R"),
    file.path(
      publication_repo_root(),
      "rebuttal_CellPress",
      "edit_covariates.R"
    ),
    file.path(
      publication_repo_root(),
      "vignettes",
      "R",
      "edit_covariates.R"
    )
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Could not find edit_covariates.R.\n",
    "Expected rebuttal_CellPress/edit_covariates.R in the repository, ",
    "or set EDIT_COVARIATES_PATH.",
    call. = FALSE
  )
}

path_cellnexus_metadata <- function() {
  env <- Sys.getenv("CELLNEXUS_METADATA", unset = "")
  if (nzchar(env)) {
    return(require_publication_file(
      env,
      "optional local cellNexus metadata parquet override",
      zenodo_hint = FALSE,
      env_override = "CELLNEXUS_METADATA"
    ))
  }
  local <- file.path(
    publication_data_root(),
    "processed",
    "metadata.1.0.12.parquet"
  )
  if (file.exists(local)) {
    return(normalizePath(local, winslash = "/", mustWork = TRUE))
  }
  NULL
}

#' Load cellNexus metadata from the public cloud URL by default.
#' Set CELLNEXUS_METADATA to force a local parquet, or CELLNEXUS_CACHE_DIR for caching.
get_publication_metadata <- function(...) {
  if (!requireNamespace("cellNexus", quietly = TRUE)) {
    stop("Package cellNexus is required for metadata download.", call. = FALSE)
  }
  local <- path_cellnexus_metadata()
  cache_dir <- Sys.getenv(
    "CELLNEXUS_CACHE_DIR",
    unset = file.path(publication_data_root(), "processed", "cellNexus_cache")
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(local)) {
    message("Using local cellNexus metadata: ", local)
    return(cellNexus::get_metadata(
      cloud_metadata = NULL,
      local_metadata = local,
      cache_directory = cache_dir,
      ...
    ))
  }

  url <- cellNexus::get_metadata_url("metadata.1.0.12.parquet")
  message("Downloading / caching cellNexus metadata from: ", url)
  cellNexus::get_metadata(
    cloud_metadata = url,
    cache_directory = cache_dir,
    ...
  )
}

path_pseudobulk_h5ad <- function() {
  env <- Sys.getenv("PSEUDOBULK_H5AD", unset = "")
  if (nzchar(env)) {
    return(require_publication_file(
      env,
      "harmonised pseudobulk SummarizedExperiment (h5ad)",
      zenodo_hint = TRUE,
      env_override = "PSEUDOBULK_H5AD"
    ))
  }
  local <- file.path(
    publication_data_root(),
    "zenodo_release",
    "pseudobulk",
    "pseudobulk_se.h5ad"
  )
  require_publication_file(
    local,
    "harmonised pseudobulk SummarizedExperiment (h5ad)",
    zenodo_hint = TRUE,
    env_override = "PSEUDOBULK_H5AD"
  )
}

path_sccomp_estimates <- function(
    filename = "estimates_age_decade___L3___disease_FALSE___immune_only_TRUE.rds"
) {
  env_dir <- Sys.getenv("SCCOMP_ESTIMATES_DIR", unset = "")
  candidates <- c(
    if (nzchar(env_dir)) file.path(env_dir, filename) else "",
    file.path(
      publication_data_root(),
      "zenodo_release",
      "sccomp_estimates",
      filename
    ),
    file.path(path_intermediate("sccomp_estimates"), filename)
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  require_publication_file(
    candidates[[length(candidates)]],
    paste0("fitted sccomp estimates: ", filename),
    zenodo_hint = TRUE,
    env_override = "SCCOMP_ESTIMATES_DIR"
  )
}

#' Optional Age_Clock analysis mirror for full recomputation.
path_age_clock_archive <- function(...) {
  env <- Sys.getenv("AGE_CLOCK_ARCHIVE", unset = "")
  if (!nzchar(env)) {
    return(NULL)
  }
  file.path(normalizePath(env, winslash = "/", mustWork = FALSE), ...)
}

#' Cache-first loader with rebuild callback.
read_or_rebuild <- function(cache_path, rebuild_fun, description = cache_path) {
  if (file.exists(cache_path)) {
    message("Loading cached intermediate: ", cache_path)
    if (grepl("\\.RData$", cache_path, ignore.case = TRUE)) {
      e <- new.env(parent = emptyenv())
      load(cache_path, envir = e)
      return(as.list(e))
    }
    return(readRDS(cache_path))
  }
  message(
    "Cache missing for ", description, " at ", cache_path,
    " — rebuilding from inputs..."
  )
  result <- rebuild_fun()
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(result) && grepl("\\.rds$", cache_path, ignore.case = TRUE)) {
    saveRDS(result, cache_path)
  }
  result
}

#' Prefer publication figure_ready plot, else Age_Clock archive, else NULL.
resolve_plot_rds <- function(..., archive_rel = NULL) {
  pub <- path_figure_ready(...)
  if (file.exists(pub)) {
    return(pub)
  }
  if (!is.null(archive_rel)) {
    arch <- path_age_clock_archive(archive_rel)
    if (!is.null(arch) && file.exists(arch)) {
      return(arch)
    }
  }
  NULL
}

source_publication_helpers <- function() {
  helper <- file.path(
    publication_vignette_root(),
    "R",
    "export_supplementary_tables.R"
  )
  if (file.exists(helper)) {
    source(helper)
  } else {
    warning("export_supplementary_tables.R not found at ", helper)
  }
}

message(
  "Publication paths ready.\n",
  "  repo root : ", publication_repo_root(), "\n",
  "  data root : ", publication_data_root()
)
