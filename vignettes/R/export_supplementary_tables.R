# Shared helpers for exporting plot statistics to Excel supplementary tables.
# Tabular data are taken from the `data` slot of ggplot objects (or `@data` for
# S4 plot objects such as tidyHeatmap / ComplexHeatmap wrappers).

suppressPackageStartupMessages({
  if (!requireNamespace("writexl", quietly = TRUE)) {
    install.packages("writexl", repos = "https://cloud.r-project.org")
  }
})

#' @param name Candidate Excel sheet name.
#' @param used Names already assigned in the current workbook.
#' @return A valid, unique Excel sheet name (max 31 characters).
sanitize_sheet_name <- function(name, used = character()) {
  nm <- gsub("Plot__|Plot___", "", name)
  nm <- gsub("^Plot_", "", nm)
  nm <- gsub("[\\\\/:?*\\[\\]]", "_", nm)
  nm <- gsub("_+", "_", nm)
  nm <- substr(nm, 1, 31)

  if (!nzchar(nm)) {
    nm <- "sheet"
  }

  if (nm %in% used) {
    i <- 2L
    while (paste0(substr(nm, 1, 28), sprintf("_%d", i)) %in% used) {
      i <- i + 1L
    }
    nm <- paste0(substr(nm, 1, 28), sprintf("_%d", i))
  }

  nm
}

#' Coerce list / complex columns so writexl can export them.
sanitize_for_excel <- function(df) {
  df <- as.data.frame(df)

  for (col in names(df)) {
    x <- df[[col]]

    if (is.list(x) && !inherits(x, "data.frame")) {
      df[[col]] <- vapply(
        x,
        function(el) {
          if (is.null(el) || length(el) == 0) {
            return(NA_character_)
          }
          if (length(el) == 1L && is.atomic(el)) {
            return(as.character(el))
          }
          paste(as.character(unlist(el)), collapse = "; ")
        },
        FUN.VALUE = character(1)
      )
    } else if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
      df[[col]] <- as.character(x)
    } else if (is.factor(x)) {
      df[[col]] <- as.character(x)
    }
  }

  df
}

#' Extract the tabular data backing a plot object.
extract_plot_data <- function(obj, name = "plot") {
  if (is.data.frame(obj)) {
    return(as.data.frame(obj))
  }

  if (inherits(obj, "ggplot")) {
    d <- obj$data
    if (!is.null(d) && nrow(d) > 0) {
      return(as.data.frame(d))
    }
  }

  if (isS4(obj) && "data" %in% methods::slotNames(obj)) {
    d <- obj@data
    if (!is.null(d) && nrow(d) > 0) {
      return(as.data.frame(d))
    }
  }

  if (is.list(obj) && !is.null(obj$data) && is.data.frame(obj$data)) {
    return(as.data.frame(obj$data))
  }

  warning("Could not extract data from: ", name, call. = FALSE)
  NULL
}

#' Stack several related statistical tables for panel-aware export.
#'
#' A panel can need several incompatible tables (for example plotted
#' observations, boxplot quantiles and hypothesis tests). This helper labels
#' each block with `record_type`. The panel workbook writer then separates those
#' blocks into worksheets sharing the same figure-panel prefix.
#'
#' @param ... Named data frames.  Names become values in `record_type`.
#' @param .descriptions Optional named character vector describing the blocks.
#' @return A named `panel_records` list for the panel workbook writer.
bind_panel_records <- function(..., .descriptions = NULL) {
  sections <- list(...)

  if (length(sections) == 1L && is.list(sections[[1]]) &&
      !is.data.frame(sections[[1]])) {
    sections <- sections[[1]]
  }

  keep <- !vapply(sections, is.null, logical(1))
  sections <- sections[keep]

  if (length(sections) == 0L) {
    return(data.frame())
  }
  if (is.null(names(sections)) || any(names(sections) == "")) {
    stop("Every panel record block must have a name.", call. = FALSE)
  }

  blocks <- lapply(names(sections), function(section_name) {
    d <- sanitize_for_excel(sections[[section_name]])
    description <- NA_character_
    if (!is.null(.descriptions) && section_name %in% names(.descriptions)) {
      description <- unname(.descriptions[[section_name]])
    }

    d$record_description <- description
    d[c("record_description", setdiff(names(d), "record_description"))]
  })
  names(blocks) <- names(sections)
  structure(blocks, class = c("panel_records", "list"))
}

#' Distribution statistics used alongside boxplots and cohort dot plots.
#'
#' @param data Input data frame.
#' @param group_cols Character vector of grouping columns.
#' @param value_col Numeric column to summarise.
#' @param value_label Prefix used for the statistic columns.
summarise_distribution <- function(
    data,
    group_cols,
    value_col,
    value_label = value_col) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("summarise_distribution() requires dplyr.", call. = FALSE)
  }

  missing_cols <- setdiff(c(group_cols, value_col), names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing distribution column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  safe_quantile <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    unname(stats::quantile(x, probs = probability, na.rm = TRUE, names = FALSE))
  }
  safe_stat <- function(x, fn) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    fn(x)
  }

  d <- as.data.frame(data)
  d$.supp_value <- suppressWarnings(as.numeric(d[[value_col]]))

  out <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_non_missing = sum(is.finite(.data$.supp_value)),
      mean = safe_stat(.data$.supp_value, base::mean),
      sd = safe_stat(.data$.supp_value, stats::sd),
      min = safe_stat(.data$.supp_value, base::min),
      q1 = safe_quantile(.data$.supp_value, 0.25),
      median = safe_quantile(.data$.supp_value, 0.50),
      q3 = safe_quantile(.data$.supp_value, 0.75),
      max = safe_stat(.data$.supp_value, base::max),
      .groups = "drop"
    ) |>
    as.data.frame()

  stat_cols <- c("mean", "sd", "min", "q1", "median", "q3", "max")
  names(out)[match(stat_cols, names(out))] <- paste0(value_label, "_", stat_cols)
  out
}

#' Per-group linear-model and prediction-error statistics.
summarise_linear_models <- function(data, group_cols, x_col, y_col) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("summarise_linear_models() requires dplyr.", call. = FALSE)
  }

  d <- as.data.frame(data)
  d$.supp_x <- suppressWarnings(as.numeric(d[[x_col]]))
  d$.supp_y <- suppressWarnings(as.numeric(d[[y_col]]))

  d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(function(.x, .y) {
      z <- .x[is.finite(.x$.supp_x) & is.finite(.x$.supp_y), , drop = FALSE]
      if (nrow(z) < 2L || length(unique(z$.supp_x)) < 2L) {
        return(data.frame(
          n = nrow(z), intercept = NA_real_, intercept_se = NA_real_,
          slope = NA_real_, slope_se = NA_real_, slope_p_value = NA_real_,
          r_squared = NA_real_, adjusted_r_squared = NA_real_,
          correlation = NA_real_, mae_identity = NA_real_, rmse_identity = NA_real_
        ))
      }

      fit <- stats::lm(.supp_y ~ .supp_x, data = z)
      fit_summary <- summary(fit)
      coefficients <- fit_summary$coefficients
      identity_error <- z$.supp_y - z$.supp_x

      data.frame(
        n = nrow(z),
        intercept = unname(coefficients[1, "Estimate"]),
        intercept_se = unname(coefficients[1, "Std. Error"]),
        slope = unname(coefficients[2, "Estimate"]),
        slope_se = unname(coefficients[2, "Std. Error"]),
        slope_p_value = unname(coefficients[2, "Pr(>|t|)"]),
        r_squared = unname(fit_summary$r.squared),
        adjusted_r_squared = unname(fit_summary$adj.r.squared),
        correlation = stats::cor(z$.supp_x, z$.supp_y),
        mae_identity = mean(abs(identity_error)),
        rmse_identity = sqrt(mean(identity_error ^ 2))
      )
    }) |>
    dplyr::ungroup() |>
    as.data.frame()
}

#' Convert named plot objects / data frames into a list of Excel sheets.
plots_to_sheets <- function(plots = list(), extra_sheets = list(), prefix = "") {
  sheets <- list()
  used <- character()

  all_items <- c(plots, extra_sheets)
  if (is.null(names(all_items)) || any(names(all_items) == "")) {
    names(all_items) <- paste0("sheet_", seq_along(all_items))
  }

  for (nm in names(all_items)) {
    d <- extract_plot_data(all_items[[nm]], nm)
    if (is.null(d)) {
      next
    }

    sheet_nm <- sanitize_sheet_name(paste0(prefix, nm), used)
    used <- c(used, sheet_nm)
    sheets[[sheet_nm]] <- sanitize_for_excel(d)
  }

  sheets
}

#' Write one figure-level supplementary workbook.
export_supplementary_workbook <- function(
    plots = list(),
    extra_sheets = list(),
    output_path,
    prefix = "",
    also_csv = FALSE,
    csv_dir = NULL) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  sheets <- plots_to_sheets(plots = plots, extra_sheets = extra_sheets, prefix = prefix)
  if (length(sheets) == 0) {
    warning("No sheets to write for: ", output_path, call. = FALSE)
    return(invisible(NULL))
  }

  writexl::write_xlsx(sheets, path = output_path)

  if (isTRUE(also_csv) && !is.null(csv_dir)) {
    dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("readr", quietly = TRUE)) {
      for (nm in names(sheets)) {
        readr::write_csv(sheets[[nm]], file.path(csv_dir, paste0(nm, ".csv")))
      }
    }
  }

  message("Wrote ", length(sheets), " sheet(s) to ", output_path)
  invisible(sheets)
}

#' Combine multiple sheet lists into one master workbook.
build_master_workbook <- function(
    sheet_lists,
    output_path,
    existing_csvs = list()) {
  master <- list()
  used <- character()

  for (item in sheet_lists) {
    for (nm in names(item)) {
      sheet_nm <- sanitize_sheet_name(nm, used)
      used <- c(used, sheet_nm)
      master[[sheet_nm]] <- item[[nm]]
    }
  }

  for (nm in names(existing_csvs)) {
    path <- existing_csvs[[nm]]
    if (!file.exists(path)) {
      warning("Missing CSV for master workbook: ", path, call. = FALSE)
      next
    }
    if (requireNamespace("readr", quietly = TRUE)) {
      d <- readr::read_csv(path, show_col_types = FALSE)
    } else {
      d <- utils::read.csv(path, check.names = FALSE)
    }
    sheet_nm <- sanitize_sheet_name(nm, used)
    used <- c(used, sheet_nm)
    master[[sheet_nm]] <- sanitize_for_excel(as.data.frame(d))
  }

  if (length(master) == 0) {
    warning("Master workbook has no sheets: ", output_path, call. = FALSE)
    return(invisible(NULL))
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(master, path = output_path)
  message("Wrote master workbook (", length(master), " sheets) to ", output_path)
  invisible(master)
}

#' Resolve the Age_Clock project root from an optional start path.
age_clock_root <- function(start = getwd()) {
  configured <- Sys.getenv("AGE_CLOCK_ROOT", unset = "")
  home_candidate <- file.path(
    path.expand("~"), "lab", "chen", "HPC_posterior", "Age_Clock"
  )

  candidates <- unique(c(
    configured,
    home_candidate,
    start,
    file.path(start, ".."),
    file.path(start, "../.."),
    if (requireNamespace("here", quietly = TRUE)) here::here() else NULL
  ))

  for (cand in candidates) {
    if (!nzchar(cand)) next
    if (dir.exists(file.path(cand, "Supplementary_files")) &&
        dir.exists(file.path(cand, "report"))) {
      return(cand)
    }
  }

  start
}
