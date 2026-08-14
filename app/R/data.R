# Loading and reshaping the cell-proportion data.

#' Read the proportion data and put its categorical columns in a usable state.
#'
#' Factors matter for more than tidiness here: DT renders a dropdown filter for a
#' factor column and a free-text box for a character one, and the age levels
#' stored in the file are alphabetical rather than chronological.
load_prop_data <- function(path = "cell_proportion_for_shiny_app.rds") {
    if (!file.exists(path)) {
        stop(
            "Could not find '", path, "' in '", getwd(), "'.\n",
            "The app expects its own directory as the working directory -- ",
            "run it with shiny::runApp(\"app\") from the repository root.",
            call. = FALSE
        )
    }

    d <- readRDS(path)

    d$cell_type_unified_ensemble <- factor(d$cell_type_unified_ensemble)
    d$tissue_groups <- factor(d$tissue_groups)
    d$sex <- factor(as.character(d$sex), levels = SEXES)
    d$ethnicity_groups <- factor(as.character(d$ethnicity_groups), levels = ETHNICITIES)
    d$age_bin_sex_specific <- factor(as.character(d$age_bin_sex_specific), levels = AGE_BINS)

    # Re-levelling turns an unrecognised value into NA. Fail loudly at startup
    # rather than shipping panels that are quietly empty.
    for (col in c("sex", "ethnicity_groups", "age_bin_sex_specific")) {
        if (anyNA(d[[col]])) {
            stop(
                "Column '", col, "' contains values not listed in R/constants.R. ",
                "Update the constants to match the data.",
                call. = FALSE
            )
        }
    }

    d
}

#' max() that cannot return -Inf.
#'
#' The colour scale is built with limits = c(0, max_val). When every value in
#' the current selection is NA -- an entirely deselected tissue picker, or a
#' filter combination with no samples -- max(na.rm = TRUE) returns -Inf and the
#' scale errors out. Fall back to a finite upper bound instead.
safe_max <- function(..., fallback = 1) {
    vals <- unlist(list(...), use.names = FALSE)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) fallback else max(vals)
}

#' Build the per-organ frame gganatogram plots.
#'
#' Filters the proportion data down to a single cell type / sex (and optionally
#' age bin and ethnicity), then expands the tissue-group values out to one row
#' per anatomical organ.
#'
#' @return A data frame of organ, tissue_groups, value, CI_lower, CI_upper,
#'   ordered so smaller organs draw on top. Organs whose tissue group has no
#'   matching row carry NA values, which gganatogram renders as unfilled.
make_plot_data <- function(data, celltype, tissue_groups = NULL, age = NULL,
                           ethnicity = NULL, sex = "male") {
    out <- data

    if (!is.null(age)) {
        out <- out[out$age_bin_sex_specific == age, ]
    }

    if (!is.null(ethnicity)) {
        out <- out[out$ethnicity_groups == ethnicity, ]
    }

    out <- out[out$cell_type_unified_ensemble == celltype & out$sex == sex, ]

    organ_map <- organ_map_for(sex)

    out_df <- data.frame(
        organ = unlist(organ_map, use.names = FALSE),
        tissue_groups = rep(names(organ_map), lengths(organ_map)),
        stringsAsFactors = FALSE
    )

    out_df$organ <- factor(out_df$organ, levels = organ_levels_for(sex))

    # The filters above leave at most one row per tissue group, so this is a
    # lookup rather than a join -- match() would silently take the first row if
    # that ever stopped holding.
    idx <- match(out_df$tissue_groups, as.character(out$tissue_groups))
    out_df$value <- out$proportion_mean[idx]
    out_df$CI_lower <- out$proportion_lower[idx]
    out_df$CI_upper <- out$proportion_upper[idx]

    if (!is.null(tissue_groups)) {
        out_df <- out_df[out_df$tissue_groups %in% tissue_groups, ]
    }

    out_df[order(out_df$organ), ]
}

#' make_plot_data() applied across a set of panels sharing the same filters.
#'
#' `by` names the argument that varies between panels ("age" or "ethnicity");
#' `values` supplies one level per panel. Replaces what used to be six and seven
#' copy-pasted calls.
make_panel_data <- function(data, values, by, ...) {
    panels <- lapply(values, function(v) {
        args <- c(list(data = data), list(...))
        args[[by]] <- v
        do.call(make_plot_data, args)
    })
    stats::setNames(panels, slugify(values))
}
