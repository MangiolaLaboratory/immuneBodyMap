# Loading, reshaping and summarising the cell-proportion data.

#' Read the proportion data and put its columns in a usable state.
#'
#' Three things happen here beyond tidying.
#'
#' Renaming: the file ships "L3", "assay_groups___altered" and
#' "dataset_id___altered". The VizModules plot controls display raw column
#' names, but these are ugly. They are renamed here.
#'
#' Ordering factor levels for proper display in the UI.
#'
#' Ethnicity: one 16-level column carries two separate facts, the group and
#' whether that group was imputed. They are split so each can be filtered on its
#' own, so drop the raw column.
load_prop_data <- function(path = DATA_FILE) {
    if (!file.exists(path)) {
        stop(
            "Could not find '", path, "' in '", getwd(), "'.\n",
            "The app expects its own directory as the working directory; ",
            "run it with shiny::runApp(\"app\") from the repository root.",
            call. = FALSE
        )
    }

    d <- as.data.frame(readRDS(path))

    names(d)[names(d) == "L3"] <- "cell_type"
    names(d)[names(d) == "assay_groups___altered"] <- "assay"
    names(d)[names(d) == "dataset_id___altered"] <- "dataset_id"

    d$cell_type <- factor(as.character(d$cell_type))
    d$tissue_groups <- factor(as.character(d$tissue_groups))
    d$assay <- factor(as.character(d$assay))
    d$sex <- factor(as.character(d$sex), levels = SEXES_ALL)

    # Re-level chronologically and relabel to the readable form in one step.
    d$age_decade <- factor(
        as.character(d$age_decade),
        levels = AGE_DECADE_LEVELS,
        labels = AGE_DECADES
    )

    raw_ethnicity <- as.character(d$ethnicity_groups_imputed)
    d$ethnicity_imputed <- grepl("_imputed$", raw_ethnicity)
    d$ethnicity <- factor(sub("_imputed$", "", raw_ethnicity), levels = ETHNICITIES)
    d$ethnicity_groups_imputed <- NULL

    # Re-levelling turns an unrecognised value into NA. Fail loudly at startup
    # rather than shipping panels that are quietly empty.
    for (col in c("sex", "ethnicity", "age_decade")) {
        if (anyNA(d[[col]])) {
            stop(
                "Column '", col, "' contains values not listed in R/constants.R. ",
                "Update the constants to match the data.",
                call. = FALSE
            )
        }
    }

    d[, c(
        "sample_id", "tissue_groups", "cell_type", "sex", "age_decade",
        "ethnicity", "ethnicity_imputed", "assay", "dataset_id",
        "proportion_mean", "proportion_lower", "proportion_upper"
    )]
}

#' Fail at startup if a tissue group can neither be drawn nor is known to be
#' undrawable.
#'
#' The companion to the factor validation in load_prop_data(): that one catches
#' a new sex or ethnicity, this one catches a new tissue group. Without it a
#' group the organ map has never heard of just renders as nothing.
check_tissue_coverage <- function(d) {
    known <- unique(c(
        names(MALE_ORGAN_MAP), names(FEMALE_ORGAN_MAP), ORGAN_LESS_GROUPS
    ))
    missing <- setdiff(levels(d$tissue_groups), known)

    if (length(missing)) {
        stop(
            "Tissue groups with no entry in R/organ_maps.R: ",
            paste(missing, collapse = ", "), ".\n",
            "Add them to BASE_ORGAN_MAP, to one of the sex-specific maps, or ",
            "to ORGAN_LESS_GROUPS if gganatogram cannot draw them.",
            call. = FALSE
        )
    }

    invisible(d)
}

#' max() that cannot return -Inf.
#'
#' The colour scale is built with limits = c(0, max_val). When every value in
#' the current selection is NA, which happens with an entirely deselected
#' tissue picker or a filter combination with no samples, max(na.rm = TRUE)
#' returns -Inf and the scale errors out. Fall back to a finite bound instead.
safe_max <- function(..., fallback = 1) {
    vals <- unlist(list(...), use.names = FALSE)
    vals <- vals[is.finite(vals)]
    if (!length(vals)) fallback else max(vals)
}

#' The values a panel actually draws.
#'
#' The shared colour bar describes the body maps, so it has to be scaled to what
#' is visible. Rows for ORGAN_LESS_GROUPS carry real values but no organ, and
#' letting one of them set the top of the scale would dim every drawn organ to
#' match something the user cannot see.
drawn_values <- function(pdata) {
    if (is.null(pdata) || !nrow(pdata)) {
        return(numeric(0))
    }
    pdata$value[!is.na(pdata$organ)]
}

#' Summarise a filtered slice to one row per tissue group.
#'
#' New data made this necessary. The file is one row per
#' real sample rather than the balanced posterior grid originally provide, 
#' so a tissue group can carry anything from 1 to 2,812 samples and a
#' lookup would silently return whichever one happened to come first.
#'
#' @return One row per tissue group present in `x`: the across-sample mean, its
#'   spread (sd), standard error, a normal 95% confidence interval on that mean,
#'   and the sample count it rests on. Note this is an interval for the *mean*,
#'   not the per-sample posterior credible interval carried by proportion_lower
#'   and proportion_upper.
tissue_summary <- function(x) {
    empty <- data.frame(
        tissue_groups = character(0), value = numeric(0), sd = numeric(0),
        se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0),
        n_samples = integer(0), stringsAsFactors = FALSE
    )
    if (is.null(x) || !nrow(x)) {
        return(empty)
    }

    by_group <- split(x$proportion_mean, as.character(x$tissue_groups))
    n <- lengths(by_group)
    means <- vapply(by_group, mean, numeric(1))
    sds <- vapply(
        by_group,
        function(v) if (length(v) > 1L) stats::sd(v) else NA_real_,
        numeric(1)
    )
    ses <- sds / sqrt(n)

    data.frame(
        tissue_groups = names(by_group),
        value = means,
        sd = sds,
        se = ses,
        # A proportion cannot leave [0, 1], and a normal interval on a small
        # sample readily would.
        ci_lower = pmax(0, means - 1.96 * ses),
        ci_upper = pmin(1, means + 1.96 * ses),
        n_samples = as.integer(n),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

#' Build the per-organ frame gganatogram plots.
#'
#' Filters the proportion data down to a single cell type / sex (and optionally
#' decades, ethnicities and imputation status), summarises the matching samples
#' per tissue group, then expands those values out to one row per anatomical
#' organ.
#'
#' `age` and `ethnicity` take vectors, because the controls are multi-select
#' pickers. NULL still means "no filter".
#'
#' @return A data frame of organ, tissue_groups, value, sd, se, ci_lower,
#'   ci_upper and n_samples, ordered so smaller organs draw on top. Organs whose
#'   tissue group has no matching sample carry NA values, which gganatogram
#'   renders as unfilled. Tissue groups with no organ at all appear with
#'   organ = NA so they still reach the tables.
make_plot_data <- function(data, celltype, tissue_groups = NULL, age = NULL,
                           ethnicity = NULL, include_imputed = TRUE,
                           sex = "male") {
    out <- data[data$cell_type == celltype & data$sex == sex, , drop = FALSE]

    if (!is.null(age)) {
        out <- out[out$age_decade %in% age, , drop = FALSE]
    }

    if (!is.null(ethnicity)) {
        out <- out[out$ethnicity %in% ethnicity, , drop = FALSE]
    }

    if (!isTRUE(include_imputed)) {
        out <- out[!out$ethnicity_imputed, , drop = FALSE]
    }

    agg <- tissue_summary(out)

    organ_map <- organ_map_for(sex)

    out_df <- data.frame(
        organ = unlist(organ_map, use.names = FALSE),
        tissue_groups = rep(names(organ_map), lengths(organ_map)),
        stringsAsFactors = FALSE
    )

    out_df <- rbind(out_df, data.frame(
        organ = NA_character_,
        tissue_groups = ORGAN_LESS_GROUPS,
        stringsAsFactors = FALSE
    ))

    out_df$organ <- factor(out_df$organ, levels = organ_levels_for(sex))

    # agg holds one row per tissue group by construction, so this is a lookup
    # rather than a join.
    idx <- match(out_df$tissue_groups, agg$tissue_groups)
    for (col in c("value", "sd", "se", "ci_lower", "ci_upper", "n_samples")) {
        out_df[[col]] <- agg[[col]][idx]
    }

    if (!is.null(tissue_groups)) {
        out_df <- out_df[out_df$tissue_groups %in% tissue_groups, ]
    }

    # order() puts the NA organs last, which suits both the plot (they are
    # dropped before drawing) and the tables.
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

#' Mean proportion per tissue group x cell type, across every sample.
#'
#' The composition dot plot and the heatmap both want the atlas summarised this
#' way. It has to be pre-aggregated, and the granularity has to match exactly
#' what the plot maps: `plotthis::DotPlot` aggregates only when size_by is unset
#' (and then takes first(fill_by), with a warning). With size_by supplied it
#' plots the rows as given, so two rows landing in one x/y/facet cell overplot
#' and whichever is drawn last wins. That is the same arbitrary single-sample
#' lookup make_plot_data() exists to avoid.
#'
#' @param by Extra grouping columns, e.g. "sex". The dot plot tab uses this to
#'   summarise at exactly the granularity it is about to facet on, so the two
#'   cannot drift apart.
#' @return One row per combination of tissue_groups, cell_type and `by`, with
#'   mean_proportion and n_samples. Grouping columns keep their factor levels,
#'   so facets and axes stay in the data's order rather than alphabetical.
tissue_celltype_summary <- function(d, by = character(0)) {
    group_cols <- unique(c("tissue_groups", "cell_type", by))

    # A carriage return cannot occur in any of these names, so it is a safe
    # compound key to split on and then take apart again.
    key <- do.call(
        paste,
        c(lapply(group_cols, function(col) as.character(d[[col]])), sep = "\r")
    )
    by_cell <- split(d$proportion_mean, key)
    parts <- do.call(rbind, strsplit(names(by_cell), "\r", fixed = TRUE))

    out <- as.data.frame(parts, stringsAsFactors = FALSE)
    names(out) <- group_cols

    # Restore the factor levels the compound key threw away. Without this an
    # age_decade facet would run 0-9, 10-19, 100-109 ... alphabetically, and sex
    # would come back female-first rather than in SEXES order.
    for (col in group_cols) {
        if (is.factor(d[[col]])) {
            out[[col]] <- factor(out[[col]], levels = levels(d[[col]]))
        }
    }

    out$mean_proportion <- unname(vapply(by_cell, mean, numeric(1)))
    out$n_samples <- unname(as.integer(lengths(by_cell)))
    row.names(out) <- NULL

    out[do.call(order, unname(out[group_cols])), ]
}

#' Reshape tissue_celltype_summary() to the wide frame ComplexHeatmap wants.
#'
#' @return One row per cell type: an id column plus one numeric column per
#'   tissue group. Combinations with no samples stay NA, which the module draws
#'   in its na_col rather than as zero.
celltype_tissue_matrix <- function(summary_df) {
    cells <- sort(unique(summary_df$cell_type))
    tissues <- sort(unique(summary_df$tissue_groups))

    m <- matrix(
        NA_real_,
        nrow = length(cells), ncol = length(tissues),
        dimnames = list(NULL, tissues)
    )
    m[cbind(
        match(summary_df$cell_type, cells),
        match(summary_df$tissue_groups, tissues)
    )] <- summary_df$mean_proportion

    data.frame(
        cell_type = cells,
        m,
        check.names = FALSE,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}
