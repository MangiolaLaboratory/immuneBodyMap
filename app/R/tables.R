# DT table construction. 

#' Display names for a frame's columns, falling back to the raw name.
pretty_colnames <- function(dat) {
    nm <- names(dat)
    idx <- match(nm, names(PRETTY_COLNAMES))
    unname(ifelse(is.na(idx), nm, PRETTY_COLNAMES[idx]))
}

#' formatRound() over whatever subset of `cols` the widget actually has..
#'
#' formatRound() stops with "You specified the columns: ..." when handed a name
#' that is not in the data.
format_round_if <- function(widget, cols, digits) {
    cols <- intersect(cols, names(widget$x$data))
    if (length(cols)) formatRound(widget, cols, digits) else widget
}

#' Options shared by every table in the app.
dt_options <- function(...) {
    utils::modifyList(
        list(
            search = list(regex = TRUE),
            pageLength = 10,
            dom = "Blfrtip",
            buttons = c("copy", "csv", "excel", "pdf", "print"),
            scrollX = TRUE,
            autoWidth = FALSE
        ),
        list(...)
    )
}

# se is dropped from the panel tables since ci_lower and ci_upper cover that
PANEL_TABLE_COLS <- c(
    "organ", "tissue_groups", "value", "ci_lower", "ci_upper", "sd", "n_samples"
)

#' Per-organ proportion table for one anatogram panel.
#'
#' Organs with no value are dropped.
props_table <- function(dat, caption = NULL) {
    dat <- dat[!is.na(dat$value), , drop = FALSE]
    dat <- dat[, intersect(PANEL_TABLE_COLS, names(dat)), drop = FALSE]

    datatable(
        dat,
        rownames = FALSE,
        colnames = pretty_colnames(dat),
        filter = "top",
        caption = caption,
        extensions = "Buttons",
        options = dt_options(
            columnDefs = list(list(width = "30%", targets = 1))
        )
    ) |>
        format_round_if(c("value", "ci_lower", "ci_upper", "sd"), 4)
}

#' The complete proportion dataset.
#'
#' Built once at startup rather than per session, because it is the same
#' 116,558-row widget for everyone and nothing about it is reactive. renderDT()
#' pages it server-side, so the row count is not a transport problem.
#'
#' The columns here are per-sample posterior quantities, not the across-sample
#' summaries the panel tables show; PRETTY_COLNAMES labels them CrI rather than
#' CI to keep the two apart.
full_data_table <- function(dat) {
    datatable(
        dat,
        rownames = FALSE,
        colnames = pretty_colnames(dat),
        filter = "top",
        extensions = "Buttons",
        options = dt_options(pageLength = 15)
    ) |>
        format_round_if(
            c("proportion_mean", "proportion_lower", "proportion_upper"), 4
        )
}
