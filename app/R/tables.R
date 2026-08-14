# DT table construction. Sixteen near-identical datatable() blocks collapsed
# into two functions, which also removes the styling drift between tabs.

#' Display names for a frame's columns, falling back to the raw name.
pretty_colnames <- function(dat) {
    nm <- names(dat)
    idx <- match(nm, names(PRETTY_COLNAMES))
    unname(ifelse(is.na(idx), nm, PRETTY_COLNAMES[idx]))
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

#' Per-organ proportion table for one anatogram panel.
#'
#' Organs with no value are dropped: they are an artefact of the tissue-group to
#' organ expansion rather than a result, and the tab-1 tables used to show them
#' while the tab-2 and tab-3 tables did not.
props_table <- function(dat, caption = NULL) {
    dat <- dat[!is.na(dat$value), , drop = FALSE]

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
        formatRound(c("value", "CI_lower", "CI_upper"), 4)
}

#' The complete proportion dataset.
#'
#' Built once at startup rather than per session -- it is the same 18,000-row
#' widget for everyone and nothing about it is reactive.
full_data_table <- function(dat) {
    datatable(
        dat,
        rownames = FALSE,
        colnames = pretty_colnames(dat),
        filter = "top",
        extensions = "Buttons",
        options = dt_options(pageLength = 15)
    ) |>
        formatRound(c("proportion_mean", "proportion_lower", "proportion_upper"), 4) |>
        formatRound(c("age_days", "age_days_scaled"), 2)
}
