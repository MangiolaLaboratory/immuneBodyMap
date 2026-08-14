# Shared data-selection controls.
#
# Each comparison tab needs a different subset of these, but the ones they share
# were previously duplicated verbatim. Choices come from constants.R or from the
# data's factor levels, so they are sorted rather than left in data order.

#' Label with a question-mark tooltip carrying the explanation.
with_help <- function(label, help) {
    tagList(
        label,
        tooltip(
            icon("circle-question", class = "help-icon"),
            help,
            placement = "right"
        )
    )
}

cell_type_input <- function(id, choices) {
    selectInput(
        id,
        with_help(
            "Cell type",
            "Immune cell type, harmonised across the atlas by ensemble annotation."
        ),
        choices = choices
    )
}

ethnicity_input <- function(id) {
    selectInput(
        id,
        "Ethnicity",
        choices = ETHNICITIES,
        selected = DEFAULT_ETHNICITY
    )
}

age_select_input <- function(id) {
    selectInput(
        id,
        "Age",
        choices = AGE_BINS,
        selected = DEFAULT_AGE_BIN
    )
}

age_slider_input <- function(id) {
    sliderTextInput(
        id,
        with_help(
            "Age",
            "Age bins have sex-specific boundaries, taken from the developmental
             thresholds used in the paper."
        ),
        choices = AGE_BINS,
        grid = TRUE,
        selected = DEFAULT_AGE_BIN
    )
}

tissue_groups_input <- function(id, choices) {
    pickerInput(
        id,
        with_help(
            "Tissue groups",
            "Each tissue group is drawn as one or more anatomical organs.
             Deselect groups to hide them from the body map."
        ),
        choices = choices,
        multiple = TRUE,
        selected = choices,
        options = pickerOptions(
            actionsBox = TRUE,
            size = 10,
            liveSearch = TRUE,
            selectedTextFormat = "count > 3"
        )
    )
}

#' Distinguish "nothing selected" from "no filter applied".
#'
#' pickerInput reports NULL when every option is deselected, but NULL means
#' "don't filter" to make_plot_data -- so hitting "Deselect all" used to reveal
#' every tissue group rather than hiding them.
none_if_null <- function(x) {
    if (is.null(x)) character(0) else x
}

sex_input <- function(id) {
    radioGroupButtons(
        id,
        "Sex",
        choices = SEXES,
        justified = TRUE,
        size = "sm"
    )
}
