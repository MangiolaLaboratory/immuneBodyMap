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

#' Options shared by every multi-select picker in the sidebar.
picker_opts <- function(...) {
    do.call(pickerOptions, utils::modifyList(
        list(
            actionsBox = TRUE,
            size = 10,
            selectedTextFormat = "count > 2"
        ),
        list(...)
    ))
}

#' Which decades of life to pool, as a filter on the other tabs.
#'
#' This replaced a sliderTextInput: ten tick labels under a 290px sidebar
#' collide (the CSS already had to drop them to 9px for six), and a slider can
#' only ever select one bin, where pooling neighbouring decades is the only way
#' some tissue groups reach a usable sample count.
age_picker_input <- function(id) {
    pickerInput(
        id,
        with_help(
            "Decade of life",
            "Samples are binned by decade of life. Select several to pool them;
             the first and last decades carry very few samples on their own."
        ),
        choices = AGE_DECADES,
        multiple = TRUE,
        selected = AGE_DECADES,
        options = picker_opts()
    )
}

ethnicity_picker_input <- function(id) {
    pickerInput(
        id,
        with_help(
            "Ethnicity",
            "Known or imputed ethnicity group."
        ),
        choices = ETHNICITY_PANELS,
        multiple = TRUE,
        selected = ETHNICITY_PANELS,
        options = picker_opts(liveSearch = TRUE)
    )
}

#' Whether samples whose ethnicity was imputed count towards the panels.
#'
#' A quarter of the atlas carries an imputed label, so this is not a fringe
#' setting: switching it off is a visibly different dataset. Default on, and
#' the per-sample flag is a filterable column on the Full data tab so the split
#' can be audited.
imputed_input <- function(id) {
    prettyCheckbox(
        id,
        with_help(
            "Include imputed ethnicity",
            "About 26% of samples have an ethnicity imputed rather than reported.
             Excluding them shrinks every panel, most of all the smaller groups."
        ),
        value = TRUE,
        status = "primary",
        shape = "curve"
    )
}

#' Which decades get a body map of their own on the "By age" tab.
#'
#' Ten panels in a row would leave each about half its current width, and the
#' thin decades would mostly be blank, so the tab opens on the well-covered ones
#' and the rest stay one click away.
decade_panels_input <- function(id) {
    pickerInput(
        id,
        with_help(
            "Decades to show",
            "Each selected decade gets its own body map. Panels share one colour
             scale, which rescales to whatever is shown."
        ),
        choices = AGE_DECADES,
        multiple = TRUE,
        selected = DEFAULT_DECADES,
        options = picker_opts()
    )
}

ethnicity_panels_input <- function(id) {
    pickerInput(
        id,
        with_help(
            "Groups to show",
            "Each selected group gets its own body map. Panels share one colour
             scale, which rescales to whatever is shown."
        ),
        choices = ETHNICITY_PANELS,
        multiple = TRUE,
        selected = DEFAULT_ETHNICITY_PANELS,
        options = picker_opts(liveSearch = TRUE)
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
        options = picker_opts(liveSearch = TRUE)
    )
}

#' Distinguish "nothing selected" from "no filter applied".
#'
#' pickerInput reports NULL when every option is deselected, but NULL means
#' "don't filter" to make_plot_data, so hitting "Deselect all" used to reveal
#' every tissue group rather than hiding them.
none_if_null <- function(x) {
    if (is.null(x)) character(0) else x
}

#' Sex, which doubles as the choice of body outline.
#'
#' gganatogram draws a male and a female body and nothing else, so the third
#' recorded value ("unknown", 136 samples) is not offered here. Those samples
#' are still reachable on the Full data and plot tabs.
sex_input <- function(id) {
    radioGroupButtons(
        id,
        "Sex",
        choices = SEXES,
        justified = TRUE,
        size = "sm"
    )
}
