# Building and launching the Immune Proportion Map.
#
# The app is assembled by a function rather than written out at the top level of
# app.r so that everything optional about it is an argument you can see and
# pass. No environment variable or global option decides what gets built.
# app.r is then just a call to immune_map_app().

#' The bslib theme shared by every tab.
#'
#' A function rather than a stored value so sourcing this file does not require
#' bslib to be attached yet.
app_theme <- function() {
    bs_theme(
        version = 5,
        preset = "shiny",
        # A font collection rather than font_google(): no network fetch at startup,
        # and it degrades to the platform UI font.
        base_font = font_collection(
            "Inter", "Segoe UI", "Helvetica Neue", "Arial", "sans-serif"
        ),
        primary = "#2a6f97",
        "border-radius" = "0.5rem"
    )
}


# ------------------------------------------------------------------ ui ----

#' A wrapping grid of anatogram panels.
#'
#' Pulled out of comparison_tab() because the age and ethnicity tabs now build
#' their grid inside renderUI(), from whichever panels the user asked for.
panel_grid <- function(panels, width = 1 / length(panels), fixed = FALSE) {
    layout_column_wrap(
        width = width,
        fixed_width = fixed,
        heights_equal = "row",
        gap = "0.4rem",
        class = "panel-grid",
        !!!panels
    )
}

#' The grid for a set of panel values, or a note if none are selected.
selected_panel_grid <- function(values) {
    if (!length(values)) {
        return(div(
            class = "panel-empty",
            "No panels selected. Choose at least one above to draw a body map."
        ))
    }

    panel_grid(anatogram_panel_row_ui(values))
}

#' A comparison tab: controls on the left, a grid of panels, a shared colour
#' bar, and the numbers underneath.
comparison_tab <- function(title, id, data_controls, panels_ui, table_ui) {
    nav_panel(
        title,
        layout_sidebar(
            sidebar = sidebar(
                width = 290,
                h6("Data selection", class = "sidebar-heading"),
                !!!data_controls,
                hr(),
                h6("Plot aesthetics", class = "sidebar-heading"),
                aesthetics_ui(id)
            ),
            panels_ui,
            plotOutput(paste0(id, "_legend"), height = "100px"),
            table_ui
        )
    )
}

#' The shared card holding a tab's panel table.
panel_table_card <- function(id) {
    card(
        min_height = 560,
        card_header("Panel data"),
        # fillable = FALSE stops the card body from becoming its own scroll
        # container, so a full page of rows is read by scrolling the page.
        card_body(
            fillable = FALSE,
            div(class = "table-wrap", DTOutput(id))
        )
    )
}

#' Build the Immune Proportion Map app.
#'
#' @param extra_tabs Whether to build the VizModules plot tabs (Composition,
#'   Distributions, Ageing trend, Heatmap). FALSE leaves the body maps, the full
#'   data table and the About page, and skips loading VizModules, plotly and the
#'   Bioconductor heatmap packages entirely, which is worth it for a lighter
#'   deployment.
#'   Ignored, with the tabs left out, when VizModules is not installed.
#' @param data_file Proportion data to read; defaults to the file named in
#'   R/constants.R.
#'
#' @return A Shiny app object, so the result can be handed straight to
#'   shiny::runApp() or printed at the end of app.r.
#'
#' @examples
#' \dontrun{
#' immune_map_app()                    # body maps plus the plot tabs
#' immune_map_app(extra_tabs = FALSE)    # body maps only
#' }
immune_map_app <- function(extra_tabs = TRUE, data_file = DATA_FILE) {
    prop_data <- load_prop_data(data_file)

    # Fails loudly if the data has gained a tissue group no organ map knows
    # about, rather than letting it render as nothing on every panel.
    check_tissue_coverage(prop_data)

    cell_types <- levels(prop_data$cell_type)
    tissue_groups <- levels(prop_data$tissue_groups)

    extra_tabs <- isTRUE(extra_tabs) && requireNamespace("VizModules", quietly = TRUE)
    if (extra_tabs) {
        library(VizModules)
    }
    viz <- if (extra_tabs) viz_data(prop_data) else NULL

    # Static widget: the same for every session, so build it once at startup
    # rather than re-serialising 116,000 rows per connection.
    full_data_widget <- full_data_table(prop_data)

    ui <- page_navbar(
        title = tags$span(class = "brand", "Immune Proportion Map"),
        window_title = "Immune Proportion Map",
        id = "nav",
        theme = app_theme(),
        lang = "en",
        fillable = FALSE,
        header = tagList(
            tags$head(
                tags$link(rel = "stylesheet", href = "styles.css"),
                # bslib's own busy indicators: a page-level pulse plus a spinner on
                # any output that is recalculating. No extra package needed.
                useBusyIndicators()
            ),
            # needed for VizModules' dynamic UI features
            if (extra_tabs) shinyjs::useShinyjs()
        ),

        comparison_tab(
            title = "By sex",
            id = "sex",
            data_controls = list(
                cell_type_input("sex_cell_type", cell_types),
                ethnicity_picker_input("sex_ethnicity"),
                imputed_input("sex_imputed"),
                age_picker_input("sex_age"),
                tissue_groups_input("sex_tissue_groups", tissue_groups)
            ),
            panels_ui = panel_grid(
                list(
                    anatogram_panel_ui("male", "Male", show_data = FALSE, height = SEX_CARD_HEIGHT),
                    anatogram_panel_ui("female", "Female", show_data = FALSE, height = SEX_CARD_HEIGHT)
                ),
                width = SEX_CARD_WIDTH,
                fixed = TRUE
            ),
            table_ui = layout_column_wrap(
                width = 1 / 2,
                heights_equal = "row",
                card(
                    card_header("Male proportions"),
                    card_body(div(class = "table-wrap", DTOutput("sex_male_props")))
                ),
                card(
                    card_header("Female proportions"),
                    card_body(div(class = "table-wrap", DTOutput("sex_female_props")))
                )
            )
        ),

        comparison_tab(
            title = "By age",
            id = "age",
            data_controls = list(
                cell_type_input("age_cell_type", cell_types),
                decade_panels_input("age_panel_decades"),
                ethnicity_picker_input("age_ethnicity"),
                imputed_input("age_imputed"),
                sex_input("age_sex"),
                tissue_groups_input("age_tissue_groups", tissue_groups)
            ),
            panels_ui = uiOutput("age_panels"),
            table_ui = panel_table_card("age_props")
        ),

        comparison_tab(
            title = "By ethnicity",
            id = "ethnicity",
            data_controls = list(
                cell_type_input("ethnicity_cell_type", cell_types),
                ethnicity_panels_input("ethnicity_panel_groups"),
                age_picker_input("ethnicity_age"),
                imputed_input("ethnicity_imputed"),
                sex_input("ethnicity_sex"),
                tissue_groups_input("ethnicity_tissue_groups", tissue_groups)
            ),
            panels_ui = uiOutput("ethnicity_panels"),
            table_ui = panel_table_card("ethnicity_props")
        ),

        !!!viz_tabs_ui(viz, cell_types, tissue_groups, enabled = extra_tabs),

        nav_panel(
            "Full data",
            card(
                card_header("Immune proportions across all samples"),
                card_body(
                    p(
                        class = "lead-note",
                        "Every modelled cell-type proportion in the atlas: one row ",
                        "per sample per cell type. Use the column filters to subset, ",
                        "and the buttons to copy or export the current selection."
                    ),
                    div(class = "table-wrap", DTOutput("full_data"))
                )
            )
        ),

        nav_panel(
            "About",
            card(
                class = "about-card",
                card_body(
                    h3("Immune Proportion Map"),
                    p(
                        class = "lead-note",
                        "Body maps of immune cell composition across human tissues, ",
                        "sex, age and ethnicity."
                    ),
                    div(
                        class = "about-columns",
                        h4("What the colours mean"),
                        p(
                            "Each organ is shaded by the ", strong("mean cell proportion"),
                            " of the selected immune cell type in that tissue group, ",
                            "averaged over every sample matching the current filters. ",
                            "Panels within a tab share one colour scale, so they are ",
                            "directly comparable. The tables report that mean alongside ",
                            "the number of samples behind it, their standard deviation, ",
                            "and a 95% confidence interval on the mean."
                        ),
                        p(
                            strong("Two different intervals appear in this app."),
                            " The panel tables give a confidence interval on the ",
                            "across-sample mean. The Full data tab gives each sample's own ",
                            "posterior credible interval from the model. ",
                            "They are labelled CI and CrI, respectively."
                        ),
                        h4("Where the numbers come from"),
                        p(
                            "Values are posterior predictions from a ",
                            a("sccomp", href = "https://github.com/MangiolaLaboratory/sccomp",
                                target = "_blank", rel = "noopener"),
                            " Bayesian compositional model fitted to the immune atlas, ",
                            "with tissue-group-varying effects of age, sex and ethnicity. ",
                            "This covers 4,483 samples and 26 immune cell types.",
                            "Age is recorded as decade of life."
                        ),
                        h4("How tissues are drawn"),
                        p(
                            "The atlas groups samples into 27 tissue groups, while ",
                            a("gganatogram", href = "https://github.com/jespermaag/gganatogram",
                                target = "_blank", rel = "noopener"),
                            " draws individual organs. One tissue group therefore shades ",
                            "several organs with the same value; 'respiratory system' ",
                            "fills the lungs, bronchi, diaphragm and pulmonary valve, for example ",
                            "Organs whose tissue group has no data for the current selection ",
                            "are left unfilled."
                        ),
                        p(
                            "Two groups cannot be drawn at all. ", strong("Thymus"),
                            " has no counterpart among gganatogram's organs, and ",
                            strong("digestive system (general)"),
                            " has none that is not already claimed by the stomach, small ",
                            "intestine or large intestine groups. Both still appear in ",
                            "the panel tables, with an empty Organ cell, and on the plot ",
                            "tabs."
                        ),
                        h4("Coverage is uneven, so read the sample counts"),
                        p(
                            "2,812 of the 4,483 samples are blood. At the other end, the ",
                            "gallbladder is represented by a single sample, the brainstem ",
                            "by three and the pancreas by four, and the 90-99 decade by ",
                            "eleven. Those tissues are still shaded, but a mean over ",
                            "one sample is not comparable to a mean over a thousand. ",
                            "Every panel header carries its sample count, and every table ",
                            "has a Samples column. Use them."
                        ),
                        h4("Sex and ethnicity"),
                        p(
                            "Ethnicity is recorded as one of nine groups, each either ",
                            "reported or imputed; about 26% of samples carry an imputed ",
                            "label. They are included by default and can be excluded with ",
                            "the sidebar toggle, which visibly shrinks the smaller groups. "
                        ),
                        p(
                            "Sex has three recorded values. gganatogram provides a male ",
                            "and a female body outline and nothing else, so the 136 ",
                            "samples with unknown sex are absent from the body maps. They ",
                            "are present in the Full data tab and on the plot tabs."
                        ),
                        h4("Built with"),
                        p(
                            "Body maps are drawn by ",
                            a("gganatogram",
                                href = "https://github.com/jespermaag/gganatogram",
                                target = "_blank", rel = "noopener"),
                            ". The interactive plot tabs are built from ",
                            a("VizModules",
                                href = "https://cran.r-project.org/package=VizModules",
                                target = "_blank", rel = "noopener"),
                            ", a collection of Shiny modules that decouple plotting from ",
                            "the underlying data and expose every aesthetic as a control."
                        ),
                        h4("Reference"),
                        # Built as one HTML string: htmltools puts each child on its own
                        # line, which the browser renders as a space before the comma.
                        p(HTML(paste0(
                            "Mangiola et al., <a href=\"", PREPRINT_URL,
                            "\" target=\"_blank\" rel=\"noopener\">A multi-organ map of ",
                            "the human immune system across age, sex and ethnicity</a>, bioRxiv."
                        )))
                    )
                )
            )
        ),

        nav_spacer(),
        nav_item(input_dark_mode(id = "mode", mode = "light")),
        nav_item(
            tags$a(
                class = "nav-link",
                href = PREPRINT_URL,
                target = "_blank",
                rel = "noopener",
                icon("file-lines"), " Preprint"
            )
        )
    )


    # ---------------------------------------------------------- server ----

    server <- function(input, output, session) {
        mode <- reactive(input$mode)

        output$full_data <- renderDT(full_data_widget)

        # -- By sex ------------------------------------------------------------
        #
        # Two panels differing only in sex, sharing one colour scale so male and
        # female are directly comparable.

        sex_opts <- aesthetics_server("sex")

        sex_data <- reactive({
            # An empty picker is a valid state and must fall through; only the cell
            # type is genuinely required before there is anything to draw.
            req(input$sex_cell_type)

            make_panel_data(
                prop_data,
                values = SEXES,
                by = "sex",
                celltype = input$sex_cell_type,
                tissue_groups = none_if_null(input$sex_tissue_groups),
                age = none_if_null(input$sex_age),
                ethnicity = none_if_null(input$sex_ethnicity),
                include_imputed = isTRUE(input$sex_imputed)
            )
        })

        sex_max <- reactive({
            safe_max(lapply(sex_data(), drawn_values))
        })

        for (s in SEXES) {
            local({
                key <- s
                anatogram_panel_server(
                    id = key,
                    panel_data = reactive(sex_data()[[key]]),
                    sex = reactive(key),
                    opts = sex_opts,
                    max_val = sex_max,
                    label = key,
                    mode = mode
                )
            })
        }

        output$sex_legend <- renderPlot(
            legend_bar(sex_max(), sex_opts(), fg = ink(mode())),
            bg = "transparent",
            alt = LEGEND_ALT
        ) |>
            bindCache(sex_max(), sex_opts(), mode())

        output$sex_male_props <- renderDT(props_table(sex_data()$male))
        output$sex_female_props <- renderDT(props_table(sex_data()$female))

        # -- By age ------------------------------------------------------------
        #
        # Ten decades is too many to draw at once, so the user picks which get a
        # panel. Every decade's module server is registered up front and only the
        # selected panels' UI is rendered. Shiny suspends outputs that are not
        # in the DOM, so an unselected decade costs nothing, and the ids stay
        # stable so the render cache survives a change of selection.

        age_opts <- aesthetics_server("age")

        age_panels <- reactive(none_if_null(input$age_panel_decades))

        output$age_panels <- renderUI(selected_panel_grid(age_panels()))

        age_data <- reactive({
            req(input$age_cell_type, input$age_sex)

            make_panel_data(
                prop_data,
                values = age_panels(),
                by = "age",
                celltype = input$age_cell_type,
                tissue_groups = none_if_null(input$age_tissue_groups),
                ethnicity = none_if_null(input$age_ethnicity),
                include_imputed = isTRUE(input$age_imputed),
                sex = input$age_sex
            )
        })

        age_max <- reactive({
            safe_max(lapply(age_data(), drawn_values))
        })

        age_selected <- reactiveVal(NULL)

        anatogram_panel_row_server(
            values = AGE_DECADES,
            panel_data = age_data,
            sex = reactive(input$age_sex),
            opts = age_opts,
            max_val = age_max,
            selected = age_selected,
            mode = mode
        )

        # A table left over from a previous selection would be misleading.
        observeEvent(age_data(), age_selected(NULL), ignoreInit = TRUE)

        output$age_legend <- renderPlot(
            legend_bar(age_max(), age_opts(), fg = ink(mode())),
            bg = "transparent",
            alt = LEGEND_ALT
        ) |>
            bindCache(age_max(), age_opts(), mode())

        output$age_props <- renderDT({
            sel <- age_selected()
            validate(need(
                sel,
                paste(
                    "No panel selected yet.",
                    "Choose “Data” beneath any body map above",
                    "to show that panel's numbers here."
                )
            ))
            props_table(sel$data, caption = paste0("Cell proportions: ", sel$label))
        })

        # -- By ethnicity ------------------------------------------------------

        ethnicity_opts <- aesthetics_server("ethnicity")

        ethnicity_panels <- reactive(none_if_null(input$ethnicity_panel_groups))

        output$ethnicity_panels <- renderUI(selected_panel_grid(ethnicity_panels()))

        ethnicity_data <- reactive({
            req(input$ethnicity_cell_type, input$ethnicity_sex)

            make_panel_data(
                prop_data,
                values = ethnicity_panels(),
                by = "ethnicity",
                celltype = input$ethnicity_cell_type,
                tissue_groups = none_if_null(input$ethnicity_tissue_groups),
                age = none_if_null(input$ethnicity_age),
                include_imputed = isTRUE(input$ethnicity_imputed),
                sex = input$ethnicity_sex
            )
        })

        ethnicity_max <- reactive({
            safe_max(lapply(ethnicity_data(), drawn_values))
        })

        ethnicity_selected <- reactiveVal(NULL)

        anatogram_panel_row_server(
            values = ETHNICITY_PANELS,
            panel_data = ethnicity_data,
            sex = reactive(input$ethnicity_sex),
            opts = ethnicity_opts,
            max_val = ethnicity_max,
            selected = ethnicity_selected,
            mode = mode
        )

        observeEvent(ethnicity_data(), ethnicity_selected(NULL), ignoreInit = TRUE)

        output$ethnicity_legend <- renderPlot(
            legend_bar(ethnicity_max(), ethnicity_opts(), fg = ink(mode())),
            bg = "transparent",
            alt = LEGEND_ALT
        ) |>
            bindCache(ethnicity_max(), ethnicity_opts(), mode())

        output$ethnicity_props <- renderDT({
            sel <- ethnicity_selected()
            validate(need(
                sel,
                paste(
                    "No panel selected yet.",
                    "Choose “Data” beneath any body map above",
                    "to show that panel's numbers here."
                )
            ))
            props_table(sel$data, caption = paste0("Cell proportions: ", sel$label))
        })

        # -- Plot tabs ---------------------------------------------------------

        viz_tabs_server(viz, input, output, enabled = extra_tabs)
    }


    shinyApp(ui = ui, server = server)
}
