# Immune Proportion Map
#
# Interactive body maps of immune cell composition across human tissues, built
# on precomputed sccomp posterior predictions. See R/ for the helpers and
# modules; this file wires them into the four tabs.

library(shiny)
library(bslib)
library(gganatogram)
library(shinyWidgets)
library(DT)
library(viridis)
library(colourpicker)
library(svglite)
library(ggplot2)

# Sourced explicitly rather than via Shiny's R/ autoloading so the app behaves
# identically under runApp() and a deployed bundle. local = TRUE keeps the
# helpers in the app's own environment instead of the global one.
for (f in c(
    "constants.R", "organ_maps.R", "data.R", "plots.R", "tables.R",
    "inputs.R", "mod_aesthetics.R", "mod_anatogram_panel.R"
)) {
    source(file.path("R", f), local = TRUE)
}

prop_data <- load_prop_data()

CELL_TYPES <- levels(prop_data$cell_type_unified_ensemble)
TISSUE_GROUPS <- levels(prop_data$tissue_groups)

# Static widget: the same for every session, so build it once at startup rather
# than re-serialising 18,000 rows per connection.
FULL_DATA_WIDGET <- full_data_table(prop_data)


# ---------------------------------------------------------------- theme ----

app_theme <- bs_theme(
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


# ------------------------------------------------------------------ ui ----

#' A comparison tab: controls on the left, a wrapping grid of panels, a shared
#' colour bar, and the numbers underneath.
comparison_tab <- function(title, id, data_controls, panels, table_ui,
                           panel_width = 1 / length(panels),
                           panel_fixed = FALSE) {
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
            # A fractional width pins the column count so every panel sits in
            # one row; styles.css reverts to wrapping on narrow viewports. The
            # tight gap is what makes seven panels fit on a laptop.
            layout_column_wrap(
                width = panel_width,
                fixed_width = panel_fixed,
                heights_equal = "row",
                gap = "0.4rem",
                class = "panel-grid",
                !!!panels
            ),
            plotOutput(paste0(id, "_legend"), height = "100px"),
            table_ui
        )
    )
}

ui <- page_navbar(
    title = tags$span(class = "brand", "Immune Proportion Map"),
    window_title = "Immune Proportion Map",
    id = "nav",
    theme = app_theme,
    lang = "en",
    fillable = FALSE,
    header = tags$head(
        tags$link(rel = "stylesheet", href = "styles.css"),
        # bslib's own busy indicators: a page-level pulse plus a spinner on any
        # output that is recalculating. No extra package needed.
        useBusyIndicators()
    ),

    comparison_tab(
        title = "By sex",
        id = "sex",
        data_controls = list(
            cell_type_input("sex_cell_type", CELL_TYPES),
            ethnicity_input("sex_ethnicity"),
            age_slider_input("sex_age"),
            tissue_groups_input("sex_tissue_groups", TISSUE_GROUPS)
        ),
        panel_width = SEX_CARD_WIDTH,
        panel_fixed = TRUE,
        panels = list(
            anatogram_panel_ui("male", "Male", show_data = FALSE, height = SEX_CARD_HEIGHT),
            anatogram_panel_ui("female", "Female", show_data = FALSE, height = SEX_CARD_HEIGHT)
        ),
        # Only two panels, so both tables fit side by side and there is no need
        # to make the user ask for them.
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
            cell_type_input("age_cell_type", CELL_TYPES),
            ethnicity_input("age_ethnicity"),
            sex_input("age_sex"),
            tissue_groups_input("age_tissue_groups", TISSUE_GROUPS)
        ),
        panels = anatogram_panel_row_ui(AGE_BINS),
        table_ui = card(
            min_height = 560,
            card_header("Panel data"),
            # fillable = FALSE stops the card body from becoming its own scroll
            # container, so a full page of rows is read by scrolling the page.
            card_body(
                fillable = FALSE,
                div(class = "table-wrap", DTOutput("age_props"))
            )
        )
    ),

    comparison_tab(
        title = "By ethnicity",
        id = "ethnicity",
        data_controls = list(
            cell_type_input("ethnicity_cell_type", CELL_TYPES),
            age_select_input("ethnicity_age"),
            sex_input("ethnicity_sex"),
            tissue_groups_input("ethnicity_tissue_groups", TISSUE_GROUPS)
        ),
        panels = anatogram_panel_row_ui(ETHNICITIES),
        table_ui = card(
            min_height = 560,
            card_header("Panel data"),
            card_body(
                fillable = FALSE,
                div(class = "table-wrap", DTOutput("ethnicity_props"))
            )
        )
    ),

    nav_panel(
        "Full data",
        card(
            card_header("Immune proportions across all samples"),
            card_body(
                p(
                    class = "lead-note",
                    "Every modelled cell-type proportion in the atlas. Use the ",
                    "column filters to subset, and the buttons to copy or export ",
                    "the current selection."
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
                h4("What the colours mean"),
                p(
                    "Each organ is shaded by the ", strong("mean cell proportion"),
                    " of the selected immune cell type in that tissue group: the ",
                    "share of cells at that anatomical site that are of the chosen ",
                    "type. Panels within a tab share one colour scale, so they are ",
                    "directly comparable; the scale is rescaled when the selection ",
                    "changes. The tables report the same value alongside the lower ",
                    "and upper bounds of its 95% credible interval."
                ),
                h4("Where the numbers come from"),
                p(
                    "Values are posterior predictions from a ",
                    a("sccomp", href = "https://github.com/MangiolaLaboratory/sccomp",
                        target = "_blank", rel = "noopener"),
                    " Bayesian compositional model fitted to the immune atlas, ",
                    "with tissue-group-varying effects of age, sex and ethnicity. ",
                    "Each prediction summarises 5,000 simulated cells, so it is a ",
                    "modelled estimate rather than a direct cell count. Age bins ",
                    "use sex-specific boundaries."
                ),
                h4("How tissues are drawn"),
                p(
                    "The atlas groups samples into 24 tissue groups, while ",
                    a("gganatogram", href = "https://github.com/jespermaag/gganatogram",
                        target = "_blank", rel = "noopener"),
                    " draws individual organs. One tissue group therefore shades ",
                    "several organs with the same value -- 'respiratory system' ",
                    "fills the lungs, bronchi, diaphragm and pulmonary valve alike. ",
                    "Organs whose tissue group has no data for the current selection ",
                    "are left unfilled and are omitted from the tables."
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


# -------------------------------------------------------------- server ----

server <- function(input, output, session) {
    mode <- reactive(input$mode)

    output$full_data <- renderDT(FULL_DATA_WIDGET)

    # -- By sex ------------------------------------------------------------
    #
    # Two panels differing only in sex, sharing one colour scale so male and
    # female are directly comparable.

    sex_opts <- aesthetics_server("sex")

    sex_data <- reactive({
        # An empty tissue selection is a valid state and must fall through;
        # a missing cell type / age / ethnicity is just an uninitialised input.
        req(input$sex_cell_type, input$sex_ethnicity, input$sex_age)

        make_panel_data(
            prop_data,
            values = SEXES,
            by = "sex",
            celltype = input$sex_cell_type,
            tissue_groups = none_if_null(input$sex_tissue_groups),
            age = input$sex_age,
            ethnicity = input$sex_ethnicity
        )
    })

    sex_max <- reactive({
        safe_max(lapply(sex_data(), `[[`, "value"))
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

    age_opts <- aesthetics_server("age")

    age_data <- reactive({
        req(input$age_cell_type, input$age_ethnicity, input$age_sex)

        make_panel_data(
            prop_data,
            values = AGE_BINS,
            by = "age",
            celltype = input$age_cell_type,
            tissue_groups = none_if_null(input$age_tissue_groups),
            ethnicity = input$age_ethnicity,
            sex = input$age_sex
        )
    })

    age_max <- reactive({
        safe_max(lapply(age_data(), `[[`, "value"))
    })

    age_selected <- reactiveVal(NULL)

    anatogram_panel_row_server(
        values = AGE_BINS,
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
        props_table(sel$data, caption = paste0("Cell proportions — ", sel$label))
    })

    # -- By ethnicity ------------------------------------------------------

    ethnicity_opts <- aesthetics_server("ethnicity")

    ethnicity_data <- reactive({
        req(input$ethnicity_cell_type, input$ethnicity_age, input$ethnicity_sex)

        make_panel_data(
            prop_data,
            values = ETHNICITIES,
            by = "ethnicity",
            celltype = input$ethnicity_cell_type,
            tissue_groups = none_if_null(input$ethnicity_tissue_groups),
            age = input$ethnicity_age,
            sex = input$ethnicity_sex
        )
    })

    ethnicity_max <- reactive({
        safe_max(lapply(ethnicity_data(), `[[`, "value"))
    })

    ethnicity_selected <- reactiveVal(NULL)

    anatogram_panel_row_server(
        values = ETHNICITIES,
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
        props_table(sel$data, caption = paste0("Cell proportions — ", sel$label))
    })
}

shinyApp(ui = ui, server = server)
