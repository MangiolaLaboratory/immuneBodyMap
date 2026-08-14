# One anatogram panel: a card holding the body map, a download button and
# optionally a button that pushes the panel's numbers into the tab's shared
# table. Replaces fifteen copy-pasted plot / render / download / table blocks.

#' Alt text for a panel's body map.
panel_alt_text <- function(label) {
    paste0(
        "Anatomical diagram of a human body for ", label,
        ", with each organ shaded according to the mean proportion of the ",
        "selected immune cell type. Exact values are in the panel's data table."
    )
}

#' @param label Panel heading, e.g. "Infancy".
#' @param show_data Whether to offer the "Data" button. Tabs with few enough
#'   panels to show their tables permanently pass FALSE.
anatogram_panel_ui <- function(id, label, show_data = TRUE,
                               height = PANEL_CARD_HEIGHT) {
    ns <- NS(id)

    # Alongside a "Data" button the shorter label reads fine; on its own the
    # button spans the footer and needs to say what it does.
    actions <- list(
        downloadButton(
            ns("download"),
            if (show_data) "Plot" else "Download plot",
            class = "btn-sm btn-outline-secondary"
        )
    )
    if (show_data) {
        actions <- c(actions, list(
            actionButton(
                ns("show"),
                "Data",
                icon = icon("table"),
                class = "btn-sm btn-outline-secondary"
            )
        ))
    }

    card(
        height = height,
        full_screen = TRUE,
        card_header(label, class = "panel-title"),
        card_body(
            padding = 0,
            plotOutput(ns("plot"), height = "100%")
        ),
        card_footer(class = "panel-actions", actions)
    )
}

#' @param panel_data Reactive returning this panel's frame from make_plot_data().
#' @param sex Reactive returning "male" or "female".
#' @param opts Reactive of aesthetic options from aesthetics_server().
#' @param max_val Reactive upper bound of the fill scale, shared across the tab.
#' @param label Used for the download filename and the shared table's caption.
#' @param selected reactiveVal owned by the parent tab. Clicking "Data" writes
#'   this panel's label and frame into it, which is what drives the shared table
#'   -- previously thirteen one-line observers in the main server function.
#' @param mode Reactive colour mode ("light" / "dark") for plot text.
anatogram_panel_server <- function(id, panel_data, sex, opts, max_val, label,
                                   selected = NULL, mode = reactive("light")) {
    moduleServer(id, function(input, output, session) {
        plot_obj <- reactive({
            anatogram_plot(
                pdata = panel_data(),
                sex = sex(),
                opts = opts(),
                max_val = max_val(),
                fg = ink(mode())
            )
        })

        # Caching the rendered image, not the ggplot object: re-selecting a
        # previously viewed combination then costs nothing. The key must name
        # every input the plot depends on.
        output$plot <- renderPlot(
            plot_obj(),
            bg = "transparent",
            alt = panel_alt_text(label)
        ) |>
            bindCache(panel_data(), sex(), opts(), max_val(), mode())

        # The on-screen panels share one colour bar rendered beside the grid; a
        # downloaded figure has to carry its own, so it is rebuilt with the
        # legend attached rather than reusing the screen plot.
        download_obj <- reactive({
            anatogram_plot(
                pdata = panel_data(),
                sex = sex(),
                opts = opts(),
                max_val = max_val(),
                show_legend = TRUE,
                fg = "#495057"
            )
        })

        output$download <- downloadHandler(
            filename = function() {
                paste0(slugify(label), "_anatogram", opts()$dl_type)
            },
            content = function(file) {
                save_plot(file, download_obj(), opts()$dl_type)
            }
        )

        if (!is.null(selected)) {
            observeEvent(input$show, {
                selected(list(label = label, data = panel_data()))
            })
        }
    })
}

#' Wire up a whole row of panels from a set of levels.
#'
#' @param values The levels being compared, e.g. AGE_BINS.
#' @param panel_data Reactive returning a named list of frames, keyed by
#'   slugify(values).
anatogram_panel_row_ui <- function(values, ...) {
    lapply(values, function(v) {
        anatogram_panel_ui(id = slugify(v), label = panel_label(v), ...)
    })
}

anatogram_panel_row_server <- function(values, panel_data, ...) {
    # lapply rather than a for loop: each call gets its own frame, so the
    # reactive below closes over this panel's key instead of the last one.
    lapply(values, function(v) {
        key <- slugify(v)
        lab <- panel_label(v)

        anatogram_panel_server(
            id = key,
            panel_data = reactive(panel_data()[[key]]),
            label = lab,
            ...
        )
    })

    invisible(NULL)
}
