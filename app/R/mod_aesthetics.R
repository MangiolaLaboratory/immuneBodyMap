# Plot aesthetics module.
#
# This block of controls was previously copy-pasted verbatim into all three
# comparison tabs. The server returns a single debounced reactive so a colour
# picker drag re-renders the tab's panels once, rather than on every frame.

aesthetics_ui <- function(id) {
    ns <- NS(id)

    tagList(
        splitLayout(
            cellWidths = c("40%", "60%"),
            selectInput(
                ns("palette"),
                with_help(
                    "Palette",
                    "All options are perceptually uniform and safe for the common
                 forms of colour vision deficiency."
                ),
                choices = PALETTES
            ),
            prettyCheckbox(
                ns("reverse"),
                "Reverse palette",
                value = FALSE,
                status = "primary",
                shape = "curve"
            )
        ),
        sliderInput(
            ns("opacity"),
            "Opacity",
            min = 0.1,
            max = 1,
            value = DEFAULT_OPACITY,
            step = 0.05,
            ticks = FALSE
        ),
        splitLayout(
            cellWidths = c("40%", "60%"),
            colourInput(ns("outline_colour"), "Outline colour", value = DEFAULT_OUTLINE_COLOUR),
            prettyCheckbox(
                ns("outline"),
                "Draw body outline",
                value = TRUE,
                status = "primary",
                shape = "curve"
            )
        ),
        selectInput(
            ns("dl_type"),
            with_help(
                "Download format",
                "Format used by the download button on each panel. PDF and SVG
                 are vector formats suitable for figures."
            ),
            choices = DL_FORMATS
        )
    )
}

#' @return A debounced reactive holding the current aesthetic options.
aesthetics_server <- function(id, debounce_ms = AESTHETICS_DEBOUNCE_MS) {
    moduleServer(id, function(input, output, session) {
        opts <- reactive({
            list(
                palette = input$palette,
                opacity = input$opacity,
                direction = if (isTRUE(input$reverse)) -1L else 1L,
                outline = isTRUE(input$outline),
                outline_colour = input$outline_colour,
                dl_type = input$dl_type
            )
        })

        debounce(opts, debounce_ms)
    })
}
