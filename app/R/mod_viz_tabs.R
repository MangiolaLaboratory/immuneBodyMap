# The optional VizModules plot tabs.
#
# The body maps answer "where is this cell type?" well and everything else
# badly: they show one cell type at a time, collapse every sample behind a
# tissue group to a single mean, and cannot show a trend. These four tabs cover
# what that leaves out: the whole composition at once, the spread behind each
# mean, change across the decades, and the cell-type-by-tissue structure.
#
# Every module here is a VizModules trio sharing one id: <mod>InputsUI() takes a
# plain data frame (it needs column names at UI-build time), <mod>OutputUI()
# places the plot, and <mod>Server() takes a reactive(). The same `defaults`
# list must reach both the UI and the server, which is why they are built by
# functions rather than written out twice.
#
# The keys below are input IDs and are NOT uniform across modules. They were
# read out of each module's own InputsUI rather than guessed. get_default()
# falls back silently on an unrecognised key or an out-of-range value, so a typo
# here looks like "the default did not work" rather than failing.

VIZ_SIDEBAR_WIDTH <- 430

# What the Composition tab can split panels by. "" is the module's own value for
# "no facet", so it doubles as the facet.by default.
DOT_SPLITS <- c(
    "None" = "",
    "Sex" = "sex",
    "Decade of life" = "age_decade",
    "Ethnicity" = "ethnicity",
    "Assay" = "assay"
)

#' Is the heatmap module usable in this installation?
#'
#' ComplexHeatmap and InteractiveComplexHeatmap are Suggests of VizModules, not
#' hard dependencies, and they come from Bioconductor. The module errors with a
#' BiocManager::install() message when they are missing, so the tab is left out
#' entirely rather than added and broken.
heatmap_available <- function() {
    requireNamespace("ComplexHeatmap", quietly = TRUE) &&
        requireNamespace("InteractiveComplexHeatmap", quietly = TRUE)
}

#' Everything the plot tabs need, built once at startup.
#'
#' These summaries are the same for every session, so they are computed with the
#' rest of the startup work rather than per connection. The dot plot's own frame
#' is reactive (its granularity follows the split control), so only the ungrouped
#' summary the heatmap needs is precomputed here.
viz_data <- function(prop_data) {
    summary_df <- tissue_celltype_summary(prop_data)

    list(
        samples = prop_data,
        dot = summary_df,
        heat = celltype_tissue_matrix(summary_df)
    )
}

#' Give a module's plot output a different default size.
#'
#' The OutputUI helpers take no width or height; they return a plotlyOutput()
#' whose inline style is width:100%;height:400px. Rewriting that inline style is
#' the only way to change the default that does not fight jqui_resizable: a
#' stylesheet rule would need !important to beat the inline value, and would
#' then also beat the inline value the resize handle writes, so dragging the
#' corner would appear to do nothing.
size_plot_output <- function(tag, width = "100%", height = "400px") {
    style <- paste0("width:", width, ";height:", height, ";")

    # The helper nests the output in tagLists inside tagLists, and how deeply
    # depends on `resizable`, so walk the tree rather than indexing into it.
    # "html-widget-output" is what identifies the div itself: jqui's companion
    # iframe and the wrapping lists carry no such class.
    walk <- function(x) {
        if (inherits(x, "shiny.tag")) {
            cls <- x$attribs$class
            if (!is.null(cls) && grepl("html-widget-output", cls, fixed = TRUE)) {
                x$attribs$style <- style
                return(x)
            }
            x$children <- lapply(x$children, walk)
            return(x)
        }

        if (is.list(x)) {
            out <- lapply(x, walk)
            # lapply drops the tagList class htmltools dispatches on.
            attributes(out) <- attributes(x)
            return(out)
        }

        x
    }

    walk(tag)
}

# Size is the proportion and colour is the evidence behind it, so a big pale dot
# reads as "large but thinly sampled" at a glance.
#
# The data is pre-aggregated AND size.by is set. Both matter: with size_by unset
# plotthis aggregates but takes first(fill_by); with it set it does not aggregate
# at all, so the frame must already be one row per x/y/facet cell. That is what
# ties `facet` here to the same column the data was summarised by.
dot_defaults <- function(facet = "") {
    list(
        x.data = "tissue_groups",
        y.data = "cell_type",
        size.by = "mean_proportion",
        fill.by = "n_samples",
        facet.by = facet,
        palette.name = "Spectral",
        palreverse = TRUE,
        size.max = 5,
        border.size = 0.25,
        axis.tickangle.x = 45
    )
}

# yPlot rather than a box plot module, because it draws the box and the
# per-sample points in one layer stack. That is the point of this tab: the body
# maps already show the mean, and what they hide is how much it varies between
# samples. `plots` is ordered, so jitter draws over the box.
ydist_defaults <- function() {
    list(
        var = "proportion_mean",
        group.by = "tissue_groups",
        color.by = "sex",
        plots = c("boxplot", "jitter"),
        jitter.size = 0.5,
        show.grid.x = FALSE,
        show.grid.y = FALSE,
        axis.tickangle.x = 45
    )
}

# linePlot() does its own mean and sd per x value when error.bar is on, so this
# tab is fed raw per-sample rows rather than anything pre-summarised.
line_defaults <- function() {
    list(
        x.value = "age_decade",
        y.value = "proportion_mean",
        group.by = "cell_type",
        error.bar = TRUE,
        plot.mode = "lines+markers"
    )
}

# "Rows", not "row", since the choices are None / Rows / Columns. Z-scoring by
# cell type is what makes the matrix readable: proportions run from under 0.001
# to over 0.5, so an unscaled heatmap is one bright row and 25 dark ones. The
# legend title has to say so, or the numbers look like proportions.
heat_defaults <- function(heat_df) {
    list(
        matrix.cols = setdiff(names(heat_df), "cell_type"),
        rowname.col = "cell_type",
        name = "zscore(mean proportion)",
        scale = "Rows",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        na_col = "grey90"
    )
}

#' Split control for the Composition tab.
#'
#' Drives the aggregation and the facet together. They cannot be two independent
#' controls: the plot is only correct while the data is summarised at exactly the
#' granularity the facet asks for.
dot_split_input <- function(id) {
    selectInput(
        id,
        with_help(
            "Split panels by",
            "Summarises the atlas within each level of the chosen variable and
             gives each one its own panel. Leave on None to pool every sample."
        ),
        choices = DOT_SPLITS,
        selected = ""
    )
}

#' One plot tab: our own filters above the module's own controls.
viz_tab <- function(title, inputs_ui, output_ui, ...) {
    nav_panel(
        title,
        layout_sidebar(
            sidebar = sidebar(
                width = VIZ_SIDEBAR_WIDTH,
                class = "viz-sidebar",
                ...,
                inputs_ui
            ),
            output_ui
        )
    )
}

#' The plot tabs, or nothing at all.
#'
#' Returns a list so app.r can splice it into page_navbar() with !!! whether or
#' not the tabs are switched on.
viz_tabs_ui <- function(vd, cell_types, tissue_groups, enabled = TRUE) {
    if (!enabled) {
        return(list())
    }

    tabs <- list(
        viz_tab(
            "Composition",
            # Built server-side: facet.by is a reactive default, and reactive
            # defaults are only reachable from a UI inside renderUI(). The
            # controls also have to be rebuilt when the split changes, because
            # the frame gains and loses a column.
            uiOutput("viz_dot_controls"),
            size_plot_output(
                plotthis_DotPlotOutputUI("viz_dot"),
                width = "min(100%, 900px)", height = "800px"
            ),
            dot_split_input("viz_dot_split"),
            hr()
        ),

        viz_tab(
            "Distributions",
            dittoViz_yPlotInputsUI("viz_dist", vd$samples, defaults = ydist_defaults()),
            dittoViz_yPlotOutputUI("viz_dist"),
            # Without this the plot pools all 26 cell types into each tissue's
            # box, which is a distribution of nothing in particular.
            cell_type_input("viz_dist_cell_type", cell_types),
            hr()
        ),

        viz_tab(
            "Ageing trend",
            linePlotInputsUI("viz_line", vd$samples, defaults = line_defaults()),
            linePlotOutputUI("viz_line"),
            selectInput(
                "viz_line_tissue",
                with_help(
                    "Tissue group",
                    "One tissue at a time; pooling them would average over very
                     different compositions. Blood carries by far the most samples."
                ),
                choices = tissue_groups,
                selected = if ("blood" %in% tissue_groups) "blood" else tissue_groups[1]
            ),
            hr()
        )
    )

    if (heatmap_available()) {
        tabs <- c(tabs, list(viz_tab(
            "Heatmap",
            ComplexHeatmap_HeatmapInputsUI(
                "viz_heat", vd$heat,
                defaults = heat_defaults(vd$heat)
            ),
            ComplexHeatmap_HeatmapOutputUI("viz_heat")
        )))
    }

    tabs
}

#' @param input The main server's input object. The tabs that carry their own
#'   filter read it directly, as the comparison tabs in app.r do.
viz_tabs_server <- function(vd, input, output, enabled = TRUE) {
    if (!enabled) {
        return(invisible(NULL))
    }

    # -- Composition -------------------------------------------------------

    dot_split <- reactive({
        s <- input$viz_dot_split
        if (is.null(s) || !nzchar(s)) character(0) else s
    })

    dot_data <- reactive(tissue_celltype_summary(vd$samples, by = dot_split()))

    # A reactive entry is the supported parent-to-child channel; driving the
    # control with updateSelectInput() from here would render the plot twice.
    dot_defs <- dot_defaults(facet = reactive({
        s <- dot_split()
        if (!length(s)) "" else s
    }))

    output$viz_dot_controls <- renderUI(
        plotthis_DotPlotInputsUI("viz_dot", dot_data(), defaults = dot_defs)
    )

    plotthis_DotPlotServer("viz_dot", data = dot_data, defaults = dot_defs)

    # -- Distributions -----------------------------------------------------

    dittoViz_yPlotServer(
        "viz_dist",
        data = reactive({
            req(input$viz_dist_cell_type)
            vd$samples[vd$samples$cell_type == input$viz_dist_cell_type, , drop = FALSE]
        }),
        defaults = ydist_defaults()
    )

    # -- Ageing trend ------------------------------------------------------

    linePlotServer(
        "viz_line",
        data = reactive({
            req(input$viz_line_tissue)
            vd$samples[vd$samples$tissue_groups == input$viz_line_tissue, , drop = FALSE]
        }),
        defaults = line_defaults()
    )

    # -- Heatmap -----------------------------------------------------------

    if (heatmap_available()) {
        ComplexHeatmap_HeatmapServer(
            "viz_heat",
            data = reactive(vd$heat),
            defaults = heat_defaults(vd$heat)
        )
    }

    invisible(NULL)
}
