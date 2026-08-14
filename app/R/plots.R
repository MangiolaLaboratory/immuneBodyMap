# Plot construction. Every anatogram in the app comes from anatogram_plot(),
# which used to be fifteen copy-pasted blocks.

#' Foreground colour for plot text, following the active colour mode.
#'
#' Plots are drawn on a transparent background so the card shows through, which
#' means fixed dark text would vanish in dark mode.
ink <- function(mode = "light") {
    if (identical(mode, "dark")) "#dee2e6" else "#495057"
}

#' Theme shared by every plot: no axes, no background, transparent throughout.
theme_anatogram <- function(fg = ink()) {
    theme_void(base_size = 12) +
        theme(
            text = element_text(colour = fg),
            plot.background = element_rect(fill = NA, colour = NA),
            panel.background = element_rect(fill = NA, colour = NA),
            legend.background = element_rect(fill = NA, colour = NA),
            plot.margin = margin(2, 2, 2, 2)
        )
}

#' Placeholder shown when a selection matches no data.
empty_plot <- function(message = "No results", fg = ink()) {
    ggplot() +
        geom_text(
            aes(x = 0.5, y = 0.5, label = message),
            inherit.aes = FALSE,
            check_overlap = TRUE,
            colour = fg,
            size = 4.2,
            alpha = 0.7
        ) +
        theme_anatogram(fg) +
        theme(plot.margin = margin(1, 1, 1, 1, "cm"))
}

#' A gganatogram body heatmap for one panel.
#'
#' @param pdata Frame from make_plot_data().
#' @param sex "male" or "female"; selects the body outline.
#' @param opts Aesthetic options from aesthetics_server().
#' @param max_val Upper limit of the fill scale, shared across a tab's panels so
#'   they stay comparable.
#' @param show_legend Panels normally hide their legend and defer to the single
#'   shared colour bar rendered once per tab.
anatogram_plot <- function(pdata, sex, opts, max_val, show_legend = FALSE, fg = ink()) {
    # No rows at all (everything deselected) or no values at all (a filter
    # combination with no samples) both mean there is nothing to draw.
    if (is.null(pdata) || !nrow(pdata) || all(is.na(pdata$value))) {
        return(empty_plot(fg = fg))
    }

    gganatogram(
        data = pdata,
        sex = sex,
        fill = "value",
        organism = "human",
        outline = isTRUE(opts$outline),
        fillOutline = opts$outline_colour
    ) +
        scale_fill_viridis(
            option = opts$palette,
            alpha = opts$opacity,
            direction = opts$direction,
            limits = c(0, max_val),
            guide = if (show_legend) "colourbar" else "none"
        ) +
        # Without this the body is stretched to fill whatever shape the card
        # happens to be -- wide cards produced visibly squat, broad figures.
        # x and y share one coordinate space, so ratio = 1 is true anatomy.
        coord_fixed(ratio = 1) +
        guides(fill = if (show_legend) colourbar_guide(fg) else "none") +
        theme_anatogram(fg) +
        theme(
            legend.position = if (show_legend) "bottom" else "none",
            legend.title.position = "top",
            legend.margin = margin(t = 8, b = 2)
        )
}

#' The colour bar carried by a downloaded figure.
#'
#' On screen the panels defer to one shared colour bar per tab, but a downloaded
#' plot travels on its own and has to carry its own scale. The sizing lives in
#' the guide's own theme -- plot-level legend.key.width does not size a
#' colourbar, which left the bar wider than the page and clipped at both ends.
colourbar_guide <- function(fg = ink(), width_in = 2.1) {
    guide_colourbar(
        title = "Mean cell proportion",
        theme = theme(
            legend.key.width = unit(width_in, "in"),
            legend.key.height = unit(0.16, "in"),
            legend.title = element_text(colour = fg, size = 10, hjust = 0.5),
            legend.text = element_text(colour = fg, size = 8.5),
            legend.ticks = element_line(colour = fg, linewidth = 0.3),
            legend.frame = element_blank()
        )
    )
}

#' The single colour bar shared by all panels on a tab.
#'
#' Each panel used to draw its own identical legend, which on a 210px-wide
#' anatogram consumed a large share of the panel.
legend_bar <- function(max_val, opts, title = "Mean cell proportion", fg = ink()) {
    ramp <- data.frame(x = seq(0, max_val, length.out = 256), y = 1)

    ggplot(ramp, aes(x = x, y = y, fill = x)) +
        geom_raster() +
        scale_fill_viridis(
            option = opts$palette,
            alpha = opts$opacity,
            direction = opts$direction,
            limits = c(0, max_val),
            guide = "none"
        ) +
        scale_x_continuous(
            expand = expansion(0),
            labels = function(v) format(signif(v, 2), scientific = FALSE, trim = TRUE)
        ) +
        scale_y_continuous(expand = expansion(0)) +
        labs(x = NULL, y = NULL, title = title) +
        theme_anatogram(fg) +
        theme(
            # Sized to sit alongside the surrounding UI text (~14px) rather than
            # the ggplot default, which read as fine print next to the cards.
            axis.text.x = element_text(colour = fg, size = 16, margin = margin(t = 5)),
            axis.ticks.x = element_line(colour = fg, linewidth = 0.4),
            axis.ticks.length.x = unit(5, "pt"),
            plot.title = element_text(colour = fg, size = 18, hjust = 0, margin = margin(b = 7)),
            # The end labels are centred on the ends of the bar, so half of each
            # sits outside the panel -- without side margins they get clipped.
            plot.margin = margin(2, 24, 4, 24)
        )
}

#' Write a plot to disk in the format the user actually asked for.
#'
#' ggsave() infers its device from the path it is handed. Shiny hands
#' downloadHandler a temp file with no meaningful extension, so without an
#' explicit device every download came out a PNG regardless of the selected
#' format -- only the filename changed.
save_plot <- function(file, plot, dl_type,
                      width = PLOT_EXPORT_WIDTH, height = PLOT_EXPORT_HEIGHT) {
    ggsave(
        filename = file,
        plot = plot,
        device = sub("^\\.", "", dl_type),
        width = width,
        height = height,
        units = "in",
        dpi = 300,
        # Screen plots are transparent so the card shows through; an exported
        # figure wants an opaque background.
        bg = "white"
    )
}
