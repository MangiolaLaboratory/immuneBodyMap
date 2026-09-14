# Immune Proportion Map

library(shiny)
library(bslib)
library(gganatogram)
library(shinyWidgets)
library(DT)
library(viridis)
library(colourpicker)
library(svglite)
library(ggplot2)

# Sourced explicitly so the app behaves
# identically under runApp() and a deployed bundle. local = TRUE keeps the
# helpers in the app's own environment instead of the global one.
for (f in c(
    "constants.R", "organ_maps.R", "data.R", "plots.R", "tables.R",
    "inputs.R", "mod_aesthetics.R", "mod_anatogram_panel.R", "mod_viz_tabs.R",
    "launch.R"
)) {
    source(file.path("R", f), local = TRUE)
}

# Set extra_tabs = FALSE here for the body maps alone
immune_map_app(extra_tabs = TRUE)
