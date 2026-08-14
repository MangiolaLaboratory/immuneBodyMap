# Categorical vocabulary of the app.
#
# These are the single source of truth: the UI builds its controls from them and
# load_prop_data() re-levels the data against them, so a mismatch surfaces as a
# startup error rather than as silently empty panels.

# Chronological. The levels stored in the .rds are alphabetical, which is why
# load_prop_data() re-levels rather than reusing them.
AGE_BINS <- c(
    "Infancy",
    "Childhood",
    "Adolescence",
    "Young Adulthood",
    "Middle Age",
    "Senior"
)

# Display order for the ethnicity panels, matching the original app's layout.
ETHNICITIES <- c(
    "European",
    "African",
    "East Asian",
    "Hispanic/Latin American",
    "Native American & Pacific Islander",
    "South Asian",
    "Other/Unknown"
)

# Two of the ethnicity names are far too long for a narrow panel header.
ETHNICITY_LABELS <- c(
    "Hispanic/Latin American" = "Hisp./Latin Amer.",
    "Native American & Pacific Islander" = "Nat. Amer./Pac. Isl."
)

SEXES <- c("male", "female")

# Every option is a perceptually uniform, colour-vision-deficiency-safe ramp.
# turbo is deliberately excluded: it is not perceptually uniform and reads
# poorly with the common forms of colour vision deficiency.
PALETTES <- c(
    "viridis", "magma", "plasma", "inferno",
    "cividis", "mako", "rocket"
)

DL_FORMATS <- c(".png", ".pdf", ".svg")

DEFAULT_AGE_BIN <- "Adolescence"
DEFAULT_ETHNICITY <- "European"

# A mid grey that reads against both the light and the dark theme.
DEFAULT_OUTLINE_COLOUR <- "#9aa3ad"

# The original default of 0.51 washed the viridis ramp out badly; 0.9 keeps the
# organ fills distinguishable while still letting the outline show through.
DEFAULT_OPACITY <- 0.9

# gganatogram's body outlines span roughly 107 x-units by 196 y-units, so a
# correctly proportioned figure is about 0.54 as wide as it is tall. gganatogram
# itself never calls coord_fixed(), so the figure stretches to whatever shape the
# plotting device happens to be -- anatogram_plot() pins it instead.
ANATOGRAM_ASPECT <- 0.544

# The age and ethnicity tabs lay their panels out in a single row, so the column
# width follows from the panel count rather than being fixed. This height suits
# the resulting column width on a typical laptop.
# Height budget per card: roughly 34px of header plus 41px of footer.
PANEL_CARD_HEIGHT <- 370

# Below this viewport width a single row would squeeze the figures too far, so
# the grid falls back to wrapping (see .panel-grid in www/styles.css).
PANEL_GRID_WRAP_BELOW <- "1200px"

# The sex tab shows only two panels, so it can afford larger figures and keeps a
# fixed column width rather than stretching them across the row.
SEX_CARD_WIDTH <- "280px"
SEX_CARD_HEIGHT <- 520

# Exported figures. The width and the body height are matched to the anatogram's
# aspect so the figure is not letterboxed; the extra inch is the colour bar the
# download carries so the numbers are readable away from the app.
PLOT_EXPORT_WIDTH <- 3.81
PLOT_EXPORT_LEGEND_IN <- 0.9
PLOT_EXPORT_HEIGHT <- round(PLOT_EXPORT_WIDTH / ANATOGRAM_ASPECT + PLOT_EXPORT_LEGEND_IN, 2)

# Aesthetic controls fire continuously (colour picker drags, slider scrubs) and
# each event re-renders every panel on the tab, so their effect is debounced.
AESTHETICS_DEBOUNCE_MS <- 400

# Raw column names are not presentable. Anything absent here falls through to
# the raw name, so adding a column to the data cannot break the tables.
PRETTY_COLNAMES <- c(
    organ = "Organ",
    tissue_groups = "Tissue group",
    value = "Mean proportion",
    CI_lower = "95% CI lower",
    CI_upper = "95% CI upper",
    sample_id = "Sample ID",
    age_bin_sex_specific = "Age bin",
    sex = "Sex",
    ethnicity_groups = "Ethnicity",
    age_days = "Age (days)",
    age_days_scaled = "Age (scaled)",
    n = "Cells simulated",
    cell_type_unified_ensemble = "Cell type",
    proportion_mean = "Mean proportion",
    proportion_lower = "95% CI lower",
    proportion_upper = "95% CI upper"
)

PREPRINT_URL <- "https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3"

LEGEND_ALT <- paste(
    "Colour bar showing the scale shared by every panel on this tab,",
    "running from zero to the largest mean cell proportion in the current",
    "selection."
)

#' Turn a label into a module id safe for use in HTML element ids.
slugify <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    gsub("^_|_$", "", x)
}

#' Display label for a panel, shortened where the full name will not fit.
panel_label <- function(x) {
    unname(ifelse(x %in% names(ETHNICITY_LABELS), ETHNICITY_LABELS[x], x))
}
