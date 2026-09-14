# Categorical vocabulary of the app.
#
# These are the single source of truth: the UI builds its controls from them and
# load_prop_data() re-levels the data against them, so a mismatch surfaces as a
# startup error rather than as silently empty panels.

DATA_FILE <- "cell_proportion_for_shiny_app_20260910.rds"

# Age is decade of life. The file stores the decade number as a factor whose
# levels happen to be in the order 5, 1, 10, 2, 3, and so on. load_prop_data()
# re-levels against AGE_DECADE_LEVELS and relabels to AGE_DECADES in one step,
# so the rest of the app only ever sees the readable form. That matters in more
# places than the panel headings: DT filter dropdowns and plotly axes show these
# strings verbatim.
AGE_DECADE_LEVELS <- as.character(1:10)

AGE_DECADES <- c(
    "0-9", "10-19", "20-29", "30-39", "40-49",
    "50-59", "60-69", "70-79", "80-89", "90-99"
)

# Which decades get a panel on the "By age" tab by default. The tails are thin
# enough to mislead: across both sexes the 90-99 decade has 11 samples, and for
# females alone it covers a single tissue group. The app therefore opens on the
# six decades with real tissue coverage and leaves the rest one click away.
DEFAULT_DECADES <- AGE_DECADES[3:8]

# The nine base groups behind the 16 ethnicity_groups_imputed levels. Each is
# stored either bare or with an "_imputed" suffix; load_prop_data() splits that
# one column into an ethnicity and a logical imputed flag.
ETHNICITIES <- c(
    "European",
    "African",
    "East Asian",
    "Hispanic/Latin American",
    "Japanese",
    "Middle Eastern & North African",
    "Native American & Pacific Islander",
    "South Asian",
    "Other/Unknown"
)

# "Middle Eastern & North African" is a declared level with no samples at all in
# this release (not even an imputed twin), so it would only ever draw an empty
# panel. It stays in ETHNICITIES so the validation in load_prop_data() keeps
# accepting it if a later release fills it in.
ETHNICITY_PANELS <- setdiff(ETHNICITIES, "Middle Eastern & North African")

DEFAULT_ETHNICITY_PANELS <- c("European", "African", "East Asian", "South Asian")

# The two body outlines gganatogram can draw, and the levels actually present in
# the data. "unknown" (136 samples) has no outline, so it is reachable on the
# Full data and plot tabs but never gets a body map.
SEXES <- c("male", "female")
SEXES_ALL <- c(SEXES, "unknown")

# Panel headings, shortened where the full name will not fit a narrow card.
PANEL_LABELS <- c(
    "Hispanic/Latin American" = "Hisp./Latin Amer.",
    "Native American & Pacific Islander" = "Nat. Amer./Pac. Isl.",
    "Middle Eastern & North African" = "M. East./N. Afr."
)

# Every option is a perceptually uniform, colour-vision-deficiency-safe ramp.
# turbo is deliberately excluded: it is not perceptually uniform and reads
# poorly with the common forms of colour vision deficiency.
PALETTES <- c(
    "viridis", "magma", "plasma", "inferno",
    "cividis", "mako", "rocket"
)

DL_FORMATS <- c(".png", ".pdf", ".svg")

# A mid grey that reads against both the light and the dark theme.
DEFAULT_OUTLINE_COLOUR <- "#9aa3ad"

# The original default of 0.51 washed the viridis ramp out badly; 0.9 keeps the
# organ fills distinguishable while still letting the outline show through.
DEFAULT_OPACITY <- 0.9

# gganatogram's body outlines span roughly 107 x-units by 196 y-units, so a
# correctly proportioned figure is about 0.54 as wide as it is tall. gganatogram
# itself never calls coord_fixed(), so the figure stretches to whatever shape the
# plotting device happens to be. anatogram_plot() pins it instead.
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
AESTHETICS_DEBOUNCE_MS <- 600

# Raw column names are not presentable. Anything absent here falls through to
# the raw name, so adding a column to the data cannot break the tables.
#
# Note the deliberate split between the two kinds of uncertainty the app now
# carries. proportion_* are per-sample posterior quantities straight from the
# model ("CrI"); value / sd / se / ci_* summarise across the samples in a
# selection ("CI"). The labels are the only thing telling them apart.
PRETTY_COLNAMES <- c(
    organ = "Organ",
    tissue_groups = "Tissue group",
    value = "Mean proportion",
    sd = "SD across samples",
    se = "SE of mean",
    ci_lower = "95% CI lower",
    ci_upper = "95% CI upper",
    n_samples = "Samples",
    sample_id = "Sample ID",
    age_decade = "Decade of life",
    sex = "Sex",
    ethnicity = "Ethnicity",
    ethnicity_imputed = "Ethnicity imputed",
    cell_type = "Cell type",
    assay = "Assay",
    dataset_id = "Dataset",
    mean_proportion = "Mean proportion",
    proportion_mean = "Posterior mean proportion",
    proportion_lower = "95% CrI lower",
    proportion_upper = "95% CrI upper"
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
    x <- gsub("^_|_$", "", x)

    # Decade labels are digits, so this would otherwise return ids like "0_9",
    # and Shiny and bslib build jQuery selectors from these ids. "#0_9-plot" is
    # not a valid CSS identifier and never matches.
    sub("^(?=[0-9])", "d", x, perl = TRUE)
}

#' Display label for a panel, shortened where the full name will not fit.
panel_label <- function(x) {
    unname(ifelse(x %in% names(PANEL_LABELS), PANEL_LABELS[x], x))
}
