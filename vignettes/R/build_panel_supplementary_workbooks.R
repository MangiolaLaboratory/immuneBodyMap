#!/usr/bin/env Rscript
# Build separate plot-data, summary and test worksheets for each figure panel.
#
# Plot caches remain the source of truth for the values actually drawn. Each
# worksheet name retains the manuscript panel prefix (for example, `Fig4_B_`),
# while incompatible logical tables are kept in separate worksheets.

root <- Sys.getenv(
  "AGE_CLOCK_ROOT",
  unset = Sys.getenv("AGE_CLOCK_ARCHIVE", unset = "")
)
if (!nzchar(root)) {
  stop(
    "Set AGE_CLOCK_ROOT or AGE_CLOCK_ARCHIVE to rebuild SI workbooks from Age_Clock plot caches.\n",
    "Publication copies of SI tables already live in vignettes/data/source_tables/.",
    call. = FALSE
  )
}
helper <- file.path(root, "report", "R", "export_supplementary_tables.R")
if (!file.exists(helper)) {
  helper <- file.path(
    Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = normalizePath(file.path("..", ".."), mustWork = FALSE)),
    "vignettes", "R", "export_supplementary_tables.R"
  )
}
source(helper)
root <- age_clock_root(root)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

fig_dir <- file.path(root, "report", "Fig_files")
tissue_dir <- Sys.getenv("TISSUE_UMAP_DIR", unset = file.path(dirname(root), "Tissue_umap"))
data_dir <- file.path(root, "data")
age_data_dir <- file.path(data_dir, "age_fig")
suppl_dir <- file.path(root, "Supplementary_files")
selected_figure <- Sys.getenv("SUPPLEMENTARY_FIGURE", "all")
build_figure <- function(number) selected_figure %in% c("all", number)

read_plot <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing required plot cache: ", path, call. = FALSE)
    return(NULL)
  }
  readRDS(path)
}

read_plot_data <- function(path, required = TRUE) {
  x <- read_plot(path, required = required)
  if (is.null(x)) return(NULL)
  extract_plot_data(x, basename(path))
}

first_existing <- function(...) {
  paths <- list(...)
  for (path in paths) if (!is.null(path) && file.exists(path)) return(path)
  NULL
}

drop_empty_columns <- function(data) {
  data <- as.data.frame(data)
  has_value <- vapply(data, function(x) {
    if (is.character(x)) return(any(!is.na(x) & nzchar(x)))
    any(!is.na(x))
  }, logical(1))
  data <- data[has_value]
  rownames(data) <- NULL
  data
}

panel_sheet_prefix <- function(panel_name) {
  matched <- regmatches(
    panel_name,
    regexpr("^Fig[0-9]+_[A-Z]", panel_name, perl = TRUE)
  )
  if (length(matched) == 0L || !nzchar(matched)) panel_name else matched
}

split_panel_worksheets <- function(panels) {
  worksheets <- list()

  for (panel_name in names(panels)) {
    panel <- panels[[panel_name]]

    if (!inherits(panel, "panel_records")) {
      sheet_name <- sanitize_sheet_name(panel_name, names(worksheets))
      worksheets[[sheet_name]] <- drop_empty_columns(sanitize_for_excel(panel))
      next
    }

    for (record_type in names(panel)) {
      section <- drop_empty_columns(sanitize_for_excel(panel[[record_type]]))
      sheet_name <- sanitize_sheet_name(
        paste(panel_sheet_prefix(panel_name), record_type, sep = "_"),
        names(worksheets)
      )
      worksheets[[sheet_name]] <- section
    }
  }

  worksheets
}

combine_panel_records <- function(..., .label_col = "panel_component") {
  record_sets <- list(...)
  labels <- names(record_sets)
  if (is.null(labels) || any(!nzchar(labels))) {
    stop("Combined panel-record sets must have descriptive names.", call. = FALSE)
  }

  record_types <- unique(unlist(lapply(record_sets, names), use.names = FALSE))
  combined <- lapply(record_types, function(record_type) {
    pieces <- Map(function(records, label) {
      if (!record_type %in% names(records)) return(NULL)
      data <- as.data.frame(records[[record_type]])
      data[[.label_col]] <- label
      data
    }, record_sets, labels)
    pieces <- pieces[!vapply(pieces, is.null, logical(1))]
    as.data.frame(dplyr::bind_rows(pieces))
  })
  names(combined) <- record_types
  structure(combined, class = c("panel_records", "list"))
}

write_figure_panels <- function(panels, filename) {
  worksheets <- split_panel_worksheets(panels)
  export_supplementary_workbook(
    extra_sheets = worksheets,
    output_path = file.path(suppl_dir, filename)
  )
  invisible(names(worksheets))
}

add_asin_sqrt <- function(data, value_col = "proportion_mean") {
  value <- suppressWarnings(as.numeric(data[[value_col]]))
  data$proportion_asin_sqrt <- asin(sqrt(pmin(1, pmax(0, value))))
  data
}

summarise_trajectory <- function(data, group_cols, age_col, value_col) {
  data |>
    mutate(.age = suppressWarnings(as.numeric(as.character(.data[[age_col]])))) |>
    filter(is.finite(.age), is.finite(.data[[value_col]])) |>
    arrange(across(all_of(group_cols)), .age) |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      n_age_points = n(),
      min_age = min(.age),
      max_age = max(.age),
      start_value = first(.data[[value_col]]),
      end_value = last(.data[[value_col]]),
      absolute_change = end_value - start_value,
      relative_change = if_else(start_value == 0, NA_real_, end_value / start_value - 1),
      minimum_value = min(.data[[value_col]]),
      maximum_value = max(.data[[value_col]]),
      age_at_maximum = .age[which.max(.data[[value_col]])][1],
      .groups = "drop"
    ) |>
    as.data.frame()
}

boxplot_panel <- function(plot_data, selection_tests, cohort_col) {
  raw <- add_asin_sqrt(as.data.frame(plot_data)) |>
    mutate(plot_role = if_else(
      version == "prediction", "boxplot distribution", "adjusted-data points"
    ))

  raw_summary <- summarise_distribution(
    raw,
    group_cols = c("L3", "L3_abbrev", cohort_col, "version", "plot_role", "c_FDR"),
    value_col = "proportion_mean",
    value_label = "proportion"
  )
  transformed_summary <- summarise_distribution(
    raw,
    group_cols = c("L3", "L3_abbrev", cohort_col, "version", "plot_role", "c_FDR"),
    value_col = "proportion_asin_sqrt",
    value_label = "asin_sqrt_proportion"
  )

  bind_panel_records(
    plotted_observations = raw,
    distribution_summary_raw = raw_summary,
    distribution_summary_asin_sqrt = transformed_summary,
    selection_tests = selection_tests,
    .descriptions = c(
      plotted_observations = "Values backing the points and boxes; prediction rows form boxes and adjusted rows form points.",
      distribution_summary_raw = "Counts and quantiles on the raw proportion scale by plotted cohort.",
      distribution_summary_asin_sqrt = "Counts and quantiles on the arcsine-square-root scale used in the panel.",
      selection_tests = "Hypothesis-test results and ranking used to choose displayed cell types."
    )
  )
}

# --- Fig 1: adjusted-count UMAPs ---------------------------------------------
if (build_figure("1")) {
tissue_plot_path <- first_existing(
  file.path(tissue_dir, "Plot___umap_tissue_fig1.rds"),
  file.path(tissue_dir, "Plot___umap_tissue_adj_fig1.rds"),
  file.path(tissue_dir, "Plot___umap_tissue_groups_fig1.rds")
)
fig1_tissue <- read_plot_data(tissue_plot_path)

# Older tissue UMAP caches did not retain cell type in the main ggplot data.
# The companion cache has identical coordinates, so restore that metadata.
cell_type_path <- file.path(tissue_dir, "Plot___umap_cell_type_for_tissue_fig1.rds")
if (!"cell_type" %in% names(fig1_tissue) && file.exists(cell_type_path)) {
  cell_type_data <- read_plot_data(cell_type_path)
  same_coordinates <- nrow(fig1_tissue) == nrow(cell_type_data) &&
    isTRUE(all.equal(fig1_tissue$umap_1, cell_type_data$umap_1)) &&
    isTRUE(all.equal(fig1_tissue$umap_2, cell_type_data$umap_2))
  if (same_coordinates) fig1_tissue$cell_type <- cell_type_data$cell_type
}

fig1_tissue_counts <- fig1_tissue |>
  count(across(any_of(c("tissue_groups", "cell_type"))), name = "n_points") |>
  as.data.frame()

fig1_age <- read_plot_data(file.path(tissue_dir, "Plot___umap_age_fig1.rds"))
fig1_age_counts <- fig1_age |>
  count(across(any_of(c("age_decade", "tissue_groups", "cell_type"))), name = "n_points") |>
  as.data.frame()
fig1_age_summary <- summarise_distribution(
  fig1_age,
  group_cols = intersect(c("tissue_groups", "cell_type"), names(fig1_age)),
  value_col = "age_years",
  value_label = "age_years"
)

fig1_panels <- list(
  Fig1_A_tissue_UMAP = bind_panel_records(
    plotted_coordinates = fig1_tissue,
    cohort_counts = fig1_tissue_counts
  ),
  Fig1_B_age_UMAP = bind_panel_records(
    plotted_coordinates = fig1_age,
    age_decade_counts = fig1_age_counts,
    cohort_age_summary = fig1_age_summary
  )
)
fig1_sheets <- write_figure_panels(fig1_panels, "SupplData_AgeClock_Fig1.xlsx")
rm(list = ls(pattern = "^fig1_"))
invisible(gc())
}

# --- Fig 2 -------------------------------------------------------------------
if (build_figure("2")) {
fig2_sample <- read_plot_data(file.path(fig_dir, "Fig2_Plot__sample_summary.rds"))
fig2_volcano <- read_plot_data(file.path(fig_dir, "Fig2_Plot__volcano.rds")) |>
  mutate(
    test_direction = case_when(c_effect > 0 ~ "increase", c_effect < 0 ~ "decrease", TRUE ~ "zero"),
    passes_selection_FDR = c_FDR < 0.05
  ) |>
  group_by(test_direction) |>
  arrange(c_FDR, .by_group = TRUE) |>
  mutate(direction_FDR_rank = row_number()) |>
  ungroup()

fig2_heat <- read_plot_data(first_existing(
  file.path(fig_dir, "Fig2_Plot__Heatmap_Binarised_age_v2.rds"),
  file.path(fig_dir, "Fig2_Plot__Heatmap_Binarised_age.rds")
)) |>
  mutate(significant = c_FDR < 0.05)

anatogram_full_path <- file.path(age_data_dir, "Fig2_Plot__anatogram_df_full.csv")
anatogram_plot_path <- file.path(age_data_dir, "Fig2_Plot__anatogram_df.csv")
if (!file.exists(anatogram_full_path) || !file.exists(anatogram_plot_path)) {
  stop("Missing compact anatogram statistics under data/age_fig.", call. = FALSE)
}
fig2_anatogram_full <- readr::read_csv(anatogram_full_path, show_col_types = FALSE)
fig2_anatogram_plot <- readr::read_csv(anatogram_plot_path, show_col_types = FALSE) |>
  distinct()

make_fig2_selection <- function(direction, selected_cell_types) {
  fig2_volcano |>
    filter(passes_selection_FDR, test_direction == direction, L3 != "cytotoxic") |>
    arrange(c_FDR) |>
    mutate(
      selection_rank = row_number(),
      selected_for_panel = L3 %in% selected_cell_types
    ) |>
    as.data.frame()
}

fig2_box_up <- read_plot_data(file.path(fig_dir, "Fig2_Plot__boxplot_up_cell_type.rds"))
fig2_box_down <- read_plot_data(file.path(fig_dir, "Fig2_Plot__boxplot_down_cell_type.rds")) |>
  mutate(L3_abbrev = if_else(
    c_FDR > 0 & grepl("\\(FDR: 0\\)$", L3_abbrev),
    paste0(
      sub("\\n\\(FDR:.*$", "", as.character(L3_abbrev)),
      "\n(FDR: ", sub("^0", "", formatC(c_FDR, format = "f", digits = 3)), ")"
    ),
    as.character(L3_abbrev)
  ))

fig2_up_tests <- make_fig2_selection("increase", unique(fig2_box_up$L3))
fig2_down_tests <- make_fig2_selection("decrease", unique(fig2_box_down$L3))

fig2_line_raw <- read_plot_data(file.path(fig_dir, "Fig2_Plot__plasma_line_chart.rds"))
fig2_line <- fig2_line_raw |>
  mutate(age_decade = as.character(age_decade)) |>
  group_by(tissue_groups, tissue_groups_short, cell_type_unified_ensemble,
           age_decade, significant) |>
  summarise(
    n_source_rows = n(),
    n_samples = n_distinct(sample_id),
    proportion_mean = first(proportion_mean),
    proportion_lower = first(proportion_lower),
    proportion_upper = first(proportion_upper),
    .groups = "drop"
  ) |>
  as.data.frame()
fig2_line_summary <- summarise_trajectory(
  fig2_line, c("tissue_groups", "tissue_groups_short", "significant"),
  "age_decade", "proportion_mean"
)
fig2_plasma_tests <- fig2_heat |>
  filter(cell_type_unified_ensemble == "plasma", tissue %in% fig2_line$tissue_groups) |>
  as.data.frame()

fig2_panels <- list(
  Fig2_A_age_range = bind_panel_records(plotted_summary = fig2_sample),
  Fig2_B_age_volcano = bind_panel_records(age_composition_tests = fig2_volcano),
  Fig2_C_anatogram = bind_panel_records(
    plotted_organ_statistics = fig2_anatogram_plot,
    lineage_and_contrast_tests = fig2_anatogram_full,
    .descriptions = c(
      plotted_organ_statistics = "Organ-level effects used for the anatogram; polygon coordinates are intentionally omitted.",
      lineage_and_contrast_tests = "Full lineage posterior summaries and myeloid-versus-lymphoid contrast statistics."
    )
  ),
  Fig2_D_age_heatmap = bind_panel_records(age_effect_tests = fig2_heat),
  Fig2_E_boxplots = combine_panel_records(
    `age-associated enrichment` = boxplot_panel(
      fig2_box_up, fig2_up_tests, "age_decade"
    ),
    `age-associated depletion` = boxplot_panel(
      fig2_box_down, fig2_down_tests, "age_decade"
    )
  ),
  Fig2_F_plasma_trajectory = bind_panel_records(
    plotted_estimates = fig2_line,
    trajectory_summary = fig2_line_summary,
    age_effect_tests = fig2_plasma_tests,
    .descriptions = c(
      plotted_estimates = "Unique tissue-by-age estimates; n_source_rows records duplicated cache rows and n_samples records the cohort size.",
      trajectory_summary = "Start, end, change and extrema across the plotted age decades.",
      age_effect_tests = "Plasma-cell tissue age-effect tests that determine highlighted trajectories."
    )
  )
)
fig2_sheets <- write_figure_panels(fig2_panels, "SupplData_AgeClock_Fig2.xlsx")
rm(list = ls(pattern = "^fig2_"))
invisible(gc())
}

# --- Manuscript Fig 4 (legacy cache prefix Fig3) -----------------------------
if (build_figure("4")) {
fig3_linear <- read_plot_data(file.path(fig_dir, "Fig3_Plot__linear_distance_all.rds"))
fig3_linear_summary <- summarise_trajectory(
  fig3_linear, "version", "age_decade_continuous", "cum_dist_scaled"
)

fig3_sccomp_selected <- read_plot_data(file.path(fig_dir, "Fig3_Plot__sccomp_1d_colored.rds")) |>
  mutate(selected_max_abs_effect = TRUE)
fig3_sccomp_all_path <- file.path(age_data_dir, "Fig3_age_split_tests_compact.rds")
fig3_sccomp_all <- if (file.exists(fig3_sccomp_all_path)) {
  readRDS(fig3_sccomp_all_path)
} else {
  warning("Missing compact Fig3 age-split tests; exporting selected tests only.", call. = FALSE)
  NULL
}

fig3_upset <- read_plot_data(file.path(fig_dir, "Fig3_Plot__upset.rds")) |>
  mutate(ct_set_label = vapply(ct_set, paste, collapse = "; ", FUN.VALUE = character(1)))
fig3_upset_summary <- fig3_upset |>
  count(ct_set_label, n_ct, name = "n_DE_genes") |>
  arrange(desc(n_DE_genes), desc(n_ct)) |>
  as.data.frame()

fig3_velocity <- read_plot_data(file.path(fig_dir, "Fig3_Plot__smoothed_age_dist_all.rds")) |>
  select(-any_of(c("data", "peaks")))
fig3_peak_summary <- fig3_velocity |>
  distinct(tissue_groups, tissue_groups_short, version, min_age_decade,
           max_age_decade, main_peak_age, secondary_peak_age, min_main_peak) |>
  arrange(main_peak_age, tissue_groups, version) |>
  as.data.frame()

fig3_clock <- read_plot_data(file.path(fig_dir, "Fig3_Plot__pred_whole.rds")) |>
  mutate(prediction_error_years = predicted_age_scaled - age_years)
fig3_clock_models <- summarise_linear_models(
  fig3_clock, "tissue_groups", "age_years", "predicted_age_scaled"
)
fig3_clock_errors <- summarise_distribution(
  fig3_clock, "tissue_groups", "prediction_error_years", "prediction_error_years"
)

fig3_donors <- readRDS(file.path(age_data_dir, "Fig3_total_pred_summary.rds")) |>
  as.data.frame()
fig3_donor_summary <- fig3_donors |>
  mutate(
    predicted_age_mean_years = predicted_age_mean / 365,
    actual_age_years = actual_age / 365,
    absolute_error_years = abs(predicted_age_mean_years - actual_age_years)
  ) |>
  group_by(donor_id) |>
  summarise(
    actual_age_years = first(actual_age_years),
    n_tissues = n_distinct(tissue_groups),
    n_samples = sum(n_s),
    mean_predicted_age_years = mean(predicted_age_mean_years),
    min_predicted_age_years = min(predicted_age_mean_years),
    max_predicted_age_years = max(predicted_age_mean_years),
    mean_absolute_error_years = mean(absolute_error_years),
    .groups = "drop"
  ) |>
  as.data.frame()

fig3_panels <- list(
  Fig4_A_linear_distance = bind_panel_records(
    plotted_distances = fig3_linear,
    distance_summary = fig3_linear_summary
  ),
  Fig4_B_age_split_effects = bind_panel_records(
    selected_plotted_tests = fig3_sccomp_selected,
    all_age_split_tests = fig3_sccomp_all
  ),
  Fig4_C_shared_DE_genes = bind_panel_records(
    gene_membership = fig3_upset,
    intersection_counts = fig3_upset_summary
  ),
  Fig4_D_ageing_velocity = bind_panel_records(
    plotted_velocity_curves = fig3_velocity,
    peak_summary = fig3_peak_summary
  ),
  Fig4_E_age_clock = bind_panel_records(
    plotted_predictions = fig3_clock,
    tissue_linear_models = fig3_clock_models,
    prediction_error_summary = fig3_clock_errors,
    .descriptions = c(
      plotted_predictions = "Actual and predicted ages used for points and tissue-specific linear smooths.",
      tissue_linear_models = "Linear-model statistics matching the geom_smooth fit within each tissue.",
      prediction_error_summary = "Prediction error (predicted minus actual age) distribution by tissue."
    )
  ),
  Fig4_F_donor_predictions = bind_panel_records(
    plotted_donor_tissue_summary = fig3_donors,
    donor_summary = fig3_donor_summary
  )
)
fig3_sheets <- write_figure_panels(fig3_panels, "SupplData_AgeClock_Fig4.xlsx")
rm(list = ls(pattern = "^fig3_"))
invisible(gc())
}

# --- Manuscript Fig 5 (legacy cache prefix Fig4) -----------------------------
if (build_figure("5")) {
fig4_counts <- read_plot_data(first_existing(
  file.path(fig_dir, "Fig4_Plot__age_sex_ethnicity_tissue_bar_consensus.rds"),
  file.path(fig_dir, "Fig4_Plot__age_sex_ethnicity_tissue_bar.rds")
))
fig5_all_counts_path <- file.path(fig_dir, "Fig5_A_donor_counts_all_tissues.rds")
if (!file.exists(fig5_all_counts_path)) {
  stop(
    "Missing corrected all-tissue Fig. 5A donor counts: ", fig5_all_counts_path,
    ". Rebuild panel A with REBUILD_FIG5_PANEL_A=true before exporting.",
    call. = FALSE
  )
}
fig4_counts_all <- readRDS(fig5_all_counts_path) |>
  as.data.frame()
fig4_counts_summary <- fig4_counts_all |>
  group_by(tissue_groups, tissue_label, sex) |>
  summarise(donor_count = sum(donor_count), .groups = "drop") |>
  as.data.frame()

fig4_box_red <- read_plot_data(file.path(fig_dir, "Fig4_Plot__boxplot_uniformly_red.rds")) |>
  mutate(selection_direction = "positive male-minus-female effect")
fig4_box_blue <- read_plot_data(file.path(fig_dir, "Fig4_Plot__boxplot_uniformly_blue.rds")) |>
  mutate(selection_direction = "negative male-minus-female effect")
fig4_box <- bind_rows(fig4_box_red, fig4_box_blue)

fig4_global_test_path <- file.path(age_data_dir, "Fig4_global_sex_tests_compact.rds")
if (!file.exists(fig4_global_test_path)) {
  stop("Missing Fig4 global sex-test cache: ", fig4_global_test_path, call. = FALSE)
}
fig4_global_tests <- readRDS(fig4_global_test_path) |>
  mutate(significant = c_FDR < 0.05)
fig4_best_global <- fig4_global_tests |>
  group_by(cell_type_unified_ensemble) |>
  arrange(c_FDR, desc(abs(c_effect)), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(
    selection_direction = case_when(
      c_effect > 0 ~ "positive male-minus-female effect",
      c_effect < 0 ~ "negative male-minus-female effect",
      TRUE ~ "zero effect"
    )
  ) |>
  group_by(selection_direction) |>
  mutate(
    # Match the panel-selection code: rank effects only among candidates that
    # pass FDR < 0.05, separately for positive and negative directions.
    direction_effect_rank = min_rank(if_else(
      significant,
      if_else(
        selection_direction == "positive male-minus-female effect",
        -c_effect,
        c_effect
      ),
      NA_real_
    )),
    selected_for_panel = cell_type_unified_ensemble %in% unique(fig4_box$L3)
  ) |>
  ungroup() |>
  as.data.frame()

fig4_heat <- read_plot_data(file.path(fig_dir, "Fig4_Plot__heatmap_male_vs_female.rds")) |>
  mutate(significant = c_FDR < 0.05)

fig4_line_raw <- read_plot_data(file.path(fig_dir, "Fig4_Plot__line_chart_plasma.rds"))
fig4_line <- fig4_line_raw |>
  mutate(Age = as.character(Age)) |>
  group_by(Tissue, tissue_groups_short, cell_type_unified_ensemble, Age, Sex) |>
  summarise(
    n_source_rows = n(),
    n_samples = n_distinct(sample_id),
    proportion_mean = first(proportion_mean),
    proportion_lower = first(proportion_lower),
    proportion_upper = first(proportion_upper),
    .groups = "drop"
  ) |>
  as.data.frame()
fig4_line_summary <- summarise_trajectory(
  fig4_line, c("Tissue", "tissue_groups_short", "Sex"), "Age", "proportion_mean"
)
fig4_line_sex_difference <- fig4_line |>
  select(Tissue, tissue_groups_short, Age, Sex, proportion_mean) |>
  pivot_wider(names_from = Sex, values_from = proportion_mean) |>
  mutate(male_minus_female_prediction = male - female) |>
  as.data.frame()
fig4_plasma_tests <- fig4_heat |>
  filter(cell_type_unified_ensemble == "plasma", Tissue %in% fig4_line$Tissue) |>
  as.data.frame()

fig4_box_raw <- add_asin_sqrt(fig4_box) |>
  mutate(plot_role = if_else(
    version == "prediction", "boxplot distribution", "adjusted-data points"
  ))
fig4_box_summary_raw <- summarise_distribution(
  fig4_box_raw,
  c("L3", "L3_abbrev", "sex", "version", "plot_role", "selection_direction", "c_FDR"),
  "proportion_mean", "proportion"
)
fig4_box_summary_transformed <- summarise_distribution(
  fig4_box_raw,
  c("L3", "L3_abbrev", "sex", "version", "plot_role", "selection_direction", "c_FDR"),
  "proportion_asin_sqrt", "asin_sqrt_proportion"
)

fig4_panels <- list(
  Fig5_A_donor_counts = bind_panel_records(
    plotted_stacked_counts = fig4_counts,
    all_tissue_donor_counts = fig4_counts_all,
    all_tissue_sex_totals = fig4_counts_summary,
    .descriptions = c(
      plotted_stacked_counts = "Donor counts for the 25 tissues displayed in Fig. 5A; signed_donor_count is used only to make the diverging plot.",
      all_tissue_donor_counts = "Unique donor counts for every represented tissue, age group and sex; donors with multiple samples or assays are counted once per tissue.",
      all_tissue_sex_totals = "Unique donor totals by tissue and sex, summed across the three displayed age groups."
    )
  ),
  Fig5_B_sex_boxplots = bind_panel_records(
    plotted_observations = fig4_box_raw,
    distribution_summary_raw = fig4_box_summary_raw,
    distribution_summary_asin_sqrt = fig4_box_summary_transformed,
    selected_cell_type_tests = fig4_best_global,
    all_global_sex_tests_by_age = fig4_global_tests
  ),
  Fig5_C_sex_heatmap = bind_panel_records(tissue_sex_effect_tests = fig4_heat),
  Fig5_D_plasma_trajectories = bind_panel_records(
    plotted_estimates = fig4_line,
    trajectory_summary = fig4_line_summary,
    predicted_sex_differences = fig4_line_sex_difference,
    plasma_tissue_sex_tests = fig4_plasma_tests
  )
)
fig4_sheets <- write_figure_panels(fig4_panels, "SupplData_AgeClock_Fig5.xlsx")
rm(list = ls(pattern = "^fig4_"))
invisible(gc())
}

message("Done. Split-table per-figure workbooks written to: ", suppl_dir)
message("Build the master workbook in a fresh process with report/R/build_master_supplementary_workbook.R")
