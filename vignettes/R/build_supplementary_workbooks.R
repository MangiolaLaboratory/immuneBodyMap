#!/usr/bin/env Rscript
# Build SI workbooks. Prefers AGE_CLOCK_ARCHIVE / AGE_CLOCK_ROOT when set so
# that the Age_Clock layout can still be used; otherwise documents publication layout.
repo <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
if (!nzchar(repo)) {
  repo <- normalizePath(file.path("..", ".."), mustWork = FALSE)
}
root <- Sys.getenv(
  "AGE_CLOCK_ROOT",
  unset = Sys.getenv("AGE_CLOCK_ARCHIVE", unset = "")
)
if (!nzchar(root)) {
  message(
    "AGE_CLOCK_ROOT / AGE_CLOCK_ARCHIVE unset. Publication SI tables are already\n",
    "shipped under vignettes/data/source_tables/. To rebuild from Age_Clock plot\n",
    "caches, set AGE_CLOCK_ARCHIVE to the analysis tree."
  )
  quit(save = "no", status = 0L)
}
active_builder <- file.path(root, "report", "R", "build_panel_supplementary_workbooks.R")
if (!file.exists(active_builder)) {
  active_builder <- file.path(repo, "vignettes", "R", "build_panel_supplementary_workbooks.R")
}
master_builder <- file.path(root, "report", "R", "build_master_supplementary_workbook.R")
if (!file.exists(master_builder)) {
  master_builder <- file.path(repo, "vignettes", "R", "build_master_supplementary_workbook.R")
}
selected_figure <- Sys.getenv("SUPPLEMENTARY_FIGURE", "all")

if (selected_figure == "all") {
  for (figure in c("1", "2", "4", "5")) {
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      active_builder,
      env = c(
        paste0("SUPPLEMENTARY_FIGURE=", figure),
        paste0("AGE_CLOCK_ROOT=", root),
        paste0("IMMUNE_HEALTHY_BODY_MAP_ROOT=", repo)
      )
    )
    if (!identical(status, 0L)) {
      stop("Supplementary workbook build failed for Fig. ", figure, call. = FALSE)
    }
  }
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    master_builder,
    env = c(paste0("AGE_CLOCK_ROOT=", root), paste0("IMMUNE_HEALTHY_BODY_MAP_ROOT=", repo))
  )
  if (!identical(status, 0L)) stop("Master workbook build failed.", call. = FALSE)
} else {
  Sys.setenv(AGE_CLOCK_ROOT = root)
  source(active_builder)
}
quit(save = "no", status = 0L)


# --- legacy Age_Clock-only body retained below for reference (unreachable after quit) ---

# Build supplementary Excel workbooks from cached plot RDS objects.
# Run after updating figure caches, or when qmd render is not practical.

suppressPackageStartupMessages({
  if (!requireNamespace("here", quietly = TRUE)) {
    install.packages("here", repos = "https://cloud.r-project.org")
  }
})

source(file.path(here::here(), "report", "R", "export_supplementary_tables.R"))

fig_dir <- file.path(here::here(), "report", "Fig_files")
tissue_dir <- normalizePath(
  "~/lab/chen/HPC_posterior/Tissue_umap",
  mustWork = FALSE
)
suppl_dir <- file.path(here::here(), "Supplementary_files")

read_plot_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

first_existing <- function(...) {
  paths <- list(...)
  for (p in paths) {
    if (!is.null(p) && file.exists(p)) {
      return(p)
    }
  }
  NULL
}

load_named_plots <- function(specs) {
  plots <- list()
  for (nm in names(specs)) {
    path <- specs[[nm]]
    if (is.null(path)) {
      message("Skipping ", nm, " (no RDS found)")
      next
    }
    message("Loading ", basename(path))
    plots[[nm]] <- read_plot_rds(path)
  }
  plots
}

# --- Fig 1 --------------------------------------------------------------------
fig1_specs <- list(
  umap_tissue = first_existing(
    file.path(tissue_dir, "Plot___umap_tissue_fig1.rds"),
    file.path(tissue_dir, "Plot___umap_tissue_groups_fig1.rds")
  ),
  umap_age = file.path(tissue_dir, "Plot___umap_age_fig1.rds")
)

fig1_plots <- load_named_plots(fig1_specs)
fig1_sheets <- plots_to_sheets(fig1_plots, prefix = "Fig1_")
export_supplementary_workbook(
  plots = fig1_plots,
  output_path = file.path(suppl_dir, "SupplData_AgeClock_Fig1.xlsx"),
  prefix = "Fig1_"
)

# --- Fig 2 --------------------------------------------------------------------
fig2_specs <- list(
  sample_summary = file.path(fig_dir, "Fig2_Plot__sample_summary.rds"),
  volcano = file.path(fig_dir, "Fig2_Plot__volcano.rds"),
  anatogram = file.path(fig_dir, "Fig2_Plot__anatogram.rds"),
  heatmap_binarised_age = first_existing(
    file.path(fig_dir, "Fig2_Plot__Heatmap_Binarised_age_v2.rds"),
    file.path(fig_dir, "Fig2_Plot__Heatmap_Binarised_age.rds")
  ),
  boxplot_up_cell_type = file.path(fig_dir, "Fig2_Plot__boxplot_up_cell_type.rds"),
  boxplot_down_cell_type = file.path(fig_dir, "Fig2_Plot__boxplot_down_cell_type.rds"),
  plasma_line_chart = file.path(fig_dir, "Fig2_Plot__plasma_line_chart.rds")
)

fig2_plots <- load_named_plots(fig2_specs)
fig2_sheets <- plots_to_sheets(fig2_plots, prefix = "Fig2_")
export_supplementary_workbook(
  plots = fig2_plots,
  output_path = file.path(suppl_dir, "SupplData_AgeClock_Fig2.xlsx"),
  prefix = "Fig2_"
)

# --- Manuscript Fig 4 (legacy cache prefix Fig3) -----------------------------
fig3_specs <- list(
  linear_distance_all = file.path(fig_dir, "Fig3_Plot__linear_distance_all.rds"),
  smoothed_age_dist_all = file.path(fig_dir, "Fig3_Plot__smoothed_age_dist_all.rds"),
  sccomp_1d_colored = file.path(fig_dir, "Fig3_Plot__sccomp_1d_colored.rds"),
  upset = file.path(fig_dir, "Fig3_Plot__upset.rds"),
  pred_whole = file.path(fig_dir, "Fig3_Plot__pred_whole.rds"),
  total_pred_line_chart = file.path(fig_dir, "Fig3_Plot__total_pred_line_chart.rds")
)

fig3_plots <- load_named_plots(fig3_specs)

# Linear-distance panel is not cached as RDS; rebuild from saved summary if present.
path_total_pred_summary <- file.path(here::here(), "data", "age_fig", "Fig3_total_pred_summary.rds")
fig3_extra <- list()
if (file.exists(path_total_pred_summary)) {
  fig3_extra$total_pred_summary <- readRDS(path_total_pred_summary)
}

fig3_sheets <- plots_to_sheets(fig3_plots, extra_sheets = fig3_extra, prefix = "Fig4_")
export_supplementary_workbook(
  plots = fig3_plots,
  extra_sheets = fig3_extra,
  output_path = file.path(suppl_dir, "SupplData_AgeClock_Fig4.xlsx"),
  prefix = "Fig4_"
)

# --- Manuscript Fig 5 (legacy cache prefix Fig4) -----------------------------
fig4_specs <- list(
  age_sex_ethnicity_tissue_bar = first_existing(
    file.path(fig_dir, "Fig4_Plot__age_sex_ethnicity_tissue_bar_consensus.rds"),
    file.path(fig_dir, "Fig4_Plot__age_sex_ethnicity_tissue_bar.rds")
  ),
  boxplot_uniformly_red = file.path(fig_dir, "Fig4_Plot__boxplot_uniformly_red.rds"),
  boxplot_uniformly_blue = file.path(fig_dir, "Fig4_Plot__boxplot_uniformly_blue.rds"),
  heatmap_male_vs_female = file.path(fig_dir, "Fig4_Plot__heatmap_male_vs_female.rds"),
  line_chart_plasma = file.path(fig_dir, "Fig4_Plot__line_chart_plasma.rds")
)

fig4_plots <- load_named_plots(fig4_specs)
fig4_sheets <- plots_to_sheets(fig4_plots, prefix = "Fig5_")
export_supplementary_workbook(
  plots = fig4_plots,
  output_path = file.path(suppl_dir, "SupplData_AgeClock_Fig5.xlsx"),
  prefix = "Fig5_"
)

# --- Master workbook ----------------------------------------------------------
existing_csvs <- list(
  S1_multitissue_donor_predictions = file.path(
    suppl_dir, "SupplData_AgeClock_S1_multitissue_donor_predictions.csv"
  ),
  S2_pertissue_accuracy_multitissue = file.path(
    suppl_dir, "SupplData_AgeClock_S2_pertissue_accuracy_multitissue.csv"
  ),
  S3_per_donor_summary = file.path(
    suppl_dir, "SupplData_AgeClock_S3_per_donor_summary.csv"
  ),
  S4_overall_pertissue_accuracy = file.path(
    suppl_dir, "SupplData_AgeClock_S4_overall_pertissue_accuracy.csv"
  )
)

# Include Fig4 venn supplementary tables if already generated.
venn_xlsx <- file.path(fig_dir, "Fig4_venn_supplementary_tables.xlsx")
if (file.exists(venn_xlsx) && requireNamespace("readxl", quietly = TRUE)) {
  venn_sheets <- readxl::excel_sheets(venn_xlsx)
  venn_data <- lapply(
    venn_sheets,
    function(s) readxl::read_xlsx(venn_xlsx, sheet = s)
  )
  names(venn_data) <- paste0("Fig5_venn_", venn_sheets)
  fig4_sheets <- c(fig4_sheets, venn_data)
}

build_master_workbook(
  sheet_lists = list(fig1_sheets, fig2_sheets, fig3_sheets, fig4_sheets),
  existing_csvs = existing_csvs,
  output_path = file.path(suppl_dir, "SupplData_AgeClock_master.xlsx")
)

message("Done. Supplementary files written to: ", suppl_dir)
