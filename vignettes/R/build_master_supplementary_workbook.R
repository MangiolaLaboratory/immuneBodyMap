#!/usr/bin/env Rscript
# Assemble the per-figure workbooks in a fresh R process to keep peak memory
# below that required by figure-cache extraction and workbook writing together.

root <- Sys.getenv(
  "AGE_CLOCK_ROOT",
  unset = Sys.getenv("AGE_CLOCK_ARCHIVE", unset = "")
)
if (!nzchar(root)) {
  stop("Set AGE_CLOCK_ROOT or AGE_CLOCK_ARCHIVE to rebuild the master workbook.", call. = FALSE)
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

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("build_master_supplementary_workbook.R requires readxl.", call. = FALSE)
}

suppl_dir <- file.path(root, "Supplementary_files")
fig_dir <- file.path(root, "report", "Fig_files")
figure_workbooks <- file.path(
  suppl_dir,
  c(
    "SupplData_AgeClock_Fig1.xlsx",
    "SupplData_AgeClock_Fig2.xlsx",
    "SupplData_AgeClock_Fig4.xlsx",
    "SupplData_AgeClock_Fig5.xlsx"
  )
)

master <- list()
used <- character()
append_workbook <- function(path, name_prefix = "") {
  if (!file.exists(path)) stop("Missing workbook: ", path, call. = FALSE)
  for (sheet in readxl::excel_sheets(path)) {
    sheet_name <- sanitize_sheet_name(paste0(name_prefix, sheet), used)
    used <<- c(used, sheet_name)
    master[[sheet_name]] <<- sanitize_for_excel(
      as.data.frame(readxl::read_xlsx(
        path,
        sheet = sheet,
        guess_max = 1000000L,
        .name_repair = "minimal"
      ))
    )
  }
  invisible(gc())
}

for (path in figure_workbooks) append_workbook(path)

venn_xlsx <- file.path(fig_dir, "Fig4_venn_supplementary_tables.xlsx")
if (file.exists(venn_xlsx)) append_workbook(venn_xlsx, "Fig5_venn_")

csvs <- c(
  S1_multitissue_donor_predictions = "SupplData_AgeClock_S1_multitissue_donor_predictions.csv",
  S2_pertissue_accuracy_multitissue = "SupplData_AgeClock_S2_pertissue_accuracy_multitissue.csv",
  S3_per_donor_summary = "SupplData_AgeClock_S3_per_donor_summary.csv",
  S4_overall_pertissue_accuracy = "SupplData_AgeClock_S4_overall_pertissue_accuracy.csv"
)
for (name in names(csvs)) {
  path <- file.path(suppl_dir, csvs[[name]])
  if (!file.exists(path)) {
    warning("Missing CSV for master workbook: ", path, call. = FALSE)
    next
  }
  sheet_name <- sanitize_sheet_name(name, used)
  used <- c(used, sheet_name)
  master[[sheet_name]] <- sanitize_for_excel(
    utils::read.csv(path, check.names = FALSE)
  )
}

output_path <- file.path(suppl_dir, "SupplData_AgeClock_master.xlsx")
writexl::write_xlsx(master, path = output_path)
message("Wrote ", length(master), " sheet(s) to ", output_path)
