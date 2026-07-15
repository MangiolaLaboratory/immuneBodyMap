#!/usr/bin/env Rscript
# Export the metadata recoding rules used by edit_covariates() as a
# manuscript-facing supplementary workbook.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

if (!requireNamespace("writexl", quietly = TRUE)) {
  stop("The writexl package is required to create the supplementary workbook.")
}

# Prefer publication package roots; fall back to AGE_CLOCK_ARCHIVE for developers.
source(file.path(
  Sys.getenv(
    "IMMUNE_HEALTHY_BODY_MAP_ROOT",
    unset = normalizePath(file.path("..", ".."), mustWork = FALSE)
  ),
  "vignettes", "R", "paths.R"
))

project_root <- Sys.getenv(
  "AGE_CLOCK_ROOT",
  unset = Sys.getenv(
    "AGE_CLOCK_ARCHIVE",
    unset = publication_repo_root()
  )
)
edit_covariates_path <- Sys.getenv(
  "EDIT_COVARIATES_PATH",
  unset = path_edit_covariates()
)
disease_override_path <- Sys.getenv("DISEASE_GROUP_OVERRIDE_PATH", unset = "")
if (!nzchar(disease_override_path)) {
  candidates <- c(
    file.path(publication_data_root(), "processed", "disease_data_grouped_further.csv"),
    file.path(publication_data_root(), "zenodo_release", "metadata", "disease_data_grouped_further.csv")
  )
  disease_override_path <- candidates[file.exists(candidates)][1]
  if (is.na(disease_override_path)) {
    stop(
      "Set DISEASE_GROUP_OVERRIDE_PATH to disease_data_grouped_further.csv ",
      "(see manifests/external_large_files.csv / Zenodo metadata).",
      call. = FALSE
    )
  }
}
output_path <- file.path(
  path_source_tables(),
  "SupplData_MultiTissueImmuneAging_MetadataMappings.xlsx"
)

stopifnot(file.exists(edit_covariates_path), file.exists(disease_override_path))
source(edit_covariates_path)

# The mapping objects are local to edit_covariates(). Evaluate only their
# assignment expressions, never the final data-mutating pipeline.
function_body <- body(edit_covariates)
mapping_environment <- new.env(parent = globalenv())
for (expression_index in 2:5) {
  eval(function_body[[expression_index]], envir = mapping_environment)
}

tissue_raw <- mapping_environment$tissue_grouped
ethnicity_raw <- mapping_environment$ethnicity_grouped
assay_raw <- mapping_environment$assay_data_grouped
disease_default <- mapping_environment$disease_data_grouped
disease_override <- readr::read_csv(
  disease_override_path,
  show_col_types = FALSE
)

# Reproduce the override direction in edit_covariates(): the base table is on
# the left, so labels found only in the override file are not joinable by the
# function and are shown explicitly as such in the audit table.
disease_effective <- disease_default |>
  left_join(
    disease_override |>
      rename(disease_groups_further = disease_groups),
    by = "disease"
  ) |>
  mutate(disease_groups = if_else(
    !is.na(disease_groups_further), disease_groups_further, disease_groups
  )) |>
  select(disease, disease_groups)

tissue_counts <- tissue_raw |>
  count(tissue_groups, name = "source_labels_in_group")
tissue_groups <- tissue_raw |>
  transmute(
    source_field = "tissue",
    source_label = tissue,
    derived_field = "tissue_groups",
    derived_group = tissue_groups,
    post_join_override = tissue %in% c(
      "gallbladder", "pancreas", "exocrine pancreas",
      "liver", "caudate lobe of liver", "hepatic cecum"
    ),
    mapping_method = "case-sensitive exact left join"
  ) |>
  left_join(tissue_counts, by = c("derived_group" = "tissue_groups")) |>
  arrange(derived_group, source_label)

assay_counts <- assay_raw |>
  count(assay_groups, name = "source_labels_in_group")
assay_groups <- assay_raw |>
  transmute(
    source_field = "assay",
    source_label = assay,
    derived_field = "assay_groups",
    derived_group = assay_groups,
    mapping_method = "case-sensitive exact left join"
  ) |>
  left_join(assay_counts, by = c("derived_group" = "assay_groups")) |>
  arrange(derived_group, source_label)

ethnicity_counts <- ethnicity_raw |>
  count(ethnicity_groups, name = "source_labels_in_group")
ethnicity_groups <- ethnicity_raw |>
  transmute(
    source_field = "self_reported_ethnicity",
    source_label = self_reported_ethnicity,
    derived_field = "ethnicity_groups",
    derived_group = ethnicity_groups,
    mapping_method = "case-sensitive exact left join"
  ) |>
  left_join(ethnicity_counts, by = c("derived_group" = "ethnicity_groups")) |>
  arrange(derived_group, source_label)

disease_groups <- full_join(
  disease_default |>
    rename(default_group = disease_groups) |>
    mutate(present_in_base_mapping = TRUE),
  disease_override |>
    rename(external_override_group = disease_groups) |>
    mutate(present_in_external_override = TRUE),
  by = "disease"
) |>
  mutate(
    present_in_base_mapping = replace_na(present_in_base_mapping, FALSE),
    present_in_external_override = replace_na(
      present_in_external_override, FALSE
    ),
    override_applied_by_function =
      present_in_base_mapping & present_in_external_override,
    effective_group_before_tissue_suffix = case_when(
      present_in_base_mapping & present_in_external_override ~
        external_override_group,
      present_in_base_mapping ~ default_group,
      TRUE ~ NA_character_
    ),
    mapping_status = case_when(
      override_applied_by_function ~ "external override applied",
      present_in_base_mapping ~ "default mapping retained",
      TRUE ~ paste(
        "override-only label; not reachable because the function left-joins",
        "the override onto the base mapping"
      )
    ),
    final_model_value_rule = paste(
      "Normal remains Normal; every other mapped group is suffixed with",
      "'_' plus tissue_groups"
    )
  ) |>
  transmute(
    source_field = "disease",
    source_label = disease,
    default_group,
    external_override_group,
    effective_group_before_tissue_suffix,
    present_in_base_mapping,
    present_in_external_override,
    override_applied_by_function,
    mapping_status,
    final_model_value_rule
  ) |>
  arrange(effective_group_before_tissue_suffix, source_label)

age_bins <- tribble(
  ~record_type, ~rule_order, ~input_field, ~output_field,
  ~lower_years, ~lower_inclusive, ~upper_years, ~upper_inclusive,
  ~output_value, ~code_rule, ~notes,
  "inclusion", 1L, "age_days", "row retained", 1, FALSE, NA, NA,
  "TRUE", "age_days > 365", "Exactly 365 days and younger are excluded",
  "derived value", 2L, "age_days", "age_years", NA, NA, NA, NA,
  NA, "age_days / 365", "A 365-day year is used",
  "age_bin", 3L, "age_years", "age_bin", 1, FALSE, 3, FALSE,
  "Infancy", "age_years < 3", "Lower bound follows the study inclusion filter",
  "age_bin", 4L, "age_years", "age_bin", 3, TRUE, 12, FALSE,
  "Childhood", "age_years < 12", NA,
  "age_bin", 5L, "age_years", "age_bin", 12, TRUE, 20, FALSE,
  "Adolescence", "age_years < 20", NA,
  "age_bin", 6L, "age_years", "age_bin", 20, TRUE, 40, FALSE,
  "Young Adulthood", "age_years < 40", NA,
  "age_bin", 7L, "age_years", "age_bin", 40, TRUE, 50, FALSE,
  "Middle Age", "age_years < 50", NA,
  "age_bin", 8L, "age_years", "age_bin", 50, TRUE, 60, FALSE,
  "Senior_50", "age_years < 60", NA,
  "age_bin", 9L, "age_years", "age_bin", 60, TRUE, 70, FALSE,
  "Senior_60", "age_years < 70", NA,
  "age_bin", 10L, "age_years", "age_bin", 70, TRUE, NA, NA,
  "Senior_70", "age_years >= 70", NA,
  "scaled value", 21L, "age_days", "age_days_scaled", NA, NA, NA, NA,
  NA, "scale(age_days, center = 50 * 365)",
  paste(
    "The default scale=TRUE makes this cohort-dependent; the numeric centre",
    "is 50 years, although the source-code comment says adolescence"
  )
)

age_decades <- tibble(
  record_type = "age_decade",
  rule_order = 10L + seq_len(10L),
  input_field = "age_years",
  output_field = "age_decade",
  lower_years = (0:9) * 10,
  lower_inclusive = FALSE,
  upper_years = (1:10) * 10,
  upper_inclusive = TRUE,
  output_value = as.character(1:10),
  code_rule = "as.character(as.integer(ceiling(age_years / 10)))",
  notes = "The formula is not capped; values above 100 years would create decade 11+"
)
age_groups <- bind_rows(age_bins, age_decades) |>
  arrange(rule_order)

derived_fields <- tribble(
  ~step_order, ~source_field, ~derived_field, ~transformation, ~notes,
  1L, "tissue", "tissue_groups",
  "Remove any existing tissue_groups column, then exact-left-join Tissue_groups",
  "Unmatched tissue labels remain missing at this step",
  2L, "assay", "assay_groups",
  "Exact-left-join Assay_groups",
  "Unmatched assay labels remain missing",
  3L, "disease", "disease_groups",
  "Exact-left-join the effective Disease_groups mapping",
  "The external override refines only labels already present in the base mapping",
  4L, "disease_groups + tissue_groups", "disease_groups",
  "paste(disease_groups, tissue_groups, sep = '_')",
  "Normal_<tissue> is reset to Normal; other diseases become tissue-specific",
  5L, "tissue", "tissue_groups",
  "Explicitly reset gallbladder, pancreas/exocrine pancreas, and liver labels",
  "These overrides currently agree with the primary tissue mapping",
  6L, "sex", "sex",
  "Replace missing values with 'unknown'",
  "Non-missing labels are retained unchanged",
  7L, "age_days", "row inclusion",
  "Retain age_days > 365",
  "Applied before all remaining age and ethnicity transformations",
  8L, "age_days", "age_years; age_bin; age_decade; age_days_scaled",
  "Apply the rules documented in Age_groups",
  "age_days_scaled depends on the full input cohort",
  9L, "self_reported_ethnicity", "ethnicity_groups",
  "Exact-left-join Ethnicity_groups",
  "Unmatched ethnicity labels remain missing",
  10L, "all mapped fields", "output tibble",
  "Return as a tibble",
  "Mappings are case-sensitive and do not trim whitespace or normalise punctuation"
)

mapping_summary <- tribble(
  ~metadata_field, ~source_column, ~derived_column, ~source_labels_or_rules,
  ~possible_groups, ~implementation_note,
  "Tissue", "tissue", "tissue_groups", nrow(tissue_groups),
  n_distinct(tissue_groups$derived_group),
  "38 possible code-defined groups; the manuscript reports 36 represented groups after analysis inclusion",
  "Assay", "assay", "assay_groups", nrow(assay_groups),
  n_distinct(assay_groups$derived_group), "Exact left join",
  "Self-reported ethnicity", "self_reported_ethnicity", "ethnicity_groups",
  nrow(ethnicity_groups), n_distinct(ethnicity_groups$derived_group),
  "Exact left join; this is label harmonisation, not genetic-ancestry inference",
  "Disease", "disease", "disease_groups", nrow(disease_groups),
  n_distinct(na.omit(disease_groups$effective_group_before_tissue_suffix)),
  "External overrides refine the base groups; non-Normal values are then made tissue-specific",
  "Age", "age_days", "age_years; age_bin; age_decade; age_days_scaled",
  nrow(age_groups), NA_integer_, "Formula- and interval-based derivations",
  "Sex", "sex", "sex", 1L, NA_integer_, "Only missing values are changed to 'unknown'"
)

source_md5 <- unname(tools::md5sum(edit_covariates_path))
override_md5 <- unname(tools::md5sum(disease_override_path))
readme <- tribble(
  ~item, ~value,
  "Workbook purpose",
  paste(
    "Metadata mappings and transformations used for the manuscript",
    "A multi-tissue map of human immune aging_2026"
  ),
  "Methods section", "Data source, metadata harmonisation and sample inclusion",
  "Methods wording",
  paste(
    "For the present study, selected metadata fields were further transformed",
    "or grouped to improve model stability and interpretability."
  ),
  "Manuscript document",
  "https://docs.google.com/document/d/1Ri2Dyr4Rhv1rvgIXeUm5NI27IrdPKJS1jI2J-SJjk9I",
  "Mapping code", edit_covariates_path,
  "Mapping code MD5", source_md5,
  "Disease override source", disease_override_path,
  "Disease override MD5", override_md5,
  "Generated on", as.character(Sys.Date()),
  "Join semantics",
  paste(
    "All categorical mappings are case-sensitive exact left joins.",
    "Labels absent from a mapping do not receive a grouped value."
  ),
  "Tissue group scope",
  paste(
    "The function defines 38 possible tissue groups from 281 source labels;",
    "the manuscript reports 36 groups represented in the analysed atlas."
  ),
  "Disease group scope",
  paste(
    "Disease_groups lists the pre-suffix mapping. In model-ready output,",
    "non-Normal disease groups are suffixed by tissue group."
  ),
  "Age scaling caution",
  paste(
    "age_days_scaled uses scale() with a 50-year centre and default scaling,",
    "so exact scaled values depend on the input cohort."
  )
)

# Validation against the effective mapping reproduced from the source function.
stopifnot(
  nrow(tissue_groups) == 281L,
  n_distinct(tissue_groups$derived_group) == 38L,
  !anyDuplicated(tissue_groups$source_label),
  nrow(assay_groups) == 24L,
  !anyDuplicated(assay_groups$source_label),
  nrow(ethnicity_groups) == 31L,
  !anyDuplicated(ethnicity_groups$source_label),
  nrow(disease_effective) == 108L,
  !anyDuplicated(disease_effective$disease)
)

disease_check <- disease_groups |>
  filter(present_in_base_mapping) |>
  select(
    disease = source_label,
    disease_groups = effective_group_before_tissue_suffix
  ) |>
  arrange(disease)
stopifnot(identical(
  disease_check,
  disease_effective |> arrange(disease)
))

sheets <- list(
  README = readme,
  Mapping_summary = mapping_summary,
  Tissue_groups = tissue_groups,
  Assay_groups = assay_groups,
  Ethnicity_groups = ethnicity_groups,
  Disease_groups = disease_groups,
  Age_groups = age_groups,
  Derived_fields = derived_fields
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
writexl::write_xlsx(sheets, output_path)
message("Wrote metadata mapping workbook to: ", output_path)
