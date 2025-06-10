
library(targets)
setwd('/home/zhanchen/From_scratch/HPC_sccomp/sccomp_on_cellNexus_1_0_10_2')

tar_make(
  # callr_function = NULL,
  script = "/scratch/pawsey1192/zhanchen/HPC_sccomp/sccomp_on_cellNexus_1_0_10_2/_targets.R", 
  store = "/scratch/pawsey1192/zhanchen/HPC_sccomp/sccomp_on_cellNexus_1_0_10_2/_targets", 
  reporter = "verbose" #, callr_function = NULL
) 
  
  
  
  



system("rclone copy /scratch/pawsey1192/zhanchen/HPC_sccomp/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins___L3___disease_TRUE___immune_only_TRUE.rds UofA_Box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_10_2/")



tar_meta(store = "/scratch/pawsey1192/zhanchen/HPC_sccomp/sccomp_on_cellNexus_1_0_10_2/_targets") |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  dplyr::select(name, error)




library(tidyverse)
library(sccomp)
library(magrittr)
library(glue)
library(forcats)
library(stringr)

library(arrow)
library(dplyr)
library(duckdb)

library(targets)

x = tar_read(input_relative, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/_targets")


tar_workspace(estimates_ae6e35523d730ab3, 
              script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/_targets.R", 
              store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/_targets"
              )

tar_meta(starts_with("estimates_"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/_targets")

# Age proportion prediciton
estimates_age_bins |> 
  sccomp_predict(
    formula_composition = ~ 1 + age_bin*sex + (1 + age_bin*sex | tissue_groups), 
    number_of_draws = 100,
    summary_instead_of_draws = TRUE
  ) |> 
  mutate(age_bin = factor(
    age_bin, 
    c("Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior"),
    ordered = TRUE
  )) |> 
  mutate(age_bin_numeric = age_bin |> as.integer())  |> 
  saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/prediction_age_bins.rds")

system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/prediction_age_bins.rds UofA_Box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

tar_read(formula_df, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/_targets")


# For Hong
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins___L3.rds")
estimate_age_bins = estimate_age_bins |> dplyr::select(-count_data)
attr(estimate_age_bins, "fit") = NULL
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins_effect_tibble_only.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins_effect_tibble_only.rds UofA_Box:/immune_map_disease/")

# Save fit
library(magrittr)
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins.rds")
estimate_age_bins |> attr("fit") %$% save_object(file = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins_FIT_FOR_PORTABILITY.rds") 
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins_FIT_FOR_PORTABILITY.rds UofA_Box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
estimate_age_bins |> attr("fit") = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins_FIT_FOR_PORTABILITY.rds")
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins.rds UofA_Box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")


# estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/estimates_age_bins.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_10_2/21_11_2024_sccomp_archive_before_factor_ordering/estimates_age_bins.rds UofA_Box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

# Benchmark
tic()
estimate_age_bins |> sccomp_test(contrasts = c(  "respiratory system" = "sexmale + `sexmale___respiratory system`",
                                                 "blood" = "sexmale + `sexmale___blood`"))
toc()


estimate_age_bins |> 
  sccomp_test()
