library(tidyverse)
library(forcats)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

args = commandArgs(trailingOnly=TRUE)
filter_blood = args[[1]]
input_file = args[[2]]
output_file_1 = args[[3]]
output_file_1_blood = args[[4]]
output_file_2 = args[[5]]
output_file_3 = args[[6]]

differential_composition_ethnicity_absolute =
  readRDS(input_file) |>

  # Drop only-immune organs
  filter(!tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone")) |>
  mutate(is_immune = as.character(is_immune)) |>

  filter(development_stage!="unknown") |>

  # Scale days
  mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>

	# Create one-to-many grouping for multilevel modelling
  unite("tissue_harmonised_ethnicity", c(tissue_harmonised , ethnicity_simplified), remove = FALSE) |>

  # Estimate
  sccomp_estimate(
  	formula_composition = ~ age_days*sex + disease + ethnicity_simplified + assay_simplified + disease + group + (1 + age_days + sex + ethnicity_simplified | tissue_harmonised),
  	formula_variability = ~ age_days*sex + disease + ethnicity_simplified,
    sample_, is_immune,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(10, 100)),
    #output_directory_fit_draws = "/vast/scratch/users/mangiola.s", 
    max_sampling_iterations = 5000
  )

differential_composition_ethnicity_absolute |> saveRDS(output_file_1)

differential_composition_ethnicity_absolute_no_outlier = 
  differential_composition_ethnicity_absolute |> 
  sccomp_remove_outliers(max_sampling_iterations = 4000) 

differential_composition_ethnicity_absolute_no_outlier |> 
  saveRDS(output_file_2)

# Counts RUV Absolute
differential_composition_ethnicity_absolute_no_outlier |>
  remove_unwanted_variation(~  ethnicity_simplified, ~  ethnicity_simplified) |>
  saveRDS(output_file_3)

