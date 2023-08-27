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
output_blood = args[[4]]
output_file_2 = args[[5]]


my_data =
  readRDS(input_file) |>

  # Scale days
  mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>

	unite("group", c(tissue_harmonised , file_id), remove = FALSE) |>
	
  unite("tissue_harmonised_ethnicity", c(tissue_harmonised , ethnicity_simplified), remove = FALSE) 

if(filter_blood=="TRUE"){


  res_relative_blood =
    my_data |>

    # nly blood
    filter(tissue_harmonised=="blood") |>

    # Estimate
  	sccomp_estimate(
    	formula_composition = ~ age_days*sex + disease + ethnicity_simplified + assay_simplified + disease + group + (1 + age_days + sex + ethnicity_simplified | tissue_harmonised),
    	formula_variability = ~ age_days*sex + disease,
      sample_, cell_type_harmonised, counts_from_tissue,
      check_outliers = F,
      approximate_posterior_inference = FALSE,
      cores = 10,
      mcmc_seed = 42,
      verbose = T,
      prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40)),
      #output_directory_fit_draws = "/vast/scratch/users/mangiola.s"
    	max_sampling_iterations = 5000
    )

  res_relative_blood |> saveRDS(output_blood)

  predicted_blood_composition =
    res_relative_blood |>
    sccomp_predict(
      formula_composition = ~ 0 + ethnicity  + sex  + age_days +  assay ,
      number_of_draws = 500,
      new_data =
        my_data |>
        distinct(sample_, sex , ethnicity , age_days, assay )
    )

  predicted_proportion_of_naive =
    predicted_blood_composition |>
    mutate(is_naive = cell_type_harmonised |> str_detect("naive|stem")) |>
    with_groups(c(sample_, is_naive), ~ .x |> summarise(predicted_proportion = sum(proportion_mean))) |>
    filter(is_naive)

  counts_to_subtract =
    my_data |>
    mutate(is_naive = cell_type_harmonised |> str_detect("naive|stem")) |>
    with_groups(sample_, ~ .x |> mutate(exposure = sum(counts_from_tissue))) |>
    with_groups(c(sample_, exposure, is_naive), ~ .x |> summarise(n = sum(counts_from_tissue))) |>
    complete(nesting(sample_, exposure), is_naive, fill = list(n = 0)) |>
    with_groups(c(sample_), ~ .x |> mutate(observed_proportion = n/sum(n))) |>
    filter(is_naive) |>
    left_join(predicted_proportion_of_naive) |>
    mutate(blood_contamination = observed_proportion / predicted_proportion) |>
    mutate(total_blood_count = exposure * blood_contamination) |>
    left_join(predicted_blood_composition) |>
    mutate(counts_from_blood = floor(proportion_mean * total_blood_count)) |>
    select(sample_, cell_type_harmonised, counts_from_blood)

  my_data =
    my_data |>
    left_join(counts_to_subtract) |>
    mutate(counts_blood_free = if_else(
      !tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone", "thymus"),
      counts_from_tissue - counts_from_blood,
      counts_from_tissue
    )) |>
    mutate(counts_blood_free = pmax(counts_blood_free, 0))


}


res_relative =
  my_data |>

  # Estimate
  sccomp_estimate(
  	formula_composition = ~ disease + age_days*sex + ethnicity_simplified + assay_simplified + group + (1 + age_days*sex + ethnicity_simplified | tissue_harmonised),
  	formula_variability = ~ disease,
    sample_, cell_type_harmonised,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40)),
    #output_directory_fit_draws = "/vast/scratch/users/mangiola.s", 
    max_sampling_iterations = 5000
  )


res_relative |>
  saveRDS(output_file_1)

# Remove unwanted variation
res_relative |>
  remove_unwanted_variation(~  ethnicity_simplified, ~  ethnicity_simplified) |>
  saveRDS(output_file_2)
