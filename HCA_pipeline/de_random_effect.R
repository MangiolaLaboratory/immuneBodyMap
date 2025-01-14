
library(tidySummarizedExperiment)
library(HPCell)
library(magrittr)
library(tibble)
library(dplyr)
# devtools::load_all("~/PostDoc/tidybulk/")


# # Dispersion 2 days calculation
# job::job({
#   library(tictoc)
#   tic()
#   glmGamPoi_subsample = 
#     pseudobulk_sample |> 
#     glmGamPoi::glm_gp(
#       on_disk = T,
#       subsample = TRUE,
#       design = 
#         pseudobulk_sample |> 
#         tidybulk::resolve_complete_confounders_of_non_interest(dataset_id, assay_groups, tissue_groups, disease_groups) |> 
#         colData() |> 
#         droplevels() |> 
#         model.matrix(~ dataset_id + assay_groups + disease_groups, data = _  ) , 
#       verbose = TRUE,
#       use_assay = "counts"
#     ) 
#   toc()
# })
# 
# # Save
# glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_3.rds")$overdispersions
# glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])


result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3"
library(targets)

tar_script({
  
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(qs)
  library(crew)
  library(crew.cluster)
  
  tar_option_set(
    
    
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    
    #cue = tar_cue(mode = "never"), 
    
    workspace_on_error = TRUE,
    format = "qs",
    
    debug = "estimates_chunk",
    
    controller = crew_controller_group(
      
      crew_controller_slurm(
        name = "elastic",
        workers = 500,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_error = 7,
        options_cluster = crew_options_slurm(
          memory_gigabytes_required = c(5, 10, 20, 40, 80, 160), 
          cpus_per_task = 2, 
          time_minutes = c(60*4, 60*4, 60*4, 60*4, 60*24, 60*24),
          verbose = T
        )
      ),
      crew_controller_slurm(
        name = "elastic_big",
        workers = 150,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_error = 5,
        options_cluster = crew_options_slurm(
          memory_gigabytes_required = c(80, 160), 
          cpus_per_task = 2, 
          time_minutes = c(60*24, 60*24),
          verbose = T
        )
      )
    )
    
    
  )
  
  
  #-----------------------#
  # Functions
  #-----------------------#  
  #' Remove Unwanted Effects from a brmsfit Model
  #'
  #' This function utilises fitted values, residuals, and partial predictions from a \code{brmsfit} model
  #' to produce adjusted outcomes that highlight the contribution of a specified factor (e.g. ethnicity) 
  #' while removing unwanted effects from other parts of the model.
  #'
  #' @param fit A \code{brmsfit} object, resulting from a model fitted by \code{\link[brms]{brm}}.
  #' @param factor A character string specifying the response variable (factor) of interest for which 
  #'   unwanted effects should be removed.
  #' @param robust A logical value indicating whether to use robust (median-based) summaries rather than means. 
  #'   Defaults to \code{FALSE}.
  #'
  #' @return A \code{tibble} containing posterior summaries of:
  #'   \itemize{
  #'     \item Adjusted outcomes (prefix: \code{adjusted___}): The adjusted counts that primarily reflect the variation due to the specified factor.
  #'     \item Residuals (prefix: \code{residuals___}): The residual terms used in the adjustment.
  #'     \item Fitted values for the factor (prefix: \code{fitted___}): The model's fitted values for the specified factor with other effects removed.
  #'   }
  #'
  #' @details
  #' This function:
  #' \enumerate{
  #'   \item Extracts the full fitted values from the model.
  #'   \item Calculates residuals by comparing the fitted values to the observed counts.
  #'   \item Normalises these residuals by the offset to ensure they are on the appropriate scale.
  #'   \item Extracts fitted values for the specified factor alone (using \code{re_formula = ~0} to remove random effects).
  #'   \item Combines the fitted factor-specific predictions with the scaled residuals to obtain adjusted outcomes that primarily reflect the factor's contribution.
  #'   \item Summarises these draws (fitted factor-specific, residuals, and adjusted outcomes) into a \code{tibble}.
  #' }
  #'
  #' This approach is useful for examining how a particular factor influences the outcome after "removing" or controlling 
  #' for the other effects included in the model, including random effects and other fixed effects.
  #'
  #' @examples
  #' \dontrun{
  #' # Suppose 'fit' is a brmsfit model object with a response count variable and a factor 'ethnicity'
  #' # We remove unwanted effects to isolate the contribution of 'ethnicity'
  #' adjusted_results <- remove_unwanted_effect(fit, factor = "ethnicity_groups")
  #' }
  #'
  #' @importFrom magrittr %>%
  #' @importFrom dplyr bind_cols
  #' @importFrom tibble as_tibble
  #' @importFrom brms posterior_summary fitted
  #'
  #' @export
  remove_unwanted_effect = function(fit, factor, robust = FALSE, correct_by_offset = T){
    
    # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
    # This places residuals on a consistent scale, making them addable to adjusted predictions later.
    fitted_residuals =   fit |> residuals(robust = robust, summary = FALSE) 
    
    # Correct by offset
    if(correct_by_offset)
      fitted_residuals = fitted_residuals |>
        sweep(2, fit$data$offset |> exp(), FUN = "/")
    
    # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
    # 'resp = factor' focuses on the selected response variable (factor)
    fitted_values_ethnicity <- fitted(fit, re_formula = ~0, resp = factor, summary = FALSE, offset=0)
    
    # Adjusted counts are obtained by adding the factor-specific fitted values and the normalised residuals
    adjusted_counts = fitted_values_ethnicity + fitted_residuals
    
    # Summarise residuals into a tibble, prefixed to denote their source
    fitted_residuals_tbl = 
      fitted_residuals |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    fitted_residuals_tbl |> colnames() = paste0("residuals___", fitted_residuals_tbl |> colnames())
    
    # Summarise the factor-only fitted values into a tibble, prefixed accordingly
    fitted_values_ethnicity_tbl = 
      fitted_values_ethnicity |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    fitted_values_ethnicity_tbl |> colnames() = paste0("fitted___", fitted_values_ethnicity_tbl |> colnames())
    
    # Summarise the adjusted counts (factor + residuals) into a tibble, prefixed for clarity
    adjusted_counts_tbl = 
      adjusted_counts |> 
      posterior_summary(robust = robust) |> 
      as_tibble()
    adjusted_counts_tbl |> colnames() = paste0("adjusted___", adjusted_counts_tbl |> colnames())
    
    # Combine all three resulting tables into one tibble
    adjusted_counts_tbl |> 
      bind_cols(fitted_residuals_tbl) |> 
      bind_cols(fitted_values_ethnicity_tbl)
  }
  
  
  #-----------------------#
  # Pipeline
  #-----------------------#
  
  list(
    
    
    # tar_target(
    #   result_directory,
    #   "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/de_ethnicity_pseudobulk_sample"
    # ),
    # tar_target(
    #   glmGamPoi_overdispersions,
    #   {
    #     glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_3.rds")$overdispersions
    #     glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])
    #     glmGamPoi_overdispersions
    #   }, 
    #   deployment = "main"
    #   
    # ),
    
    tar_target(
      pseudobulk_sample,
      {
        se = 
          loadHDF5SummarizedExperiment("/vast/projects/cellxgene_curated/cellNexus/pseudobulk_sample_is_immune") |> 
          filter(is_gene_shared) |> 
          
          #---------------------------------#
          # Edit or add more filters here for analyses
          #---------------------------------#
          filter(is_immune & do_analyse) 
        
        
        
        samples_with_right_number_of_detected_genes = 
          (se |> assay() > 0) |> 
          colSums() |> 
          divide_by(nrow(se)) |> 
          between(0.3, 0.7)
        
        se = se[,samples_with_right_number_of_detected_genes] 
        
        se = 
          se |> 
          keep_abundant(design = 
                          se |> 
                          
                          # Discretise the age for the following operation
                          mutate(is_old_individual = age_days > 50*365) |> 
                          resolve_complete_confounders_of_non_interest(tissue_groups, sex, ethnicity_groups, is_old_individual) |> 
                          colData() |> 
                          droplevels() |> 
                          model.matrix(~ tissue_groups + sex + ethnicity_groups + is_old_individual, data = _  ), 
                        minimum_counts = 100
          ) |> 
          
          # Get scaling factor
          scale_abundance(method = "TMMwsp", reference_sample = "0edf00b9d5cd39b046f90be198fb07db___1") |> 
          
          # Drop sex unknown as causes problem during fit
          filter(sex != "unknown") |> 
          
          # Eliminate complete confounders
          tidybulk:::resolve_complete_confounders_of_non_interest(assay_groups, dataset_id, disease_groups) |> 
          
          # sibrary size factor is the reciproque of the multiplier (correction factor)
          mutate(offset = log(1/multiplier))
        
        # # Add dispersion
        # rowData(se)  = 
        #   rowData(se) |> 
        #   as_tibble(rownames = ".feature") |> 
        #   left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |> 
        #   data.frame(row.names = ".feature") |> DataFrame()
        
        se
        
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big")),
      memory = "persistent"
    ),
    
    # Split in gene chunks
    tar_target(
      feature_df, 
      pseudobulk_sample |> 
        distinct(.feature)|> 
        group_by(.feature) |> 
        tar_group(), 
      iteration = "group",
      packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    tar_target(
      se_df, 
      feature_df |> mutate(se = map(.feature, ~ 
                                      pseudobulk_sample[.x, , drop=FALSE]
      ))  , 
      pattern = map(feature_df),
      packages = c( "brms", "glue"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    tar_target(
      estimates_chunk, 
      
      se_df |> mutate(brms_fit = map(se, ~ {
        
        data = 
          .x |>
          as_tibble() |> 
          mutate(counts = counts |> as.integer()) |> 
          droplevels()
        
        # Drop NA counts. Not sure why they are there. E.g.:
        # $ .feature             <chr> "ENSG00000134419"
        # $ .sample              <chr> "3bfa31867cc1c823e0cb2f1ff24df26b___1"
        # $ counts               <int> NA
        # $ gene_presence        <int> 25
        # $ counts_scaled        <dbl> 38316.51
        # $ sample_id            <chr> "3bfa31867cc1c823e0cb2f1ff24df26b"
        # $ is_immune            <int> 1
        # $ do_analyse           <lgl> TRUE
        # $ donor_id             <chr> "one_Ten"
        # $ title                <chr> "Individual Single-Cell RNA-seq PBMC Data from Schulte-Schrepping et al."
        # $ dataset_id           <chr> "5e717147-0f75-4de1-8bd2-6fda01b8d75f"
        # $ collection_id        <chr> "b9fc3d70-5a72-4479-a046-c2cc1ab19efc"
        # $ age_days             <int> 10585
        # $ sex                  <chr> "male"
        # $ ethnicity_groups     <fct> European
        # $ tissue_groups        <chr> "blood"
        # $ assay_groups         <fct> 10x Genomics 3
        # $ disease_groups       <fct> Normal
        # $ age_bin_sex_specific <fct> Young Adulthood
        # $ TMM                  <dbl> 2.538842
        # $ multiplier           <dbl> 1.394764e-05
        # $ offset               <dbl> -11.1802
        # $ is_gene_shared       <lgl> TRUE
        # $ .abundant            <lgl> TRUE
        # $ dispersion           <dbl> 0.4960158
        n_NAs = data |> filter(counts |> is.na()) |> nrow()
        if(n_NAs > 0){
          warning(glue("You have {n_NAs} NAs in counts. They have been filtered out"))
          data = 
            data |> 
            filter(!counts |> is.na())
        }
        
        # Check if dispersion estimation has failed
        if(data |> pull(dispersion) |> unique() |> is.na()){
          warning("The dispersion calculation has failed. 1 is given as default prior.")
          data = data |> mutate(dispersion = 1)
        }
        
        # Define the model formula
        formula <- bf(
          
          # Formula for counts
          counts ~ 1 + offset(offset) + age_bin_sex_specific*sex + disease_groups + ethnicity_groups + assay_groups + 
            (1 | dataset_id) + 
            (1 + age_bin_sex_specific*sex + ethnicity_groups | tissue_groups),
          
          # Formula for dispersion
          shape ~ 1 + disease_groups + assay_groups + ethnicity_groups + (1 | tissue_groups)  # Model 'shape' as a function of scaled 'disp'
          
          # Using the externally, eBayes inferred overdispersion
          # shape ~ 1 + offset(log(1/dispersion))
         )
        
        prior = c(
          prior(normal(i, 5), class = Intercept),
          prior(normal(0, 2), class = Intercept, dpar = shape),
          prior(normal(0, 5), class = b),
          prior(normal(0, 2), class = b, dpar = shape)
        ) |> 
          substitute(env = list(i = mean(log1p(data$counts / exp(data$offset))))) |> 
          eval()
        
        chains = 2
        inits <- list(Intercept = mean(log1p(data$counts / exp(data$offset))))
        inits <- replicate(chains, inits, simplify = FALSE)
        
        
        brm(
          formula = formula,
          data = data,
          family = zero_inflated_negbinomial(),
          prior = prior,
          chains = chains,
          cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)), #, threads = 2,
          warmup = 300, 
          refresh = 10,
          backend = "cmdstanr", 
          #sparse = TRUE,
          #save_model = glue("{external_directory}~/temp.rds"),
          #algorithm = "pathfinder",
          init = inits,
          iter = 400  # Increase iterations for better convergence
        )
        
      })) |> 
        
        # Drop data because it is withn the brms object
        select(-se), 
      pattern = map(se_df),
      packages = c( "brms", "glue", "dplyr", "purrr", "SummarizedExperiment", "tidySummarizedExperiment"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      cue = tar_cue(mode = "never")
      
    ),
    tar_target(
      summary, 
      estimates_chunk |> 
        mutate(summary = map(brms_fit, ~ .x |> hypothesis(
          c(
            "Europeans" = "(ethnicity_groupsAfrican
    + ethnicity_groupsEastAsian
    + ethnicity_groupsHispanicDLatinAmerican
    + ethnicity_groupsSouthAsian
    + `ethnicity_groupsNativeAmerican&PacificIslander`) / 5 = 0",
            "EastAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsNativeAmerican&PacificIslander`
     - 5 * ethnicity_groupsEastAsian
     ) / 5 = 0",
            "SouthAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsEastAsian
     + `ethnicity_groupsNativeAmerican&PacificIslander`
     - 5 * ethnicity_groupsSouthAsian 
     ) / 5 = 0",
            "African" = "(
       ethnicity_groupsEastAsian
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsNativeAmerican&PacificIslander`
     - 5 * ethnicity_groupsAfrican 
     ) / 5 = 0",
            "HispanicDLatinAmerican" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsEastAsian
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsNativeAmerican&PacificIslander`
     - 5 * ethnicity_groupsHispanicDLatinAmerican 
     ) / 5 = 0",
            "NativeAmericanPacificIslander" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + ethnicity_groupsEastAsian
     - 5 * `ethnicity_groupsNativeAmerican&PacificIslander` 
     ) / 5 = 0"
      ),
      
      # Median instead and mad of mean and sd
      robust=TRUE)
        )) |> 
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    tar_target(
      effect_removed, 
      estimates_chunk |> 
        mutate(brms_fit_adjusted = map(brms_fit, ~ .x |> remove_unwanted_effect("ethnicity_groups", robust = TRUE) )) |> 
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    )
    
  )
  
  
}, ask = FALSE, script = glue::glue("{result_directory}/_targets.R"))


job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue::glue("{result_directory}/_targets.R"),
    store = glue::glue("{result_directory}/_targets")
  )
  
})


pseudobulk_sample = tar_read(pseudobulk_sample, store = glue::glue("{result_directory}/_targets"))

lib_size = pseudobulk_sample |> assay() |> colSums()
plot(log(lib_size), log(colData(pseudobulk_sample)$multiplier))



effect_removed = 
  tar_read(
    effect_removed,
    store = glue::glue("{result_directory}/_targets")
  ) 

effect_removed = 
  effect_removed |>
  filter(map_int(brms_fit_adjusted, nrow) == 4926 ) |> 
  mutate(brms_fit_adjusted = map(brms_fit_adjusted, ~ .x |> 
                                   select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                                   mutate(sample_id = colnames(pseudobulk_sample)), 
                                 .progress = TRUE
                                )) |>
  unnest(brms_fit_adjusted)

m = 
  effect_removed |> 
  pivot_wider(names_from = sample_id, values_from = adjusted___Estimate) |> 
  select(-tar_group) |> 
  tidybulk:::as_matrix(rownames = ".feature") |> 
  as("sparseMatrix")  |> 
  Matrix::Matrix(sparse = T)

pseudobulk_sample = pseudobulk_sample[rownames(m),, drop=FALSE ] 
assay(pseudobulk_sample, "counts_adjusted_ethnicity") = m
pseudobulk_sample = pseudobulk_sample[((assay(pseudobulk_sample, "gene_presence") > 0) |> rowSums() > (ncol(pseudobulk_sample) * 0.8)),,drop=FALSE ]

# Cap infinite
max_rm_infinite = 
  pseudobulk_sample |> 
  assay("counts_adjusted_ethnicity") |> 
  _[!pseudobulk_sample |> assay("counts_adjusted_ethnicity") |> is.infinite()] |> 
  quantile(0.999)
  
pseudobulk_sample |> 
  assay("counts_adjusted_ethnicity") |> 
  _[pseudobulk_sample |> assay("counts_adjusted_ethnicity") > max_rm_infinite] = 
  max_rm_infinite

pseudobulk_sample |> 
  assay("counts_adjusted_ethnicity") |> 
  _[pseudobulk_sample |> assay("counts_adjusted_ethnicity") < 0] = 
  0

summaries = 
  tar_read(
    summary,
    store = glue::glue("{result_directory}/_targets")
  ) |>
  mutate(summary = map(summary, ~ .x %$% hypothesis |> as_tibble())) |>
  unnest(summary) |>
  filter(Star == "*") |>
  filter(.feature %in% rownames(pseudobulk_sample)) |> 
  mutate(closest_to_zero = pmin(abs(CI.Lower), abs(CI.Upper))) |>
  add_count(.feature) |> 
  filter(n < 4) |> 
  with_groups(Hypothesis, ~ .x |> arrange(desc(closest_to_zero)) |> dplyr::slice(1:50))

pseudobulk_sample_for_PCA = pseudobulk_sample

pseudobulk_sample_for_PCA |> 
  as("SingleCellExperiment") |> 
  zellkonverter::writeH5AD(file = "~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity.h5ad", compression = "gzip")
system("~/bin/rclone copy ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity.h5ad box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/")


pseudobulk_sample_for_PCA = 
  pseudobulk_sample_for_PCA |> 
  select(-contains("PC"), -contains("tSNE"), -contains("UMAP")) |> 
  # _[summaries |> pull(.feature) |> unique(), , drop=FALSE] |> 
  filter(ethnicity_groups != "Other/Unknown") |> 
  tidybulk::reduce_dimensions(method = "PCA", .abundance = counts_adjusted_ethnicity, .dims = 20 ) |> 
  tidybulk::reduce_dimensions(method = "tSNE", .abundance = counts_adjusted_ethnicity, initial_dims = 10, .dims = 3)  |> 
  tidybulk::reduce_dimensions(method = "UMAP", .abundance = counts_adjusted_ethnicity, pca = 10, .dims = 3, calculate_for_pca_dimensions = NULL)  

pseudobulk_sample_for_PCA  = pseudobulk_sample_for_PCA |> filter(PC1 < 60)
pseudobulk_sample_for_PCA  = pseudobulk_sample_for_PCA |> filter(PC4 > -20)



colData(pseudobulk_sample)$sum = pseudobulk_sample |> assay("counts_adjusted_ethnicity") |> colSums()

pseudobulk_sample_for_PCA |> 
  pivot_sample() |> 
  ggplot(aes(tSNE1, tSNE2, color = ethnicity_groups)) +
  geom_point() +
  guides(color = "none")


pseudobulk_sample_for_PCA |> 
  pivot_sample() |> 
  select(contains("PCA"), everything()) |> 
GGally::ggpairs(columns = 19:38, ggplot2::aes(colour=`ethnicity_groups`))

pseudobulk_sample_for_PCA |>
  pivot_sample() |> 
  plot_ly(
    x = ~`tSNE1`,
    y = ~`tSNE2`,
    z = ~`tSNE3`,
    color = ~disease_groups
  ) %>%
  add_markers(size = I(10))


# tSNE coloured by PCA
pseudobulk_sample_for_PCA |> 
  pivot_sample() |> 
  select(sample_id, contains("PC"), contains("tSNE")) |> 
  pivot_longer(contains("PC"), names_to = "PC_number", values_to = "PC_value") |> 
  filter(PC_number %in% paste0("PC", 1:10)) |> 
  mutate(
    PC_value = if_else(PC_value > quantile(PC_value, 0.9), quantile(PC_value, 0.9),PC_value ),
    PC_value = if_else(PC_value < quantile(PC_value, 0.1), quantile(PC_value, 0.1),PC_value )
  ) |> 
  ggplot(aes(tSNE1, tSNE2, color = PC_value)) +
  geom_point(size = 0.3) +
  facet_wrap(~ PC_number) +
  scale_color_distiller(palette = "Spectral") +
  guides(color = "none") +
  theme_bw()

x |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_chunk_062dd2621e3cb7e1.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_chunk_062dd2621e3cb7e1.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/")


tar_workspace(
  summary_6e9d63f17424cc2f, 
  store = glue::glue("{result_directory}/_targets"),
  script = glue::glue("{result_directory}/_targets.R")
)



tar_meta(store = glue::glue("{result_directory}/_targets")) |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  select(name, error)

library(tidybayes)
library(brms)
library(magrittr)
library(ggallin)




fit = tar_read_raw("estimates_chunk_34ae3b7786cdde9a", store = glue::glue("{result_directory}/_targets"), branches = 1) |> 
  pull(brms_fit) |> 
  _[[1]] 

fit |> fitted() |> head()

fit |> 
  
  # Remove unwanted effects
  remove_unwanted_effect("ethnicity_groups", robust = TRUE) |> 
  mutate(o = fit$data$offset ) |> 
 # filter(residuals___Estimate<1e20) |> 
  
  # Plot
  ggplot(aes(o,adjusted___Estimate)) + 
 geom_errorbar(aes(ymin = residuals___Q2.5, ymax = residuals___Q97.5), alpha = 0.01) +
  geom_point(shape = ".") + 
  scale_y_continuous(trans = pseudolog10_trans) +
  #scale_y_log10() + scale_x_log10() + 
  geom_hline(yintercept = 0, colour = "red") 

pp_check(fit) + scale_x_continuous(trans = pseudolog10_trans) 


fit2 = fit |> update(
  threads = 4, 
  cores = 2,
  prior = c(
    prior(normal(6, 5), class = Intercept),
    prior(normal(0, 2), class = Intercept, dpar = shape),
    prior(normal(0, 5), class = b),
    prior(normal(0, 2), class = b, dpar = shape)
  ) 
 
)

check_brms <- function(model,             # brms model
                       integer = FALSE,   # integer response? (TRUE/FALSE)
                       plot = TRUE,       # make plot?
                       ...                # further arguments for DHARMa::plotResiduals 
) {
  
  mdata <- brms::standata(model)
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model)),
    observedResponse = model$data$counts, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, re.form = NA)),
      1,
      mean),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
  
}

check_brms(fit, integer = T)

# 1 c NA   
# 2 estimates_chunk_77429ffc4ef9e584 NA   
# 3 estimates_chunk_7cbe0084f123ee96 NA   
# 4 estimates_chunk_b1bd8ace2583e01e NA   
# 5 estimates_chunk_5ef9151358ee4d99 NA   
# 6 estimates_chunk_ec9ec278f52ec218 NA   
# 7 estimates_chunk_c0bf30ef83e8e21d NA   
# 8 estimates_chunk_34ae3b7786cdde9a NA   
# 9 estimates_chunk_e150400caf575755 NA   
# 10 estimates_chunk_f0801034799af6ad NA 