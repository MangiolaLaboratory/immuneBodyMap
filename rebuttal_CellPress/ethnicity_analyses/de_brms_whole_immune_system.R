
library(tidySummarizedExperiment)
library(HPCell)
library(magrittr)
library(tidyverse)
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
# glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_6.rds")$overdispersions
# glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])


result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6"
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
  #' This function calculates posterior residuals from a \code{brmsfit} model and combines them with 
  #' factor-specific fitted values (potentially excluding random effects or other parts of the model), 
  #' thereby producing adjusted outcomes that highlight the contribution of a specified factor or subset 
  #' of model terms.
  #'
  #' @param fit A \code{brmsfit} object, resulting from a model fitted by \code{\link[brms]{brm}}.
  #' @param newdata A data frame or list containing new data. Passed to \code{\link[brms]{fitted}} 
  #'   to obtain factor-specific fitted values at specified covariate levels.
  #' @param robust A logical value indicating whether to use robust (median-based) summaries rather 
  #'   than means. Defaults to \code{FALSE}.
  #' @param correct_by_offset A logical value indicating whether to divide the residuals by 
  #'   \code{exp(offset)} (from \code{fit$data$offset}). Defaults to \code{TRUE}.
  #' @param re_formula A formula specifying which random effects (if any) to include when generating 
  #'   fitted values. Defaults to \code{~0}, which removes random effects and thus isolates the 
  #'   contribution of fixed effects in the new data.
  #'
  #' @return A \code{tibble} containing posterior summaries of:
  #'   \itemize{
  #'     \item Adjusted outcomes (prefix: \code{adjusted___}): The combined values of the specified 
  #'     factor's fitted counts and the residuals.
  #'     \item Residuals (prefix: \code{residuals___}): The model's posterior residuals, possibly 
  #'     normalised by the offset.
  #'     \item Fitted values for the factor (prefix: \code{fitted___}): The model's fitted values based 
  #'     on the \code{re_formula} and provided \code{newdata}.
  #'   }
  #'
  #' @details
  #' The function proceeds as follows:
  #' \enumerate{
  #'   \item Extracts posterior residuals (via \code{\link[brms]{residuals}}).
  #'   \item (Optionally) divides these residuals by the exponential of the offset, if \code{correct_by_offset = TRUE}.
  #'   \item Obtains new fitted values from the model (via \code{\link[brms]{fitted}}), usually excluding random effects 
  #'         by specifying \code{re_formula = ~0}.
  #'   \item Adds these residuals to the factor-specific fitted values to obtain adjusted outcomes 
  #'         that highlight the contribution of the factor of interest.
  #'   \item Summarises all these draws (residuals, fitted values, adjusted outcomes) and returns them 
  #'         in a single \code{tibble}.
  #' }
  #'
  #' This method is particularly useful for examining how a factor or other subset of the model 
  #' affects the outcome when other model components (e.g., random intercepts) are removed. 
  #' It can assist in visualising or quantifying the partial contribution of certain terms.
  #'
  #' @examples
  #' \dontrun{
  #' # Suppose 'fit' is a brmsfit model object predicting a 'counts' outcome
  #' # We create a new data frame 'some_data' for which we want partial predictions
  #' adjusted_results <- remove_unwanted_effect(
  #'   fit,
  #'   newdata = some_data,
  #'   robust = TRUE,
  #'   correct_by_offset = TRUE,
  #'   re_formula = ~0
  #' )
  #' }
  #'
  #' @importFrom magrittr %>%
  #' @importFrom dplyr bind_cols
  #' @importFrom tibble as_tibble
  #' @importFrom brms posterior_summary fitted residuals
  #'
  #' @export
  remove_unwanted_effect = function(fit, newdata, robust = FALSE, correct_by_offset = T, re_formula = ~0){
    
    # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
    # This places residuals on a consistent scale, making them addable to adjusted predictions later.
    fitted_residuals =   fit |> residuals(robust = robust, summary = FALSE) 
    
    # Correct by offset
    if(correct_by_offset)
      fitted_residuals = fitted_residuals |>
        sweep(2, fit$data$offset |> exp(), FUN = "/")
    
    # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
    # 'resp = factor' focuses on the selected response variable (factor)
    fitted_values_ethnicity <- fitted(fit, newdata = newdata, re_formula = re_formula, summary = FALSE, offset=0)
    
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
  
  remove_unwanted_effect_new = function(fit, newdata, robust = FALSE, correct_by_offset = T, re_formula = ~0){
    
    # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
    # This places residuals on a consistent scale, making them addable to adjusted predictions later.
    fitted_residuals =   fit |> predictive_error(robust = robust, summary = FALSE, offset = 0) 
    
    # Correct by offset
    if(correct_by_offset)
      fitted_residuals = fitted_residuals |>
        sweep(2, fit$data$offset |> exp(), FUN = "/")
    
    # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
    # 'resp = factor' focuses on the selected response variable (factor)
    fitted_values_ethnicity <- posterior_epred(fit, newdata = newdata, re_formula = re_formula,  offset=0)
    
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
    #     glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_6.rds")$overdispersions
    #     glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])
    #     glmGamPoi_overdispersions
    #   }, 
    #   deployment = "main"
    #   
    # ),
    
    # This target loads and processes the pseudobulk sample data. It imports a HDF5 SummarizedExperiment, 
    # applies filters to retain shared genes, immune cells, and samples marked for analysis, integrates age metadata,
    # filters for common genes and samples with an appropriate number of detected genes, computes the mean library size, 
    # selects a reference sample, and performs normalisation and scaling.
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
        
        # TEMPORARY BECAUSE I FORGOT TO INTEGRATE AGE BINS
        se = se |> 
          left_join(
            readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_6/cell_metadata_1_0_6_sccomp_input_counts.rds") |> 
              distinct(sample_id,  age_days, age_bin) 
          )
        
        # Filter common genes
        se = se[((assay(se, "gene_presence") > 0) |> rowSums() > (ncol(se) * 0.95)),,drop=FALSE ]
        
        # Filter samples that have enough genes > 0 but not too many
        samples_with_right_number_of_detected_genes = 
          (se |> assay() > 0) |> 
          colSums() |> 
          divide_by(nrow(se)) |> 
          dplyr::between(0.3, 1)
        
        se = se[,samples_with_right_number_of_detected_genes] 
        
        # Compute mean library size
        mean_library_size <- se |>
          assay("counts") |>
          _[nrow(se) |> seq_len() |> sample(size = 2000), ] |> 
          colSums() |>
          mean()
        
        # Optional: retrieve the sample name (column name in the SummarizedExperiment)
        reference_sample <- colnames(se)[
          se |>
            assay("counts") |>
            colSums() |>
            {\(x) abs(x - mean_library_size)}() |>  # Calculate absolute difference from the mean
            which.min()                             # Identify the smallest difference
        ]
        
        
        se = 
          se |> 
          keep_abundant(design = 
                          se |> 
                          
                          # Discretise the age for the following operation
                          mutate(is_old_individual = age_days > 50*365) |> 
                          
                          # This is to resolve some confounders to preserve the genes.
                          # In this case we care about data variability, not the actual meaning of the variables
                          resolve_complete_confounders_of_non_interest(tissue_groups, sex, ethnicity_groups, is_old_individual) |> 
                          colData() |> 
                          droplevels() |> 
                          model.matrix(~ tissue_groups + sex___altered + ethnicity_groups___altered + is_old_individual___altered, data = _  ), 
                        minimum_counts = 100
          ) |> 
          
          # Get scaling factor
          scale_abundance(method = "TMMwsp", reference_sample = reference_sample) |> 
          
          # Drop sex unknown as causes problem during fit
          mutate(
            sex = if_else(sex |> is.na(), "unknown", sex),
            ethnicity_groups = if_else(ethnicity_groups |> is.na(), "Other/Unknown", ethnicity_groups)
          ) |> 
          filter(sex != "unknown") |> 
          filter(!age_bin |> is.na()) |> 
          
          # Eliminate complete confounders
          tidybulk:::resolve_complete_confounders_of_non_interest(assay_groups, dataset_id, disease_groups) |> 
          
          # sibrary size factor is the reciproque of the multiplier (correction factor)
          mutate(offset = log(1/multiplier)) |> 
          
          # Set intercept
          mutate(
            ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
            assay_groups___altered = fct_relevel(assay_groups___altered, "10x Genomics 3"),
            disease_groups___altered = fct_relevel(disease_groups___altered, "Normal"),
            age_bin = fct_relevel(age_bin, "Adolescence")
          ) 
        
        # # Add dispersion
        # rowData(se)  = 
        #   rowData(se) |> 
        #   as_tibble(rownames = ".feature") |> 
        #   left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |> 
        #   data.frame(row.names = ".feature") |> DataFrame()
        
        se
        
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble", "forcats"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big")),
      memory = "persistent", 
      error = "stop"
    ),
    
    # This target extracts unique features from the pseudobulk sample and groups them into 
    # chunks for parallel processing.
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
    
    # This target creates a list-column of SummarizedExperiment objects,
    # with each object corresponding to a distinct feature.
    tar_target(
      se_df, 
      feature_df |> mutate(se = map(.feature, ~ 
                                      pseudobulk_sample[.x, , drop=FALSE]
      ))  , 
      pattern = map(feature_df),
      packages = c( "brms", "glue"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # This target fits Bayesian models on chunks of the data. It processes each feature's data, handles missing values,
    # defines the model specification with priors, and runs the Bayesian inference using the brm function.
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
        # $ age_bin <fct> Young Adulthood
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
            filter(!counts |> is.na()) |> 
            droplevels()
        }
        
        # # Check if dispersion estimation has failed
        # if(data |> pull(dispersion) |> unique() |> is.na()){
        #   warning("The dispersion calculation has failed. 1 is given as default prior.")
        #   data = data |> mutate(dispersion = 1)
        # }
        
        # Define the model formula
        formula <- bf(
          
          # Formula for counts
          counts ~ 1 + offset(offset) + age_bin*sex + disease_groups___altered + ethnicity_groups + assay_groups___altered + 
            (1 | dataset_id___altered) + 
            (1 + age_bin*sex + ethnicity_groups | tissue_groups),
          
          # Formula for dispersion
          shape ~ 1 + disease_groups___altered + assay_groups___altered + ethnicity_groups + (1 | tissue_groups)  # Model 'shape' as a function of scaled 'disp'
          
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
     
    # This target summarises the fitted Bayesian models by performing hypothesis tests for ethnicity contrasts 
    # and extracting convergence diagnostics (Rhat) for the ethnicity parameters.
    tar_target(
      summary,
      estimates_chunk |>
        mutate(summary = map(brms_fit, ~ .x |> hypothesis(
          c(
            "Europeans" = "(ethnicity_groupsAfrican
    + ethnicity_groupsEastAsian
    + ethnicity_groupsHispanicDLatinAmerican
    + ethnicity_groupsSouthAsian
    + `ethnicity_groupsJapanese`) / 5 = 0",
            "EastAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsEastAsian
     ) / 5 = 0",
            "SouthAsian" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsEastAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsSouthAsian
     ) / 5 = 0",
            "African" = "(
       ethnicity_groupsEastAsian
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsAfrican
     ) / 5 = 0",
            "HispanicDLatinAmerican" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsEastAsian
     + ethnicity_groupsSouthAsian
     + `ethnicity_groupsJapanese`
     - 5 * ethnicity_groupsHispanicDLatinAmerican
     ) / 5 = 0",

      "Japanese" = "(
       ethnicity_groupsAfrican
     + ethnicity_groupsHispanicDLatinAmerican
     + ethnicity_groupsSouthAsian
     + ethnicity_groupsEastAsian
     - 5 * `ethnicity_groupsJapanese`
     ) / 5 = 0"
      ),
   #        c(
   #          "African" = "(ethnicity_groupsEuropean
   #  + ethnicity_groupsEastAsian
   #  + ethnicity_groupsHispanicDLatinAmerican
   #  + ethnicity_groupsSouthAsian
   #  + `ethnicity_groupsJapanese`) / 5 = 0",
   #          
   #          "Europeans" = "(
   #   ethnicity_groupsEastAsian
   # + ethnicity_groupsHispanicDLatinAmerican
   # + ethnicity_groupsSouthAsian
   # + `ethnicity_groupsJapanese`
   # - 5 * ethnicity_groupsEuropean
   # ) / 5 = 0",
   #          
   #          "EastAsian" = "(
   #   ethnicity_groupsEuropean
   # + ethnicity_groupsHispanicDLatinAmerican
   # + ethnicity_groupsSouthAsian
   # + `ethnicity_groupsJapanese`
   # - 5 * ethnicity_groupsEastAsian
   # ) / 5 = 0",
   #          
   #          "SouthAsian" = "(
   #   ethnicity_groupsEuropean
   # + ethnicity_groupsHispanicDLatinAmerican
   # + ethnicity_groupsEastAsian
   # + `ethnicity_groupsJapanese`
   # - 5 * ethnicity_groupsSouthAsian
   # ) / 5 = 0",
   #          
   #          "HispanicDLatinAmerican" = "(
   #   ethnicity_groupsEuropean
   # + ethnicity_groupsEastAsian
   # + ethnicity_groupsSouthAsian
   # + `ethnicity_groupsJapanese`
   # - 5 * ethnicity_groupsHispanicDLatinAmerican
   # ) / 5 = 0",
   #          
   #          "Japanese" = "(
   #   ethnicity_groupsEuropean
   # + ethnicity_groupsHispanicDLatinAmerican
   # + ethnicity_groupsSouthAsian
   # + ethnicity_groupsEastAsian
   # - 5 * `ethnicity_groupsJapanese`
   # ) / 5 = 0"
   #        ),

      # Median instead and mad of mean and sd
      robust=TRUE,
      alpha = 0.1
      )
        )) |>
        
      mutate(Rhat = map_dbl(brms_fit, 
                            ~ summary(.x)$fixed |> 
                              as_tibble(rownames = "par") |> 
                              filter(par |> str_detect("ethnicity")) |> 
                              pull(Rhat) |>
                              max()
                          )) |> 
        select(-brms_fit),

      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
   # This target generates adjusted model estimates by removing unwanted effects from the fitted Bayesian models,
   # thereby isolating the effects of interest. Here, nuisance covariates are set to NA and removed from the predictions.
    tar_target(
      effect_removed, 
      estimates_chunk |> 
        mutate(brms_fit_adjusted = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups=NA, sex = NA, age_bin = NA, disease_groups = NA, dataset_id = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = TRUE, 
          re_formula = ~ 0
        ))) |> 
        mutate(brms_fit_adjusted_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          newdata = .x$data |> mutate(assay_groups=NA, sex = NA, age_bin = NA, disease_groups = NA, dataset_id = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = FALSE, correct_by_offset = FALSE,
          re_formula = ~ 0
        ))) |> 
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
   # This target produces adjusted estimates from the Bayesian models, removing unwanted effects while retaining 
   # the tissue group random effect, thus preserving variability associated with tissue-specific factors.
    tar_target(
      effect_removed_keep_tissue, 
      estimates_chunk |> 
        mutate(brms_fit_adjusted = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups=NA, sex = NA, age_bin = NA, disease_groups = NA, ethnicity_groups = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = TRUE, 
          re_formula = ~ (1 | tissue_groups)
        ))) |> 
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    )
    
  )
  
  
}, ask = FALSE, script = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets.R"))


job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets.R"),
    store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")
  )
  
})



# Get genes with wrong intercept
meta_to_speed_up = 
  tar_meta(starts_with("estimates_chunk_"), store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")) |> 
  filter(name |> str_detect("keep_tissue", negate = T)) |> 
  filter(
    error |> is.na(), 
    !data |> is.na() 
  )

names_to_drop = 
  
  # FOR INCOMPLETE PIPELINE
  meta_to_speed_up |> 
  pull(name) |> 
  enframe() |> 
  mutate(keep = 
    map_lgl(value, 
      ~ .x |> tar_read_raw(
      meta = meta_to_speed_up, 
      store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")
      ) |> 
        pull(brms_fit) |> 
        _[[1]] |>
        summary() %$%
        fixed |> 
        rownames() |> 
        str_detect("Afric") |> 
        any(), 
      .progress = TRUE
    ))

pseudobulk_sample = tar_read(pseudobulk_sample, store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets"))

lib_size = pseudobulk_sample |> assay() |> colSums()
plot(log(lib_size), log(colData(pseudobulk_sample)$multiplier))

# FOR INCOMPLETE PIPELINE
meta_to_speed_up = 
  tar_meta(starts_with("effect_removed_"), store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")) |> 
  filter(name |> str_detect("keep_tissue", negate = T)) |> 
  filter(
    error |> is.na(), 
    !data |> is.na() 
  )

effect_removed = 
  
  # FOR INCOMPLETE PIPELINE
  meta_to_speed_up |> 
  pull(name) |> 
  map_dfr(
    tar_read_raw,
    meta = meta_to_speed_up, 
    store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets"), 
    .progress = TRUE
  ) |> 

  # FOR COMPLETE PIPELINE
  # tar_read(
  #   effect_removed,
  #   store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")
  # )  |>
  filter(map_int(brms_fit_adjusted, nrow) == 5230 ) |> 
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

# FOR INCOMPLETE PIPELINE
meta_to_speed_up = 
  tar_meta(starts_with("summary_"), store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")) |> 
  filter(name |> str_detect("keep_tissue", negate = T)) |> 
  filter(
    error |> is.na(), 
    !data |> is.na() 
  )

summaries = 
  
  # FOR INCOMPLETE PIPELINE
  meta_to_speed_up |> 
  pull(name) |> 
  map_dfr(
    tar_read_raw,
    meta = meta_to_speed_up, 
    store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets"), 
    .progress = TRUE
  ) |> 
  
  # # For complete pipelines
  # tar_read(
  #   summary,
  #   store =  glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")
  # ) |>
  
  mutate(summary = map(summary, ~ .x %$% hypothesis |> as_tibble())) |>
  unnest(summary) |>
  filter(Rhat |> dplyr::between(0.90, 1.1)) |> 
  filter(Star == "*") |>
  filter(.feature %in% rownames(pseudobulk_sample)) |> 
  mutate(closest_to_zero = pmin(abs(CI.Lower), abs(CI.Upper))) |>
  add_count(.feature) |> 
  filter(n < 5) |> 
  with_groups(Hypothesis, ~ .x |> arrange(desc(closest_to_zero)) |> dplyr::slice(1:50))

# Save the unknown ethnicities
pseudobulk_sample |> 
  select(-contains("PC"), -contains("tSNE"), -contains("UMAP")) |> 
  filter(ethnicity_groups == "Other/Unknown") |> 
  as("SingleCellExperiment") |> 
  zellkonverter::writeH5AD(file = "~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity_unknown_ethnicity_1_0_6.h5ad", compression = "gzip")
system("~/bin/rclone copy ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity_unknown_ethnicity_1_0_6.h5ad box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/")


pseudobulk_sample_for_PCA = pseudobulk_sample

pseudobulk_sample_for_PCA = 
  pseudobulk_sample_for_PCA |> 
  select(-contains("PC"), -contains("tSNE"), -contains("UMAP")) |> 
   _[summaries |> pull(.feature) |> unique(), , drop=FALSE] |> 
  filter(ethnicity_groups != "Other/Unknown") |> 
  tidybulk::reduce_dimensions(method = "PCA", .abundance = counts_adjusted_ethnicity, .dims = 20 ) |> 
  tidybulk::reduce_dimensions(method = "tSNE", .abundance = counts_adjusted_ethnicity, initial_dims = 10, .dims = 2)  |> 
  tidybulk::reduce_dimensions(method = "UMAP", .abundance = counts_adjusted_ethnicity, pca = 10, .dims = 2, calculate_for_pca_dimensions = NULL)  

pseudobulk_sample_for_PCA  = pseudobulk_sample_for_PCA |> filter(PC1 < 60)
pseudobulk_sample_for_PCA  = pseudobulk_sample_for_PCA |> filter(PC4 > -20)

pseudobulk_sample_for_PCA |> 
  as("SingleCellExperiment") |> 
  zellkonverter::writeH5AD(file = "~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity_1_0_6.h5ad", compression = "gzip")
system("~/bin/rclone copy ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity_1_0_6.h5ad box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/")

# pseudobulk_sample_for_PCA = 
#   zellkonverter::readH5AD(
#     file = "~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_sample_for_PCA_adjusted_ethnicity.h5ad",
#     use_hdf5 = TRUE, 
#     reader = "R"
#   )


colData(pseudobulk_sample)$sum = pseudobulk_sample |> assay("counts_adjusted_ethnicity") |> colSums()

pseudobulk_sample_for_PCA |> 
  pivot_sample() |> 
  left_join(cellNexus::get_metadata() |> distinct(sample_id, self_reported_ethnicity), copy = T) |>
  ggplot(aes(tSNE1, tSNE2, fill = ethnicity_groups)) +
  geom_point(shape = 21,  size = 0.8, stroke = 0) 

pseudobulk_sample_for_PCA |> 
  pivot_sample() |> 
  left_join(cellNexus::get_metadata() |> distinct(sample_id, self_reported_ethnicity), copy = T) |>
  ggplot(aes(UMAP1, UMAP2, fill = ethnicity_groups)) +
  geom_point(shape = 21,  size = 0.8, stroke = 0) 


pseudobulk_sample_for_PCA |> 
  pivot_sample() |> 
  select(contains("PCA"), everything()) |> 
GGally::ggpairs(columns = 19:38, ggplot2::aes(colour=`ethnicity_groups`))


library(plotly)
library(tidyomics)
pseudobulk_sample_for_PCA |>
  pivot_sample() |> 
  plot_ly(
    x = ~`tSNE1`,
    y = ~`tSNE2`,
    z = ~`tSNE3`,
    color = ~ ethnicity_groups
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

# x |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_chunk_062dd2621e3cb7e1.rds")
# system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_chunk_062dd2621e3cb7e1.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/")

tar_meta( store = glue::glue("/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets")) |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |>
  select(name, error)

tar_workspace(
  summary_6f59d31740151e3c, 
  store = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets",
  script = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets.R"
)



library(tidybayes)
library(brms)
library(magrittr)
library(ggallin)




fit = tar_read_raw("estimates_chunk_1ae83911a81cabae", store = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_6/_targets", branches = 1) |> 
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
# 10 estimates_chunk_f0801034799af6ad NA 
# 10 estimates_chunk_f0801034799af6ad NA 
# 10 estimates_chunk_f0801034799af6ad NA 


