

job::job({
  
  library(targets)
  
  tar_script({
    
    library(tidyverse)
    library(targets)
    library(tarchetypes)
    library(glue)
    library(qs)
    library(crew)
    library(crew.cluster)
  
    # # Set file path 
    # hdf5_path = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseudobulk_sample_is_immune"
    # metadata_path = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/cell_metadata_1_0_6_sccomp_input_counts.rds"
    # 
    tar_option_set(
      
      
      memory = "transient", 
      garbage_collection = 100, 
      storage = "worker", 
      retrieval = "worker", 
      error = "continue", 
      
      #cue = tar_cue(mode = "never"), 
      
      workspace_on_error = TRUE,
      format = "qs", 
      
      workspaces = "estimates_chunk_6e3989dd8c3e9257", 
      debug = "estimates_chunk",
      
      controller = crew_controller_group(
        
        crew_controller_slurm(
          name = "elastic",
          workers = 500,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_max = 7,
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
          crashes_max = 5,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = c(80, 160), 
            cpus_per_task = 5, 
            time_minutes = c(60*24, 60*24),
            verbose = T
          )
        ),
        crew_controller_slurm(
          name = "elastic_big_30_cores",
          workers = 150,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_max = 5,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = c(160), 
            cpus_per_task = 30, 
            verbose = T
          )
        ),
        crew_controller_slurm(
          name = "elastic_8_cores",
          workers = 100,
          tasks_max = 20,
          seconds_idle = 30,
          crashes_max = 5,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = c(10, 20, 40), 
            cpus_per_task = 8, 
            verbose = T
          )
        )
      )
      
      
    )
    
    
    #-----#
    # Functions
    #-----#  
  
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
    
    get_adjusted_matrix = function(summary_df, column_adjusted){
      
      column_adjusted = enquo(column_adjusted)
      
      m = 
        summary_df |>
        unnest(!!column_adjusted) |> 
        # dplyr::filter(analysis == "observed_proportion") |> 
        select(.feature, adjusted___Estimate, sample_id) |> 
        pivot_wider(names_from = sample_id, values_from = adjusted___Estimate) |> 
        tidybulk:::as_matrix(rownames = ".feature") |> 
        as("sparseMatrix")  |> 
        Matrix::Matrix(sparse = T)
      
      # Cap infinite
      max_rm_infinite = 
        m |> 
        _[!m |> is.infinite()] |> 
        quantile(0.999)
      
      m |> 
        _[m > max_rm_infinite] = 
        max_rm_infinite
      
      m |> 
        _[m < 0] = 
        0
      
      return(m)
    }
    
    check_rclone_installation = function(){
      rclone_path <- path.expand("~/bin/rclone")
      if (!file.exists(rclone_path)) {
        stop("rclone was not found in the expected location '~/bin/rclone'.")
      }
    }
    
    # Edit covariate
    source("https://gist.githubusercontent.com/stemangiola/8fe6c45b79dd95a200c0fd2314ec57d0/raw/14c753da54df886c29a02ac5fbf464684002956c/gistfile1.txt")
  
    monotone_ageing_hypothesis_testing = function(){
        
        cell_type <- cell_type_df$cell_type[[1]] %>% make.names()
        run_path = glue('{target_path}/V1_{cell_type}/')
        
        pseudobulk_sample_id = qs_read(glue('{run_path}_targets/objects/pseudobulk_sample_id')) 
        
        plan(multisession, workers = 16)  
        
        age_sex_summary <- 
          read_delim(glue('{run_path}_targets/meta/meta'), delim = '|') %>% 
          filter(name %>% str_starts('estimates_chunk')) %>%
          filter(type == 'branch' & is.na(error)) %>%
          # head(10) %>% 
          pull(name) %>%
          future_map_dfr(function(x, run_path) {
            
            tmp = qs_read(glue('{run_path}_targets/objects/{x}'))
            
            res = tmp$brms_fit[[1]]$fit %>% summary() %>% 
              as.data.frame() %>% 
              filter(
                str_detect(rownames(.), "age_bin|sex") &
                  str_detect(rownames(.), "^(b_|r_)")  
              ) %>% select(summary.mean) %>% t %>% as.data.frame()
            
            rownames(res) = tmp$.feature
            return(res)
          }, run_path = run_path, .progress = TRUE)
        
        saveRDS(age_sex_summary, file = glue('{run_path}age_sex_summary.rds'))
        
      }
  
    fit_to_age_monotonic_changes = function(fit) {
      
      full_bins <- c(
        "age_binInfancy",      "age_binChildhood",   "age_binAdolescence",
        "age_binYoungAdulthood","age_binMiddleAge",  "age_binSenior_60",
        "age_binSenior_70"
      )
      
      vars <- brms::variables(fit)            # all parameter names in the model
      
      ## a bin is present if its *population* coefficient exists
      present_bins <- full_bins[paste0("b_", full_bins) %in% vars]
      
      ## need at least two bins to make a split
      if (length(present_bins) < 2) {
        return(tibble())                      # empty result → nothing to contrast
      }
      
      ## helper to build one contrast for an arbitrary vector of bins
      build_contrast <- function(k, bins) {
        younger <- bins[1:k]
        older   <- bins[(k + 1):length(bins)]
        glue::glue("({paste(older,   collapse = ' + ')})/{length(older)} - ",
                   "({paste(younger, collapse = ' + ')})/{length(younger)} > 0")
      }
      
      ## ------------------------------------------------------------------
      ## 2. loop over split points *within the available bins*
      ## ------------------------------------------------------------------
      purrr::map_dfr(seq_len(length(present_bins) - 1), function(k) {
        
        h_txt <- build_contrast(k, present_bins)
        
        ## ---------------- total = population + random ----------------
        h_tot <- brms::hypothesis(
          fit, h_txt,
          scope = "coef", group = "tissue_groups"
        )$hypothesis %>% 
          tibble::as_tibble() %>% 
          dplyr::transmute(
            component    = "total",
            tissue       = Group,
            split_after  = present_bins[k],
            younger_bins = paste(present_bins[1:k],                 collapse = ","),
            older_bins   = paste(present_bins[(k + 1):length(present_bins)], collapse = ","),
            estimate     = Estimate,
            ci_lower     = CI.Lower,
            ci_upper     = CI.Upper,
            post_prob    = Post.Prob
          )
        
        ## ---------------- fixed = population-level only --------------
        h_fix <- brms::hypothesis(
          fit, h_txt,
          scope = "standard"
        )$hypothesis %>% 
          tibble::as_tibble() %>% 
          dplyr::transmute(
            component    = "fixed",
            tissue       = "population",
            split_after  = present_bins[k],
            younger_bins = paste(present_bins[1:k],                 collapse = ","),
            older_bins   = paste(present_bins[(k + 1):length(present_bins)], collapse = ","),
            estimate     = Estimate,
            ci_lower     = CI.Lower,
            ci_upper     = CI.Upper,
            post_prob    = Post.Prob
          )
        
        dplyr::bind_rows(h_tot, h_fix)
      })
    }
    
    ## helper: drop coefficients that are absent
    side_expr <- function(bins) paste(bins, collapse = " + ")
    
    has_param <- function(expr) expr != "0"                       # TRUE if at least one real term
    safe_hyp  <- purrr::possibly(brms::hypothesis, NULL)     
    
    build_contrast <- function(k, age_bins){
      glue("({side_expr(age_bins[(k+1):length(age_bins)])})/{length(age_bins)-k} - 
       ({side_expr(age_bins[1:k])})/{k} > 0")
    }
    
    #-----#
    # Pipeline
    #-----#
    
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
        cellchat_file,
        "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/cellNexus_lr_signaling_pathway_strength.duckdb",
        format = "file"
      ),
      tar_target(
        ethnicity_imputed,
        {
          check_rclone_installation()
          temp_path = tempdir()
          system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reports/ning/data/All_pseudobulk_1_0_6_ethnicity_imputed_colData.csv {temp_path}/"))
          
          read_csv(glue("{temp_path}/All_pseudobulk_1_0_6_ethnicity_imputed_colData.csv")) |> 
            select(sample_id, ethnicity_groups, ethnicity_groups_imputed = finalEthnicity_groups) |> 
            mutate(ethnicity_groups_imputed = ethnicity_groups_imputed |> str_replace("_imp$", "_imputed"))
          
        }, packages = c("glue", "readr", "dplyr", "stringr")
      ),
      tar_target(
        cellchat,
          duckdb::duckdb() |> 
          dbConnect(dbdir = cellchat_file) |> 
          tbl("lr_pathway_table") |> 
          filter(source == "plasma" | target == "plasma") |>  
          
          distinct(sample_id, source, target, pathway_name, interaction_weight, interaction_count, pathway_prob, annotation ) |> 
          inner_join(
            get_metadata() |> 
              distinct(
                sample_id, donor_id, dataset_id, tissue, age_days, 
                sex, self_reported_ethnicity, disease, assay, title, collection_id, 
                cell_type_unified_ensemble, cell_type, is_immune
              ) |> 
              edit_covariates() |> 
              distinct(sample_id, donor_id, dataset_id, tissue_groups, age_days, 
                       sex, ethnicity_groups, disease_groups, assay_groups, age_decade, age_bin) |> 
              
              # Add imputed ethnicities, and assign original if not present (cmposition and DE might have different samples because of filtering)
              left_join(
                ethnicity_imputed,
                by = join_by(sample_id, ethnicity_groups)
              ) |> 
              mutate(ethnicity_groups_imputed = if_else(ethnicity_groups_imputed |> is.na(), ethnicity_groups, ethnicity_groups_imputed)) |> 
              
              # Remove confounders of non interest
              tidybulk:::.resolve_complete_confounders_of_non_interest_df(dataset_id, assay_groups, disease_groups) |> 
              
              # Set intercept
              mutate(
                ethnicity_groups_imputed = fct_relevel(ethnicity_groups_imputed, "European"),
                assay_groups___altered = fct_relevel(assay_groups___altered, "10x Genomics 3"),
                disease_groups___altered = fct_relevel(disease_groups___altered, "Normal"),
                age_bin = fct_relevel(age_bin, "Senior_50"),
                age_decade = fct_relevel(age_decade, "5")
              ) , 
            copy = TRUE
          ) |> 
          as_tibble() |> 
        
          nest(data = -c(source, target, pathway_name)) |> 
          group_by(source, target, pathway_name) |> 
          tar_group(),
        iteration = "group",
        packages = c("dbplyr", "duckdb", "cellNexus"),
        resources = tar_resources(crew = tar_resources_crew("elastic_8_cores")),
        
      ),
    
      
      # estimates_chunk 
      # This target fits Bayesian models on chunks of the data. It processes each feature's data, handles missing values,
      # defines the model specification with priors, and runs the Bayesian inference using the brm function.
      tar_target(
        estimates_chunk, 
        
        cellchat |> mutate(brms_fit = map(data, ~ {
          
          data = .x |> 
            
            # Causes some problems
            filter(sex != "unknown") |> 
            
            droplevels() 
          
          # Skip if not enought samples
          if(data |> distinct(sample_id) |> nrow() < 100) return(NULL)
          
          # Manually revise data colnames to suit brms bug
          colnames(data) = colnames(data) |> stringr::str_replace_all("_+", "_")
          
          # # Check if dispersion estimation has failed
          # if(data |> pull(dispersion) |> unique() |> is.na()){
          #   warning("The dispersion calculation has failed. 1 is given as default prior.")
          #   data = data |> mutate(dispersion = 1)
          # }
          
          # Define the model formula
          formula <- bf(
            # mean model
            sqrt(interaction_weight) ~ 1 + age_bin*sex + disease_groups_altered +
              ethnicity_groups_imputed + assay_groups_altered +
              (1 | dataset_id_altered) +
              (1 + age_bin*sex  | tissue_groups), # + ethnicity_groups_imputed
            
            # dispersion model: use 'sigma' rather than 'sd'
            sigma ~ 1 # + disease_groups_altered + assay_groups_altered 
            # + ethnicity_groups_imputed #+ (1 | tissue_groups)
          )
          
          # HPC pipeline: param V2:
          prior = c(
            prior(student_t(3, i, s), class = Intercept),
            prior(student_t(3, 0, s), class = Intercept, dpar = "sigma"),
            prior(student_t(3, 0, 5), class = b)
            #,
            #prior(student_t(3, 0, 1), class = b, dpar = "sigma")
          ) |>
            substitute(env = list(
              i = mean(sqrt(data$interaction_weight)),
              s = sd(sqrt(data$interaction_weight)*2)
              
            )) |>
            eval()

          
          # # dynamically extract param from stan data
          # # code used from brm
          # bterms <- brmsterms(
          #   formula = brms:::validate_formula(
          #     formula, data = data, family = gaussian(),
          #     autocor = NULL, sparse = NULL, cov_ranef = NULL
          #   )
          # )
          # bframe <- brms:::brmsframe(bterms, data)
          # sdata <- brms:::.standata(
          #   bframe, data = data, prior = prior,
          #   data2 = NULL, stanvars = NULL, threads = NULL
          # )
          # 
          # Kc <- sdata$Kc
          # Kc_sigma <- sdata$Kc_sigma
          # M_1 <- sdata$M_1; N_1 <- sdata$N_1
          # M_2 <- sdata$M_2; N_2 <- sdata$N_2
          # M_3 <- sdata$M_3; N_3 <- sdata$N_3
          # 
          # inits <- lapply(1:chains, function(i) {
          #   list(
          #     # mean submodel
          #     b         = rnorm(Kc,       0, 5),
          #     Intercept = rnorm(1, mean(sqrt(data$interaction_weight)), 1.5),
          #     # sigma submodel
          #     b_sigma         = rnorm(Kc_sigma, 0, 2),
          #     Intercept_sigma = rnorm(1,            0, 1)
          #   )
          # })
          
          chains = 2
          
          brm(
            formula = formula,
            data = data,
            family = gaussian(),
            prior = prior,
            chains = chains,
            cores = pmin(as.numeric(parallelly::availableCores()), chains), 
            threads = threading(threads = (as.numeric(parallelly::availableCores()) / chains) |> floor()),
            warmup = 500, 
            refresh = 10,
            backend = "cmdstanr", 
            #sparse = TRUE,
            #save_model = glue("{external_directory}~/temp.rds"),
            #algorithm = "pathfinder",
            # sample_prior = TRUE, 
            #init = inits,
            iter = 700  # Increase iterations for better convergence
          )
          
        })) |> 
          
        # Drop data because it is withn the brms object
        select(-data), 
        pattern = map(cellchat),
        packages = c( "brms", "glue", "stringr", "dplyr", "purrr"),
        resources = tar_resources(crew = tar_resources_crew("elastic_8_cores"))
        #,  cue = tar_cue(mode = "never")
        
      ),
        
        tar_target(
          summaries,
          estimates_chunk |> 
            mutate(summary = map(
              brms_fit,
              ~ {
                if(.x |> is.null()) return(NULL)
                
                .x |> 
                  posterior::summarise_draws() |> 
                  rename(parameter = variable)
              }
            )) |> 
            select(-brms_fit),
         pattern = map(estimates_chunk),
         packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "tibble", "purrr", "posterior"),
         resources = tar_resources(crew = tar_resources_crew("elastic"))
          
        ),
      
      tar_target(
        hypothesis_age_monotonic,
        estimates_chunk |> 
          mutate(hypothesis_age_monotonic = map( brms_fit, fit_to_age_monotonic_changes )) |> 
          select(-brms_fit),
        pattern = map(estimates_chunk),
        packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "tibble", "purrr"),
        resources = tar_resources(crew = tar_resources_crew("elastic"))
        
      )
  
      
     # ## effect_removed 
     #  # This target generates adjusted model estimates by removing unwanted effects from the fitted Bayesian models,
     #  # thereby isolating the effects of interest. Here, nuisance covariates are set to NA and removed from the predictions.
     #  # This target produces adjusted estimates from the Bayesian models, removing unwanted effects while retaining 
     #  # the tissue group random effect, thus preserving variability associated with tissue-specific factors.
     #  tar_target(
     #    effect_removed, 
     #    estimates_chunk |> 
     #      mutate(brms_fit_adjusted_ethnicity = map(brms_fit, ~ .x |> remove_unwanted_effect(
     #        newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
     #        robust = TRUE, 
     #        re_formula = ~ 0
     #      ))) |> 
     #      mutate(brms_fit_adjusted_ethnicity_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
     #        newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
     #        robust = FALSE, correct_by_offset = FALSE,
     #        re_formula = ~ 0
     #      ))) |> 
     #      mutate(brms_fit_adjusted_tissue = map(brms_fit, ~ .x |> remove_unwanted_effect(
     #        newdata = .x$data |> mutate(assay_groups=NA, sex = NA, age_bin = NA, disease_groups = NA, ethnicity_groups = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
     #        robust = TRUE, 
     #        re_formula = ~ (1 | tissue_groups)
     #      ))) |> 
     #      mutate(brms_fit_adjusted_tissue_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
     #        newdata = .x$data |> mutate(assay_groups_altered=NA, ethnicity_groups = NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
     #        robust = FALSE, correct_by_offset = FALSE,
     #        re_formula = ~ (1 | tissue_groups)
     #      ))) |> 
     #      mutate(brms_fit_adjusted_ethnicity_estimate = map(brms_fit_adjusted_ethnicity, ~ {
     #        
     #        df = .x |> as_tibble()
     #        if (nrow(df) == length(pseudobulk_sample_id)){
     #          return(
     #            df |>
     #              select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
     #              mutate(sample_id = pseudobulk_sample_id)
     #          )
     #        }else{
     #          return(NULL)
     #        }
     #        
     #      })) |> 
     #      mutate(brms_fit_adjusted_ethnicity_new_estimate = map(brms_fit_adjusted_ethnicity_new, ~ {
     #        
     #        df = .x |> as_tibble()
     #        if (nrow(df) == length(pseudobulk_sample_id)){
     #          return(
     #            df |>
     #              select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
     #              mutate(sample_id = pseudobulk_sample_id)
     #          )
     #        }else{
     #          return(NULL)
     #        }
     #        
     #      })) |> 
     #      mutate(brms_fit_adjusted_tissue_estimate = map(brms_fit_adjusted_tissue, ~ {
     #        
     #        df = .x |> as_tibble()
     #        if (nrow(df) == length(pseudobulk_sample_id)){
     #          return(
     #            df |>
     #              select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
     #              mutate(sample_id = pseudobulk_sample_id)
     #          )
     #        }else{
     #          return(NULL)
     #        }
     #        
     #      })) |>
     #      mutate(brms_fit_adjusted_tissue_new_estimate = map(brms_fit_adjusted_tissue_new, ~ {
     #        
     #        df = .x |> as_tibble()
     #        if (nrow(df) == length(pseudobulk_sample_id)){
     #          return(
     #            df |>
     #              select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
     #              mutate(sample_id = pseudobulk_sample_id)
     #          )
     #        }else{
     #          return(NULL)
     #        }
     #        
     #      })) |>
     #      select(-brms_fit),
     #    
     #    pattern = map(estimates_chunk),
     #    packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
     #    resources = tar_resources(crew = tar_resources_crew("elastic"))
     #  ),
     
     # # adjusted_matrix -----
     # tar_target(
     #    adjusted_assay_ethnicity,
     #    get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_estimate),
     #    packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
     #    resources = tar_resources(
     #      crew = tar_resources_crew("elastic_big_30_cores")
     #    )
     #  ),
     #  
     #  tar_target(
     #    adjusted_assay_ethnicity_new,
     #    get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_new_estimate),
     #    packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
     #    resources = tar_resources(
     #      crew = tar_resources_crew("elastic_big_30_cores")
     #    )
     #  ),
     #  
     #  tar_target(
     #    adjusted_assay_tissue,
     #    get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_estimate),
     #    packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
     #    resources = tar_resources(
     #      crew = tar_resources_crew("elastic_big_30_cores")
     #    )
     #  ),
     #  
     #  tar_target(
     #    adjusted_assay_tissue_new,
     #    get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_new_estimate),
     #    packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
     #    resources = tar_resources(
     #      crew = tar_resources_crew("elastic_big_30_cores")
     #    )
     #  )
      
    ) # end ) of all target list
    
    
  },
  ask = FALSE, 
  script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets.R"
  )

  tar_make(
    # callr_function = NULL,
    script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets.R", 
    store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets", 
    reporter = "summary"
  )

})

library(targets)
tar_workspace(
  "hypothesis_age_monotonic_7e6b87ed1b81cb83",
  script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets.R", 
  store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets"
)

tar_read(
  hypothesis_age_monotonic,
  store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets", 
  branches = 1
)

tar_read(
  estimates_chunk,
  store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets", 
  branches = 1
)

tar_poll(store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets")

summary_and_convergence = 
  tar_read(
    summaries,
    store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets"
  ) |> 
  unnest(summary) |> 
  select(source, target, pathway_name, parameter, mean, median, rhat)

x = tar_read(
  hypothesis_age_monotonic,
  store ="/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/cellchat_brms_1_6_12/_targets"
  #, branches = 1
) |> 
  filter(map_int(hypothesis_age_monotonic, nrow) > 0) |>
  unnest(hypothesis_age_monotonic) |> 
  mutate(star = sign(ci_lower) == sign(ci_upper)) 
  
data_for_plot =  x |> 
  left_join(summary_and_convergence |>   filter(parameter == "b_Intercept") |> rename(intercept = mean) ) |> 
  
  # convergence
  left_join(
    summary_and_convergence |> 
      mutate(converged = rhat |> between(0.95, 1.05)) |> 
      filter(parameter |> str_detect("age"), parameter |> str_detect("sex", negate = TRUE), parameter |> str_detect("^b_|^r_")) |> 
      summarise(converged = all(converged), .by = c(source, target, pathway_name))
  ) |> 
  mutate(
    direction = if_else(source=="plasma", "out", "in"),
    other_cell_type  = if_else(source=="plasma", target, source)
  ) |> 
  filter(!other_cell_type %in% c("immune", "blood")) |> 
  filter(converged) |> 
  mutate(tissue = if_else(star, tissue, NA))

plotly::ggplotly(


  data_for_plot |> 
    filter(star) |> 
  ggplot(aes(intercept, estimate, label = glue("{pathway_name} {tissue}"))) +
    geom_point(color = "grey", shape = ".", data =  data_for_plot |> filter(!star)) +
    
  geom_point(aes(color = tissue), size = 0.5) +
    geom_hline(yintercept = 0) +
    facet_wrap(~other_cell_type)
)
