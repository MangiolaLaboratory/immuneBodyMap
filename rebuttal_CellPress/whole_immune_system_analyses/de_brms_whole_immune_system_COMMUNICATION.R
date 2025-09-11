

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
    
    # Helper (optional) to avoid repetition
    new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, backup = NULL) {
      crew_controller_slurm(
        name = name,
        workers = workers,
        crashes_max = crashes_max,
        seconds_idle = 30,
        options_cluster = crew_options_slurm(
          memory_gigabytes_required = mem_gb,
          cpus_per_task = 8,
          time_minutes = time_min
        ),
        backup = backup
      )
    }
    
    # Small → large, with fallbacks to the next size up
    elastic_160 <- new_elastic("elastic_160", 160, 60 * 24, workers = 8,  crashes_max = 2)
    elastic_80  <- new_elastic("elastic_80",   80,  60 * 4,  workers = 16, crashes_max = 1, backup = elastic_160)
    elastic_40  <- new_elastic("elastic_40",   40,  60 * 4,  workers = 24, crashes_max = 1, backup = elastic_80)
    elastic_20  <- new_elastic("elastic_20",   20,  60 * 4,  workers = 32, crashes_max = 1, backup = elastic_40)
    elastic_10  <- new_elastic("elastic_10",   10,  60 * 4,  workers = 48, crashes_max = 1, backup = elastic_20)
    elastic_5   <- new_elastic("elastic_5",     5, 60 * 4,  workers = 64, crashes_max = 6, backup = elastic_10)
    
    # Group for targets (small → large)
    controllers <- crew_controller_group(
      elastic_5, elastic_10, elastic_20, elastic_40, elastic_80, elastic_160
    )
    
    tar_option_set(
      
      
      memory = "transient", 
      garbage_collection = 100, 
      storage = "worker", 
      retrieval = "worker", 
      error = "continue", 
      
      #cue = tar_cue(mode = "never"), 
      
      workspace_on_error = TRUE,
      format = "qs", 
      
      workspaces = "estimates_chunk_633f6596569029e6", 
      debug = "estimates_chunk_8da8c24dbc0993bb",
      
      
      controller = controllers
      
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
    
    ## helper to build one contrast for an arbitrary vector of bins
    build_contrast <- function(k, bins, minimum_difference = 0) {
      
      if(minimum_difference == 0) sign = "="
      else if(minimum_difference < 0) sign = "<"
      else if(minimum_difference > 0) sign = ">"
      
      younger <- bins[1:k]
      older   <- bins[(k + 1):length(bins)]
      glue::glue("({paste(older,   collapse = ' + ')})/{length(older)} - ",
                 "({paste(younger, collapse = ' + ')})/{length(younger)} {sign} {minimum_difference}")
    }
    
    # 
    #  helper: absolute-difference test, returns a brms_hypothesis-like table
    # 
    hypothesis_abs <- function(fit, k, bins,
                               abs_threshold = 0.2,
                               scope = "coef",
                               group = "",
                               re_formula = NA,
                               resp = NULL) {
      
      browser()
      # build the two complementary algebraic strings
      h_pos <- build_contrast(k, bins,  abs_threshold)          # > +τ
      h_neg <- build_contrast(k, bins, -abs_threshold) |>
        stringr::str_replace(">", "<")                   # < –τ
      
      ## evaluate ------------------------------------------------------------
      tbl_pos <- brms::hypothesis(fit, h_pos,
                                  scope      = scope,
                                  group      = group,
                                  re_formula = re_formula,
                                  resp       = resp)$hypothesis |>
        tibble::as_tibble() |>
        dplyr::mutate(direction = "positive",
                      prob = Post.Prob)
      
      tbl_neg <- brms::hypothesis(fit, h_neg,
                                  scope      = scope,
                                  group      = group,
                                  re_formula = re_formula,
                                  resp       = resp)$hypothesis |>
        tibble::as_tibble() |>
        dplyr::mutate(direction = "negative",
                      prob = Post.Prob)
      
      long_tbl <- dplyr::bind_rows(tbl_pos, tbl_neg)
      
      if(!"Group" %in% colnames(long_tbl)) long_tbl = long_tbl |> mutate(Group = "population")
      
      ## pick winner & compute absolute evidence ----------------------------
      out <- long_tbl |>
        dplyr::group_by(Group) |>
        dplyr::mutate(
          prob_abs   = sum(prob),                       # P(|Δ| > τ)
          p_two_tail = 2 * min(prob),                   # doubled smaller tail
          p_two_tail = pmin(p_two_tail, 1)
        ) |>
        dplyr::slice_max(prob, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::transmute(
          Group, Hypothesis,
          Estimate, CI.Lower, CI.Upper,
          Post.Prob = p_two_tail,        # absolute tail probability
          direction                    # “positive” or “negative”
        )
      
      out
    }
    
    
    
    fit_to_age_monotonic_changes = function(fit) {
      
      if(fit |> is.null()) return(NULL)
      
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
      
      
      
      ## 2. loop over split points *within the available bins*
      seq_len(length(present_bins) - 1) |> 
        purrr::map_dfr(function(k) {
          
          h_txt <- build_contrast(k, present_bins)
          
          h_tot = 
            #   brms::hypothesis(
            #   fit, h_txt,
            #   scope = "coef", group = "tissue_groups"
            # )$hypothesis 
            
            # ---------------- total = population + random ----------------
          hypothesis_abs(
            fit, k, present_bins,
            abs_threshold = 0.001,
            scope = "coef", 
            group = "tissue_groups"
          )
          
          h_tot = h_tot %>% 
            tibble::as_tibble() %>% 
            dplyr::transmute(
              component    = "total",
              tissue       = Group,
              split_after  = present_bins[k],
              younger_bins = paste(present_bins[1:k],  collapse = ","),
              older_bins   = paste(present_bins[(k + 1):length(present_bins)], collapse = ","),
              estimate     = Estimate,
              ci_lower     = CI.Lower,
              ci_upper     = CI.Upper,
              post_prob    = Post.Prob,
              Hypothesis = Hypothesis,
              # p_two_tail = p_two_tail,
              direction = direction
            )
          
          
          # ---------------- fixed = population-level only --------------
          h_fix <- 
            hypothesis_abs(
              fit, k, present_bins,
              abs_threshold = 0.001,
              scope = "standard"
            )
          
          # brms::hypothesis(
          #   fit, h_txt,
          #   scope = "standard"
          # )$hypothesis
          
          h_fix = h_fix %>% 
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
              post_prob    = Post.Prob,
              Hypothesis = Hypothesis,
              # p_two_tail = p_two_tail,
              direction = direction
            )
          
          dplyr::bind_rows(h_tot, h_fix)
        })
    }
    
    ## helper: drop coefficients that are absent
    side_expr <- function(bins) paste(bins, collapse = " + ")
    
    has_param <- function(expr) expr != "0"                       # TRUE if at least one real term
    safe_hyp  <- purrr::possibly(brms::hypothesis, NULL)     
    
    
    prepare_database = function(tbl, ethnicity_imputed){
      tbl |> 
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
        ) 
    }
    
    get_pairs_to_consider = function(all_cell_types){
      myeloid_lymphoid_pairs <- tribble(
        ~source,           ~target,                 ~mechanism,                                                                              ~doi,
        
        # cDC → T (broad priming; costimulation & cytokines)
        "cdc",           "cd4 naive",         "Antigen presentation + CD80/86–CD28 costimulation (priming)",                          "10.1038/32588",
        "cdc",           "cd4 tcm",           "Recall priming / re-stimulation by cDC",                                               "10.1038/32588",
        "cdc",           "cd4 tem",           "Effector re-stimulation by cDC",                                                       "10.1038/32588",
        "cdc",           "cd4 fh em",         "ICOSL/IL-6-driven Tfh differentiation by APCs incl. DCs",                              "10.1016/j.immuni.2014.10.004",
        "cdc",           "cd4 th1 em",        "IL-12 skews Th1 during DC-T priming",                                                  "10.1084/jem.184.2.741",
        "cdc",           "cd4 th2 em",        "DC-guided Th2 under appropriate cues (e.g., OX40L/TSLP contexts)",                     "10.1038/32588",
        "cdc",           "cd4 th17 em",       "IL-23/IL-6 axis promotes Th17 from DC priming",                                        "10.1084/jem.20041257",
        "cdc",           "cd4 th1/th17 em",   "Mixed Th1/Th17 polarisation from DC cytokine milieu",                                  "10.1084/jem.20041257",
        "cdc",           "treg",              "Specialised DCs (e.g., CD103+ intestinal DCs) induce peripheral FoxP3+ Tregs",         "10.1093/intimm/dxae042",
        
        "cdc",           "cd8 naive",         "Cross-presentation primes CD8+ T cells",                                               "10.1038/35100512",
        "cdc",           "cd8 tcm",           "Recall responses via cross-presenting DCs",                                            "10.1038/35100512",
        "cdc",           "cd8 tem",           "Effector re-stimulation via cross-presentation",                                       "10.1038/35100512",
        "cdc",           "cytotoxic",         "Licensing/priming of cytotoxic T cells by cross-presentation",                         "10.1038/35100512",
        "cdc",           "tgd",               "DC antigen presentation & costimulation supports γδ T activation",                     "10.1038/32588",
        
        # pDC → lymphoid
        "pdc",           "nk",                "Type I IFN from pDC activates NK cells",                                               "10.1126/science.284.5421.1835",
        "pdc",           "b naive",           "pDC (±viral trigger) induce B-cell activation/differentiation via IFN-I + IL-6",      "10.1016/S1074-7613(03)00208-5",
        "pdc",           "b memory",          "pDC (±viral trigger) induce B-cell activation/differentiation via IFN-I + IL-6",      "10.1016/S1074-7613(03)00208-5",
        "pdc",           "plasma",            "pDC drive plasmablast/plasma-cell differentiation",                                    "10.1016/S1074-7613(03)00208-5",
        
        # Monocytes (CD14/CD16/“monocytic”) → lymphoid
        "cd14 mono",     "nk",                "IL-12/IL-18 family from myeloid cells activates NK (myeloid–NK cytokine crosstalk)",   "10.3389/fimmu.2021.739220",
        "cd16 mono",     "nk",                "IL-12/IL-18 family from myeloid cells activates NK (myeloid–NK cytokine crosstalk)",   "10.3389/fimmu.2021.739220",
        "monocytic",     "nk",                "IL-12/IL-18 family from myeloid cells activates NK (myeloid–NK cytokine crosstalk)",   "10.3389/fimmu.2021.739220",
        
        "cd14 mono",     "mait",              "Cytokine-driven (IL-12/IL-18) MAIT activation by monocytes/APCs",                      "10.3389/fimmu.2020.01014",
        "cd16 mono",     "mait",              "Cytokine-driven (IL-12/IL-18) MAIT activation by monocytes/APCs",                      "10.3389/fimmu.2020.01014",
        "monocytic",     "mait",              "Cytokine-driven (IL-12/IL-18) MAIT activation by monocytes/APCs",                      "10.3389/fimmu.2020.01014",
        
        # Macrophage → lymphoid
        "macrophage",    "cd4 th1 em",        "Macrophage IL-12 promotes Th1 differentiation/maintenance",                            "10.1084/jem.184.2.741",
        "macrophage",    "cd4 th17 em",       "Macrophage/monocyte IL-23/IL-6 supports Th17",                                         "10.1084/jem.20041257",
        "macrophage",    "cd8 naive",         "Myeloid APCs (incl. macrophages) prime/boost CD8 under some settings",                 "10.3389/fimmu.2013.00389",
        "macrophage",    "cd8 tem",           "Myeloid APCs (incl. macrophages) prime/boost CD8 under some settings",                 "10.3389/fimmu.2013.00389",
        "macrophage",    "cytotoxic",         "Boosting/maintenance of cytotoxic CD8 via macrophage APC function",                    "10.3389/fimmu.2013.00389",
        "macrophage",    "nk",                "Reciprocal cytokines (IL-12/IL-18 from myeloid ↔ IFN-γ/TNF-α from NK)",                "10.3389/fimmu.2021.739220",
        "macrophage",    "nkt",               "CD1d+ macrophages present glycolipids to iNKT (APC role)",                             "10.1016/j.coi.2007.03.014",
        "macrophage",    "mait",              "MR1/cytokine-dependent MAIT activation by macrophages",                                "10.3389/fimmu.2020.01014",
        
        # Granulocyte (neutrophil-centric) → lymphoid
        "granulocyte",   "cd4 tcm",           "Neutrophils can acquire APC features for memory CD4+",                                 "10.1182/blood-2016-10-744441",
        "granulocyte",   "cd4 tem",           "Neutrophils can acquire APC features for memory CD4+",                                 "10.1182/blood-2016-10-744441",
        "granulocyte",   "nk",                "Neutrophil–NK functional crosstalk in inflammation/cancer",                            "10.3389/fimmu.2020.570380",
        
        # Mast cell → lymphoid
        "mast",          "cd4 tem",           "Mast cells enhance T-cell activation via TNF/costims incl. OX40L",                     "10.3389/fimmu.2015.00394",
        "mast",          "treg",              "Bidirectional mast cell–Treg modulation (OX40–OX40L/cytokines)",                       "10.3389/fimmu.2015.00394",
        "mast",          "tgd",               "Direct immune synapse with γδ T cells (antiviral context)",                            "10.1172/JCI122530",
        "mast",          "b naive",           "Mast-cell mediators (e.g., IL-6) promote B-cell proliferation/differentiation",        "10.1182/blood-2009-10-250126",
        "mast",          "b memory",          "Mast-cell mediators (e.g., IL-6) promote B-cell proliferation/differentiation",        "10.1182/blood-2009-10-250126",
        
        # ILC (innate lymphoid) ↔ myeloid
        "ilc",           "macrophage",        "ILC2-derived IL-13 drives M2-like macrophage polarisation",                             "10.1155/2020/5018975",
        
        # NK ↔ myeloid (directionality kept explicit)
        "nk",            "cdc",               "DC↔NK licensing/activation loops (IL-12/IL-18, IL-15 trans-presentation)",              "10.3389/fimmu.2014.00159",
        "nk",            "cd14 mono",         "Myeloid–NK cytokine crosstalk (myeloid IL-12/18 → NK; NK IFN-γ/TNF-α → myeloid)",      "10.3389/fimmu.2021.739220",
        "nk",            "macrophage",        "Myeloid–NK cytokine crosstalk (myeloid IL-12/18 → NK; NK IFN-γ/TNF-α → myeloid)",      "10.3389/fimmu.2021.739220",
        "nk",            "pdc",               "pDC type I IFN promotes NK activation",                                                "10.1126/science.284.5421.1835",
        
        # MAIT / NKT ↔ myeloid
        "mait",          "cdc",               "MAIT cells can mature DCs via CD40L/GM-CSF (feedback)",                                "10.4049/jimmunol.1700615",
        "mait",          "monocytic",         "Cytokine-dependent activation by monocytes; MR1-dependent with macrophages",           "10.3389/fimmu.2020.01014",
        "nkt",           "cdc",               "iNKT rapidly license DCs and shape adaptive responses",                                 "10.1016/j.coi.2007.03.014",
        "nkt",           "b naive",           "iNKT provide cognate B-cell help (iNKTfh; IL-21-dependent)",                           "10.1038/ni.2172",
        "nkt",           "b memory",          "iNKT provide cognate B-cell help (iNKTfh; IL-21-dependent)",                           "10.1038/ni.2166",
        "nkt",           "plasma",            "iNKT help promotes antibody-secreting cell generation",                                 "10.1038/ni.2166"
      )
      
      
      # Define non-immune explicitly (as per earlier discussion)
      nonimmune <- c(
        "endocrine","endothelial","epithelial","fat","glial","muscle",
        "myoepithelial","neuron","pericyte","pneumocyte","progenitor",
        "renal","secretory","stromal", "reproductive", "mesothelial", "lens", "sensory", "epidermal", "cartilage", "liver", "bone"
      )
      
      # Everything else is immune
      immune <- all_cell_types |>  setdiff(nonimmune)
      
      # Cartesian product: all non-immune with all immune
      nonimmune_immune_pairs <- 
        tidyr::crossing(
          source = factor(nonimmune, levels = nonimmune),   # keep your order
          target = factor(immune,    levels = immune)
        ) |> 
        bind_rows(
          tidyr::crossing(
            target = factor(nonimmune, levels = nonimmune),   # keep your order
            source = factor(immune,    levels = immune)
          )
        )
      
      pairs_to_consider = 
        myeloid_lymphoid_pairs |> 
        bind_rows(nonimmune_immune_pairs)
      
      pairs_to_consider |> 
        select(source, target) |> 
        bind_rows(
          pairs_to_consider |> 
            select(source, target) |> 
            set_names(c("target", "source"))
        ) |> 
        distinct()
      
    }
    
    
    #-----#
    # Pipeline
    #-----#
    
    list(
      
      
      
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
      
      tar_target(cell_to_exclude, c("immune", "blood", "t", "erythrocyte", "other", "cytotoxic", "b", "monocytic", "dc", "t cd8", "t cd4")),
      tar_target(
        cell_pathway_combination,
        {
          comb = 
            duckdb::duckdb() |> 
            dbConnect(dbdir = cellchat_file)    |> 
            tbl("lr_pathway_table") |> 
            
            # drop cells I dont want
            filter(!source %in% cell_to_exclude & !target %in% cell_to_exclude) |>  
            
            distinct(source, target, pathway_name) |> 
            as_tibble() 
          
          comb |> 
            
            # Filter plausible pairs
            inner_join(
              get_pairs_to_consider(
                comb |> 
                  select(source, target) |> 
                  pivot_longer(everything()) |>
                  pull(value) |> 
                  unique()
              ) |> 
                select(source, target) |> 
                distinct()
            ) |> 
            
            group_by(source, target) |>
            arrange(source, target) |>                     # set a deterministic order if needed
            mutate(.chunk = ceiling(row_number() / 10)) |> # 1,2,3,... every 10 rows
            group_by(source, target, .chunk) |>
            tar_group()
        }, 
        iteration = "group",
        packages = c("dbplyr", "duckdb", "cellNexus", "tarchetypes"),
        resources = tar_resources(crew = tar_resources_crew("elastic_5"))
      ),
      # estimates_chunk 
      # This target fits Bayesian models on chunks of the data. It processes each feature's data, handles missing values,
      # defines the model specification with priors, and runs the Bayesian inference using the brm function.
      tar_target(
        estimates_chunk, 
        
        cell_pathway_combination |> mutate(brms_fit = pmap(list(source, target, pathway_name), \(s, t, p){
          
          con = 
            duckdb::duckdb() |> 
            dbConnect(dbdir = cellchat_file, read_only = TRUE)
          
          data =  
            con |> 
            tbl("lr_pathway_table") |> 
            
            # Filter data
            filter(source == s, target == t, pathway_name == p) |> 
            
            as_tibble()
          
          data = data |> 
            prepare_database(ethnicity_imputed) |> 
            
            # Causes some problems
            filter(sex != "unknown") |> 
            
            as_tibble() |> 
            droplevels() 
          
          dbDisconnect(con, shutdown = TRUE)
          
          # Skip if not enought samples
          if(data |> distinct(sample_id) |> nrow() < 100) return(NULL)
          
          # Manually revise data colnames to suit brms bug
          colnames(data) = colnames(data) |> stringr::str_replace_all("_+", "_")
          
          # # Check if dispersion estimation has failed
          # if(data |> pull(dispersion) |> unique() |> is.na()){
          #   warning("The dispersion calculation has failed. 1 is given as default prior.")
          #   data = data |> mutate(dispersion = 1)
          # }
          
          formula_chr = 
            "log10(pathway_prob) ~ 1 + age_bin*sex + disease_groups_altered + ethnicity_groups_imputed + assay_groups_altered +
              (1 | dataset_id_altered) +
              (1 + age_bin*sex  | tissue_groups)" # + ethnicity_groups_imputed"
          
          if(data |> distinct(disease_groups_altered) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("+ disease_groups_altered"))
          if(data |> distinct(ethnicity_groups_imputed) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("+ ethnicity_groups_imputed"))
          if(data |> distinct(assay_groups_altered) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("+ assay_groups_altered"))
          if(data |> distinct(dataset_id_altered) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("(1 | dataset_id_altered) +"))
          if(data |> distinct(sex) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("*sex"))
          
          # Define the model formula
          formula <- bf(
            # mean model
            as.formula(formula_chr),
            
            # dispersion model: use 'sigma' rather than 'sd'
            sigma ~ 1 # + disease_groups_altered + assay_groups_altered 
            # + ethnicity_groups_imputed #+ (1 | tissue_groups)
          )
          
          # HPC pipeline: param V2:
          prior = c(
            prior(student_t(3, i, s), class = Intercept),
            prior(student_t(3, 0, s), class = Intercept, dpar = "sigma"),
            prior(student_t(3, 0, 0.5), class = b),
            prior(exponential(1), class = sd),
            prior(lkj(2), class = cor)
            #,
            #prior(student_t(3, 0, 1), class = b, dpar = "sigma")
          ) |>
            substitute(env = list(
              i = mean(sqrt(data$interaction_weight)),
              s = sd(sqrt(data$interaction_weight)*2)
              
            )) |>
            eval()
          
          chains = 2
          
          brm(
            formula = formula,
            data = data,
            family = gaussian(),
            prior = prior,
            sample_prior = "yes",
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
            iter = 2400  # Increase iterations for better convergence
          )
          
        }, .progress = TRUE)), 
        pattern = map(cell_pathway_combination),
        packages = c( "brms", "glue", "stringr", "dplyr", "purrr"),
        resources = tar_resources(crew = tar_resources_crew("elastic_5"))
        , cue = tar_cue(mode = "never")
        
      ),
      
      # tar_target(
      #   estimates_chunk_limma, 
      #   
      #   cell_pathway_combination |>
      #     nest(pathways = pathway_name) |> 
      #     mutate(brms_fit = pmap(list(source, target, pathways), \(s, t, p){
      #     browser()
      #     con = 
      #       duckdb::duckdb() |> 
      #       dbConnect(dbdir = cellchat_file, read_only = TRUE)
      #     
      #     data =  
      #       con |> 
      #       tbl("lr_pathway_table") |> 
      #       
      #       # Filter data
      #       filter(source == s, target == t, pathway_name %in% p$pathway_name) |> 
      #       
      #       as_tibble()
      #     
      #     data = data |> 
      #       prepare_database(ethnicity_imputed) |> 
      #       
      #       # Causes some problems
      #       filter(sex != "unknown") |> 
      #       
      #       as_tibble() |> 
      #       droplevels() 
      #     
      #     dbDisconnect(con, shutdown = TRUE)
      #     
      #     # Skip if not enought samples
      #     if(data |> distinct(sample_id) |> nrow() < 100) return(NULL)
      #     
      #     # Manually revise data colnames to suit brms bug
      #     colnames(data) = colnames(data) |> stringr::str_replace_all("_+", "_")
      #     
      #     # # Check if dispersion estimation has failed
      #     # if(data |> pull(dispersion) |> unique() |> is.na()){
      #     #   warning("The dispersion calculation has failed. 1 is given as default prior.")
      #     #   data = data |> mutate(dispersion = 1)
      #     # }
      #     
      #     formula_chr = 
      #       "sqrt(interaction_weight) ~ 1 + age_bin*sex*tissue_groups + disease_groups_altered + ethnicity_groups_imputed + assay_groups_altered + dataset_id_altered" # + ethnicity_groups_imputed"
      #     
      #     if(data |> distinct(disease_groups_altered) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("+ disease_groups_altered"))
      #     if(data |> distinct(ethnicity_groups_imputed) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("+ ethnicity_groups_imputed"))
      #     if(data |> distinct(assay_groups_altered) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("+ assay_groups_altered"))
      #     if(data |> distinct(dataset_id_altered) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("(+ dataset_id_altered) +"))
      #     if(data |> distinct(sex) |> nrow() ==1) formula_chr = formula_chr |> str_remove_all(fixed("*sex"))
      #     
      #     data |> 
      #       as_SummarizedExperiment(.sample = sample_id, .transcript = pathway_name, .abundance = interaction_weight) |> 
      #       impute_missing_abundance(.formula = ~age_bins, force_scaling = TRUE )
      #       test_differential_expression(
      #         abundance = "interaction_weight",
      #         method = "edger_robust_likelihood_ratio",
      #         .formula = as.formula(formula_chr)
      #       )
      #     
      #     # Define the model formula
      #     formula <- bf(
      #       # mean model
      #       as.formula(formula_chr),
      #       
      #       # dispersion model: use 'sigma' rather than 'sd'
      #       sigma ~ 1 # + disease_groups_altered + assay_groups_altered 
      #       # + ethnicity_groups_imputed #+ (1 | tissue_groups)
      #     )
      #     
      #     # HPC pipeline: param V2:
      #     prior = c(
      #       prior(student_t(3, i, s), class = Intercept),
      #       prior(student_t(3, 0, s), class = Intercept, dpar = "sigma"),
      #       prior(student_t(3, 0, 5), class = b)
      #       #,
      #       #prior(student_t(3, 0, 1), class = b, dpar = "sigma")
      #     ) |>
      #       substitute(env = list(
      #         i = mean(sqrt(data$interaction_weight)),
      #         s = sd(sqrt(data$interaction_weight)*2)
      #         
      #       )) |>
      #       eval()
      #     
      #     chains = 2
      #     
      #     brm(
      #       formula = formula,
      #       data = data,
      #       family = gaussian(),
      #       prior = prior,
      #       chains = chains,
      #       cores = pmin(as.numeric(parallelly::availableCores()), chains), 
      #       threads = threading(threads = (as.numeric(parallelly::availableCores()) / chains) |> floor()),
      #       warmup = 500, 
      #       refresh = 10,
      #       backend = "cmdstanr", 
      #       #sparse = TRUE,
      #       #save_model = glue("{external_directory}~/temp.rds"),
      #       #algorithm = "pathfinder",
      #       # sample_prior = TRUE, 
      #       #init = inits,
      #       iter = 1400  # Increase iterations for better convergence
      #     )
      #     
      #   })), 
      #   pattern = map(cell_pathway_combination),
      #   packages = c( "brms", "glue", "stringr", "dplyr", "purrr", "tidybulk", "tidySummarizedExperiment"),
      #   resources = tar_resources(crew = tar_resources_crew("elastic_5"))
      #   , cue = tar_cue(mode = "never")
      #   
      # ),
      
      tar_target(
        hypothesis_age_monotonic,
        estimates_chunk |> 
          
          # Hypothesis
          mutate(hypothesis_age_monotonic = map( brms_fit, fit_to_age_monotonic_changes )) |> 
          
          # Simple summary
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
        packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "tibble", "purrr", "posterior", "stringr"),
        resources = tar_resources(crew = tar_resources_crew("elastic_5"))
        #, cue = tar_cue(mode = "never")
        
      )
      
      
      
    )
    
    
  },
  ask = FALSE, 
  script = "/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets.R"
  )
  
  tar_make(
    callr_function = NULL,
    script = "/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets.R", 
    store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets", 
    reporter = "verbose" # "balanced"
  )
  
})

library(targets)
library(tidyverse)
library(duckdb)

tar_poll(store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets")


tar_workspace(
  "estimates_chunk_04da6012b937653b",
  script = "/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets.R", 
  store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets"
)

tar_read_raw(
  "estimates_chunk_633f6596569029e6",
  store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets", 
  branches = 1
)

tar_read(
  cell_pathway_combination,
  store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets", 
  branches = 1
)

tar_read_raw(
  "estimates_chunk_80604f0079ad66c0",
  store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets", 
  branches = 1
)



job::job({
  
  # meta = 
  #   tar_meta(
  #     starts_with("hypothesis_age_monotonic_"), 
  #     store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets"
  #   ) |> 
  #   filter(!data |> is.na())
  
  data_for_plot = tar_read(
    hypothesis_age_monotonic,
    #meta = meta,
    store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_d/_targets"
  ) |> 
    rowid_to_column() |> 
    filter(!map_lgl(hypothesis_age_monotonic, is.null, .progress = TRUE)) |>
    filter(map_int(hypothesis_age_monotonic, nrow, .progress = TRUE) > 0) |> 
    filter(map_lgl(summary, 
                   ~ .x |> 
                     filter(parameter |> str_detect("b_age")) |> 
                     filter(!rhat |> between(0.9, 1.1)) |> 
                     nrow() == 0, 
                   .progress = TRUE)) |>
    select(-summary) |>
    unnest(hypothesis_age_monotonic) |> 
    
    left_join(
      duckdb::duckdb() |>
        dbConnect(dbdir =
                    tar_read(
                      cellchat_file,
                      store ="/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_c/_targets"
                    ),
                  read_only = TRUE
        )    |>
        tbl("lr_pathway_table") |> 
        summarise(
          sample_size = n(), 
          average_interaction_count = mean(interaction_count), 
          .by = c(source, target, pathway_name)
        ) , 
      copy = TRUE
    ) |> 
    
    filter(!source %in% c("immune", "blood", "t", "erythrocyte", "other") & !target %in% c("immune", "blood", "t", "erythrocyte", "other")) |>
    mutate(star = (1-post_prob)<0.05) |> # & abs(estimate) > 0.2) |> 
    
    # multitest adjustment
    mutate(BH = p.adjust(1-post_prob, method="BH"), .by = split_after) |> 
    mutate(star = BH<0.05) 
  
})

data_for_plot |> saveRDS("../rebuttal_CellPress/comunication_for_plot.rds")


data_for_plot = readRDS("rebuttal_CellPress/comunication_for_plot.rds")

library(tidyverse)




# 3. Apply them as a manual discrete scale

#head(10000) |>
# filter(source == "plasma") |>
# filter(pathway_name == "GAP") |>

source("/vast/scratch/users/mangiola.s/cellchat_brms_1_6_12_c/_targets.R")

# Plot for compensation

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(forcats)

plot_dat <-
  data_for_plot |>
  arrange(post_prob) |> 
  filter(star) |>
  
  # Filter plausible pairs
  inner_join(get_pairs_to_consider(
    data_for_plot |> 
      select(source, target) |> 
      pivot_longer(everything()) |>
      pull(value) |> 
      unique()
  )) |> 
  
  mutate(direction = if_else(estimate > 0, "up", "down")) |>
  
  # ---- keep per tissue (no averaging across tissues) ------------------- #
  summarise(
    max_effect =  (function(x) {
      if (length(x) == 0L) return(NA_real_)
      a <- abs(x)
      a[is.na(a)] <- -Inf        # drop NAs from contention
      idx <- which.max(a)
      if (is.infinite(a[idx])) NA_real_ else x[idx]
    })(estimate),             # still average within a tissue, if multiple rows exist
    n_tissues   = 1L,                         # <<< changed: one row == one tissue
    sample_size = sum(sample_size),
    .by         = c(source, target, pathway_name, tissue, direction)  # <<< changed: include tissue
  ) |>
  
  ## expand for faceting -------------------------------------------------
mutate(cell_type = map2(source, target, ~ tibble(cell_type = c(.x, .y)))) |>
  unnest(cell_type) |>
  
  summarise(
    n_raw       = n(),
    sample_size = mean(sample_size),
    .by         = c(cell_type, direction)
  ) |>
  
  ## totals per cell type ------------------------------------------------
mutate(tot = sum(n_raw), .by = cell_type) |>
  
  ## predictors ----------------------------------------------------------
mutate(
  log_sample        = log(sample_size),
  log_sample_scaled = as.numeric(scale(log_sample)),
  log_tot           = log(tot)
) |>
  
  ## adjust tot ----------------------------------------------------------
mutate(
  resid_tot    = resid(lm(log_tot ~ log_sample_scaled, data = cur_data())),
  tot_adjusted = exp(resid_tot + mean(log_tot))
) |>
  
  ## adjust n ------------------------------------------------------------
mutate(
  sign_n         = if_else(direction == "down", -1, 1),
  log_n          = log(n_raw),
  resid_n        = resid(lm(log_n ~ log_sample_scaled, data = cur_data())),
  n_raw_signed   = sign_n * n_raw,
  n_adjusted     = sign_n * exp(resid_n + mean(log_n))
) |>
  
  ## value to order by: total magnitude of adjusted counts ---------------
mutate(order_adj = sum(abs(n_adjusted)), .by = cell_type) |>
  
  ## pivot for plotting --------------------------------------------------
select(-n_raw) |>
  pivot_longer(
    c(n_raw_signed, n_adjusted),
    names_to  = "series",
    values_to = "n"
  ) |> 
  
  # Filter immune cell types
  filter(cell_type %in% c("non_immune", "cd8 naive", "cd16 mono", "cd4 tcm", "cd4 th17 em", "granulocyte", "cd4 th1/th17 em", "treg", "b memory", "b naive", "nk", "plasma", "cd4 th2 em", "mast", "cd4 th1 em", "cd8 tem", "mait", "tgd", "cdc", "cd4 fh em", "cd4 naive", "nkt", "macrophage", "cd8 tcm", "cd14 mono", "pdc", "ilc"))  


plot_overall = 
  ggplot(plot_dat, aes(x = n, y = fct_reorder(cell_type, order_adj))) +
  
  ## raw counts first (white fill, thin outline) -------------------------
geom_col(
  data      = subset(plot_dat, series == "n_raw_signed"),
  fill      = "white",
  aes(colour = direction),
  width     = 0.8,
  linewidth = 0.3
) +
  
  ## adjusted counts overlay (filled) ------------------------------------
geom_col(
  data   = subset(plot_dat, series == "n_adjusted"),
  aes(fill = direction),
  width  = 0.6,
  colour = NA
) +
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1", guide = "none") +
  labs(
    x    = "Count",
    y    = NULL,
    fill = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Single cell type plots
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(glue)
library(forcats)
library(patchwork)
library(RColorBrewer)

# # 1. Count how many cell-type levels you actually have
# n_types <- 
#   data_for_plot |>
#   filter(star) |>
#   distinct(other_cell_type) |> 
#   nrow()

# 2. Generate an “extended” Set2 of exactly that many colours
extended_cols <- colorRampPalette(brewer.pal(8, "Set2"))(30)

plot_df <-
  
  # Filter 100% of the rows with significant "star" and sample 5% of rows with non significant
  bind_rows(
    # Keep all significant rows (star = TRUE)
    data_for_plot |> filter(star),
    # Sample 5% of non-significant rows (star = FALSE)
    data_for_plot |> filter(!star) |> sample_frac(0.15)
  ) |>
  
  # Filter the row with max absolute effect within source, target, pathway_name, tissue
  filter(
    .by = c(source, target, pathway_name, tissue),
    abs(estimate) == max(abs(estimate))
  ) |> 
  
  mutate(sample_size = sum(sample_size), .by = c(source, target, pathway_name, tissue)) |> 
  
  
  
  # Filter plausible pairs
  inner_join(get_pairs_to_consider(
    data_for_plot |> 
      select(source, target) |> 
      tidyr::pivot_longer(everything()) |>
      pull(value) |> 
      unique()
  )) |> 
  
  ## expand for faceting -------------------------------------------------
mutate(cell_type = purrr::map2(source, target,
                               ~ tibble::tibble(cell_type = c(.x, .y)))) |>
  tidyr::unnest(cell_type) |>
  
  ## build one panel per focal cell -------------------------------------
tidyr::nest(data = -cell_type) |>
  mutate(
    cell_types_to_keep_for_dotplot = 
      c( 
        data_for_plot |> dplyr::distinct(source) |> dplyr::pull(source), 
        data_for_plot |> dplyr::distinct(target) |> dplyr::pull(target)
      ) |> 
      unique() |> 
      list()
  ) |> 
  mutate(
    cell_types_to_keep_for_dotplot =
      purrr::map2(
        cell_type,
        cell_types_to_keep_for_dotplot,
        ~ if (.x == "pdc") {
          c("cd4 th2 em", "cd8 naive", "cd4 fh em", "cd4 th1 em", "treg")
        } else if(.x == "nkt") { 
          c("endothelial", "plasma", "cdc", "cd4 th1 em")
        } else if(.x == "cdc") { 
          c("nkt", "mait", "cd4 th2 em", "cd4 th1/th17 em", "cd4 th17 em", "tgd")
        } else {
          .y
        }
      )
  ) |> 
  mutate(
    ## compute top partner by adjusted signed counts (reusable by any plot)
    top_partner = purrr::map2(
      data, cell_type,
      ~ {
        focal <- .y
        raw_df <-
          .x |>
          dplyr::filter(star) |>
          dplyr::mutate(
            direction   = dplyr::if_else(source == focal, "out", "in"),
            other_cell  = dplyr::if_else(source == focal, target, source),
            effect_sign = dplyr::if_else(estimate >= 0, "positive", "negative")
          ) |>
          dplyr::distinct(pathway_name, other_cell, direction, effect_sign, sample_size)
        
        if (nrow(raw_df) == 0) return(NA_character_)
        
        agg_df <- raw_df |>
          dplyr::group_by(other_cell, direction, effect_sign) |>
          dplyr::summarise(
            n_raw       = dplyr::n(),
            sample_size = mean(sample_size),
            .groups     = "drop"
          ) |>
          dplyr::mutate(
            sign_n = dplyr::if_else(effect_sign == "positive", 1, -1)
          )
        
        adj_df <- agg_df |>
          dplyr::mutate(
            log_sample        = log(sample_size),
            log_sample_scaled = as.numeric(scale(log_sample)),
            log_n             = log(n_raw)
          )
        
        if (nrow(adj_df) > 1 && !all(is.na(adj_df$log_sample_scaled))) {
          model             <- stats::lm(log_n ~ log_sample_scaled, data = adj_df)
          adj_df$resid_n    <- stats::resid(model)
          adj_df$n_adjusted <- exp(adj_df$resid_n + mean(adj_df$log_n))
        } else {
          adj_df$n_adjusted <- adj_df$n_raw
        }
        
        top_tbl <- adj_df |>
          dplyr::mutate(n_adjusted_signed = sign_n * n_adjusted) |>
          dplyr::summarise(total_abs = sum(abs(n_adjusted_signed)), .by = other_cell) |>
          dplyr::arrange(dplyr::desc(total_abs), other_cell)
        
        if (nrow(top_tbl) == 0) NA_character_ else top_tbl$other_cell[1]
      }
    )
  ) |>
  mutate(
    ## 1. Dot plot (colour by tissue) ------------------------------------
    dot_plot = purrr::pmap(
      list(data, cell_type, cell_types_to_keep_for_dotplot), ~{
        
        ..1 |>
          
          
          mutate(
            n_tissues = n_distinct(tissue),
            .by = c(source, target, pathway_name)
          ) |> 
          
          dplyr::mutate(other_cell = dplyr::if_else(source == ..2, target, source)) |>
          dplyr::distinct() |>
          
          # filter relevant cell types
          dplyr::filter(other_cell %in% ..3) |> 
          
          ggplot2::ggplot(ggplot2::aes(n_tissues, estimate)) +
          ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
          ggplot2::geom_point(ggplot2::aes(size = sample_size, colour = tissue)) +  # <<< changed
          ggrepel::geom_text_repel(
            ggplot2::aes(label = glue::glue("{source}\n{target}\n{pathway_name}")),
            colour = "grey20"
          ) +
          ggplot2::facet_wrap(~other_cell) + 
          ggplot2::theme_minimal(base_size = 14) +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::labs(
            colour = "Tissue",              # <<< changed
            size   = "Total sample\nsize",
            x      = "Number of tissues",
            y      = "Mean effect"
          ) +
          ggplot2::scale_colour_discrete(drop = FALSE) +   # <<< changed: generic discrete scale
          ggplot2::ggtitle(..2)
      }
    ),
    
    ## 2. Pyramid plot (unchanged logic; still aggregates across tissues) -
    pyramid_plot = purrr::map2(
      data, cell_type,
      ~ {
        focal <- .y
        
        raw_df <-
          .x |>
          filter(star) |>
          dplyr::mutate(
            direction   = dplyr::if_else(source == focal, "out", "in"),
            other_cell  = dplyr::if_else(source == focal, target, source),
            effect_sign = dplyr::if_else(estimate >= 0, "positive", "negative")
          ) |>
          dplyr::distinct(pathway_name, other_cell, direction, effect_sign, sample_size) |>
          dplyr::group_by(other_cell, direction, effect_sign) |>
          dplyr::summarise(
            n_raw       = dplyr::n(),
            sample_size = mean(sample_size),
            .groups     = "drop"
          ) |>
          dplyr::mutate(
            sign_n       = dplyr::if_else(effect_sign == "positive", 1, -1),
            n_raw_signed = sign_n * n_raw
          )
        
        adj_df <- raw_df |>
          dplyr::mutate(
            log_sample        = log(sample_size),
            log_sample_scaled = as.numeric(scale(log_sample)),
            log_n             = log(n_raw)
          )
        
        if (nrow(adj_df) > 1 && !all(is.na(adj_df$log_sample_scaled))) {
          model                 <- stats::lm(log_n ~ log_sample_scaled, data = adj_df)
          adj_df$resid_n        <- stats::resid(model)
          adj_df$n_adjusted     <- exp(adj_df$resid_n + mean(adj_df$log_n))
        } else {
          adj_df$n_adjusted     <- adj_df$n_raw
        }
        
        adj_df <- adj_df |>
          dplyr::mutate(
            n_adjusted_signed = sign_n * n_adjusted,
            order_stat        = sum(abs(n_adjusted)),
            .by               = other_cell
          ) |>
          dplyr::select(other_cell, direction, order_stat, n_raw_signed, n_adjusted_signed) |>
          tidyr::pivot_longer(
            cols      = c(n_raw_signed, n_adjusted_signed),
            names_to  = "series",
            values_to = "n"
          )
        
        ggplot2::ggplot(adj_df,
                        ggplot2::aes(x = n,
                                     y = forcats::fct_reorder(other_cell, order_stat),
                                     fill = direction)) +
          ggplot2::geom_col(
            data      = dplyr::filter(adj_df, series == "n_raw_signed"),
            fill      = "white",
            colour    = "grey20",
            linewidth = 0.3,
            width     = 0.8
          ) +
          ggplot2::geom_col(
            data      = dplyr::filter(adj_df, series == "n_adjusted_signed"),
            colour    = NA,
            width     = 0.6
          ) +
          ggplot2::scale_x_continuous(labels = abs, expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
          ggplot2::scale_fill_manual(values = c("out" = "#377eb8", `in` = "#e41a1c"),
                                     name   = "Direction",
                                     breaks = c("out", "in"),
                                     labels = c("Out (source)", "In (target)")) +
          ggplot2::labs(x = "Number of pathways", y = NULL) +
          ggplot2::theme_minimal(base_size = 14) +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::ggtitle(focal)
      }
    ),
    
    ## 3) Volcano plot (effect vs post_prob) ------------------------------
    volcano_plot = purrr::pmap(
      list(data, cell_type, cell_types_to_keep_for_dotplot), ~{
        ..1 |>
          dplyr::mutate(other_cell = dplyr::if_else(source == ..2, target, source)) |>
          dplyr::filter(other_cell %in% ..3) |>
          
          # map to the template's aesthetics
          dplyr::mutate(
            logFC    = estimate,
            # use (1 - posterior) as a p-value analogue; guard against zeros
            PValue   = post_prob,
            sig_flag = star
          ) |>
          
          ggplot2::ggplot(ggplot2::aes(x = logFC, y = PValue)) +
          #ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lower, xmax = ci_upper, colour = sig_flag)) +
          
          ggplot2::geom_point(ggplot2::aes(colour = sig_flag, size = sig_flag)) +
          ggplot2::scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
          ggplot2::scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black")) +
          ggplot2::scale_size_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.1)) +
          ggplot2::facet_wrap(~other_cell) +
          ggplot2::labs(
            title = paste("Volcano:", ..2),
            x     = "Effect (signed max |estimate|)",
            y     = "Posterior (shown as 1 - post_prob)"
          ) +
          ggplot2::theme_minimal(base_size = 14) +
          ggplot2::theme(legend.position = "bottom")
      }
    ),
    
    ## 4) Error bar plot for top partner (effect with CI per pathway) ------
    top_partner_errorbar_plot = purrr::pmap(
      list(data, cell_type, top_partner), ~{
        focal    <- ..2
        top_cell <- ..3
        
        if (is.na(top_cell) || length(top_cell) == 0) {
          return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(glue::glue("{focal}: no data")))
        }
        
        df_top <- ..1 |>
          dplyr::mutate(other_cell = dplyr::if_else(source == focal, target, source)) |>
          dplyr::filter(other_cell == top_cell, star) |>
          dplyr::summarise(
            abs_est = abs(estimate),
            estimate = estimate[which.max(abs_est)],
            ci_lower = ci_lower[which.max(abs_est)],
            ci_upper = ci_upper[which.max(abs_est)],
            .by = pathway_name
          ) |>
          dplyr::arrange(estimate)
        
        if (nrow(df_top) == 0) {
          return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(glue::glue("{focal} vs {top_cell}: no pathways")))
        }
        
        ggplot2::ggplot(df_top, ggplot2::aes(x = forcats::fct_reorder(pathway_name, estimate), y = estimate)) +
          ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.25, colour = "grey40") +
          ggplot2::geom_point(colour = "black", size = 1.2) +
          ggplot2::coord_flip() +
          ggplot2::labs(
            x = NULL,
            y = "Effect (estimate with CI)",
            title = glue::glue("{focal} vs {top_cell} — significant pathways")
          ) +
          ggplot2::theme_minimal(base_size = 14) +
          ggplot2::theme(legend.position = "bottom")
      }
    )
  )

## example assembly ------------------------------------------------------
combined_plot <-
  plot_df |>
  filter(cell_type %in% "cdc") |>
  mutate(
    panel = pmap(
      list(volcano_plot, pyramid_plot, top_partner_errorbar_plot),
      ~ ..1 + ..2 + ..3 + plot_layout(widths = c(4, 1, 1))   # 4-to-1 ratio
    )
  ) |>
  pull(panel) |>
  wrap_plots(ncol = 1, guides = "collect") &       # one row per cell-type, shared guides
  theme(legend.position = "none")  

plot_df |>
  pull(top_partner_errorbar_plot) |> 
  _[[1]]

library(dplyr)
library(stringr)

data_for_plot |>
  filter(
    (source == "cdc" & target == "nkt") | (source == "nkt" & target == "cdc"),
    str_detect(pathway_name, "CNTN")
  ) |>
  mutate(direction = if_else(source == "cdc", "cdc->nkt", "nkt->cdc")) |>
  transmute(
    pathway_name, tissue, direction,
    Hypothesis,           # from hypothesis_age_monotonic
    estimate, ci_lower, ci_upper,
    post_prob, star
  ) |>
  arrange(direction, tissue, pathway_name)
