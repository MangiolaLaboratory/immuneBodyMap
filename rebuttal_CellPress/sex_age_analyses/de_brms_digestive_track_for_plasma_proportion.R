
library(Seurat)
library(purrr)
library(magrittr)
library(cellNexus)
library(dplyr)
library(readr)
library(forcats)
#devtools::load_all("~/labHead/cellNexus/")

#--------------------#
# QUERY FOR DHARMESH
#--------------------#

library(readr)
library(forcats)
library(glue)

# result_directory = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1") 
system(glue("~/bin/rclone copy box_adelaide:/minh_immune_map_disease/disease_data_grouped_further.csv ./"))


# REMOVE OLD CACHE which is here get_default_cache_dir()
get_metadata() |> 
  edit_covariates(
    read_csv(glue("./disease_data_grouped_further.csv"))
  ) |>
  
  # filter( your filtering ) |>
  get_single_cell_experiment(atlas = "cellxgene")


#-----------------------#
# BRMS estimation
#-----------------------#




# result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_no_ncell_factor"
#result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_scaled_ncell_factor"
# result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_scaled_ncell_factor_MAFA"
# result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_scaled_ncell_factor_digestive_tract"
# result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_scaled_ncell_factor_digestive_tract_3_datasets"
# result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_scaled_ncell_factor_digestive_tract_multitissue"
result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3_plasma_cell_study_scaled_ncell_factor_digestive_tract_multi_celltype_all"

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
    garbage_collection = 500, 
    storage = "worker", 
    retrieval = "worker", 
    
    workspace_on_error = TRUE, workspaces = "effect_removed_e305355bab6b41ed",
    format = "qs",
    
   #  debug = "feature_df",
    
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
          cpus_per_task = 2, 
          time_minutes = c(60*24, 60*24),
          verbose = T
        )
      ),
      crew_controller_slurm(
        name = "elastic_big_10_cores",
        workers = 150,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_max = 5,
        options_cluster = crew_options_slurm(
          memory_gigabytes_required = c(80, 160), 
          cpus_per_task = 10, 
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
  
  age_bin <- function(age_days, sex) {
    # Convert age in days to age in years
    age_years <- age_days / 365.25
    
    # Initialise an empty vector to store the results
    age_bins <- vector("character", length(age_years))
    
    # Define average thresholds for "unknown" sex based on midpoint between male and female stages
    unknown_thresholds <- c(3, 13, 20, 38, 52)
    
    # Loop through each element to assign the appropriate bin based on sex and age
    for (i in seq_along(age_years)) {
      if (sex[i] == "male") {
        age_bins[i] <- dplyr::case_when(
          age_years[i] < 3 ~ "Infancy",
          age_years[i] < 13 ~ "Childhood",
          age_years[i] < 21 ~ "Adolescence",
          age_years[i] < 40 ~ "Young Adulthood",
          age_years[i] < 55 ~ "Middle Age",
          age_years[i] >= 55 ~ "Senior",
          TRUE ~ NA_character_
        )
      } else if (sex[i] == "female") {
        age_bins[i] <- dplyr::case_when(
          age_years[i] < 3 ~ "Infancy",
          age_years[i] < 13 ~ "Childhood",
          age_years[i] < 19 ~ "Adolescence",
          age_years[i] < 36 ~ "Young Adulthood",
          age_years[i] < 50 ~ "Middle Age",
          age_years[i] >= 50 ~ "Senior",
          TRUE ~ NA_character_
        )
      } else if (sex[i] == "unknown") {
        age_bins[i] <- dplyr::case_when(
          age_years[i] < unknown_thresholds[1] ~ "Infancy",
          age_years[i] < unknown_thresholds[2] ~ "Childhood",
          age_years[i] < unknown_thresholds[3] ~ "Adolescence",
          age_years[i] < unknown_thresholds[4] ~ "Young Adulthood",
          age_years[i] < unknown_thresholds[5] ~ "Middle Age",
          age_years[i] >= unknown_thresholds[5] ~ "Senior",
          TRUE ~ NA_character_
        )
      } else {
        stop("Each element of 'sex' must be either 'male', 'female', or 'unknown'.")
      }
    }
    
    return(age_bins)
  }
  
  edit_covariates = function(tbl, disease_tbl){
    
    
    ethnicity_grouped <- tribble(
      ~self_reported_ethnicity, ~ethnicity_groups,
      "unknown", "Other/Unknown",
      "European", "European",
      "Korean", "East Asian",
      "Asian", "East Asian",
      "Japanese", "Japanese",
      "African American", "African",
      "Hispanic or Latin American", "Hispanic/Latin American",
      "Singaporean Chinese", "East Asian",
      "Han Chinese", "East Asian",
      "Singaporean Indian", "South Asian",
      "Singaporean Malay", "Other/Unknown",
      "British", "European",
      "African", "African",
      "South Asian", "South Asian",
      "European American", "European",
      "East Asian", "East Asian",
      "American", "Other/Unknown",
      "African American or Afro-Caribbean", "African",
      "Oceanian", "Native American & Pacific Islander",
      "Jewish Israeli", "Middle Eastern & North African",
      "Chinese", "East Asian",
      "South East Asian", "Other/Unknown",
      "Greater Middle Eastern  (Middle Eastern or North African or Persian)", "Middle Eastern & North African",
      "Native American", "Native American & Pacific Islander",
      "Pacific Islander", "Native American & Pacific Islander",
      "Finnish", "European",
      "Bangladeshi", "South Asian",
      "Native American,Hispanic or Latin American", "Hispanic/Latin American",
      "Irish", "European",
      "Iraqi", "Middle Eastern & North African",
      "European,Asian", "European"
    )
    
    assay_data_grouped <- tribble(
      ~assay, ~assay_groups,
      "10x 3' v2", "10x Genomics 3",
      "10x 3' v3", "10x Genomics 3",
      "10x 5' v2", "10x Genomics 5",
      "10x 5' v1", "10x Genomics 5",
      "MARS-seq", "Plate based Technologies",
      "10x 3' transcription profiling", "10x Genomics 3",
      "10x 5' transcription profiling", "10x Genomics 5",
      "Smart-seq2", "Smart seq",
      "microwell-seq", "Microwell Technologies",
      "TruDrop", "TruDrop",
      "Drop-seq", "Drop based Technologies",
      "Seq-Well S3", "Microwell Technologies",
      "GEXSCOPE technology", "Other Technologies",
      "Seq-Well", "Microwell Technologies",
      "sci-RNA-seq", "Other Technologies",
      "10x 3' v1", "10x Genomics 3",
      "BD Rhapsody Whole Transcriptome Analysis", "Other Technologies",
      "BD Rhapsody Targeted mRNA", "Other Technologies",
      "CEL-seq2", "Plate based Technologies",
      "SPLiT-seq", "Other Technologies",
      "STRT-seq", "Plate based Technologies",
      "inDrop", "Drop based Technologies",
      "Smart-seq v4", "Smart seq",
      "ScaleBio single cell RNA sequencing", "Other Technologies"
    )
    
    
    disease_data_grouped <- tribble(
      ~disease, ~disease_groups,
      
      # Normal control
      "normal", "Normal",
      
      # Isolated Diseases
      "COVID-19", "COVID-19 related",
      "post-COVID-19 disorder", "COVID-19 related",
      "long COVID-19", "COVID-19 related",
      "glioblastoma", "Glioblastoma",
      "lung adenocarcinoma", "Lung Adenocarcinoma",
      "systemic lupus erythematosus", "Systemic Lupus Erythematosus",
      
      # Infectious and Immune-related Diseases (other than COVID-19)
      "Crohn disease", "Infectious and Immune-related Diseases",
      "Crohn ileitis", "Infectious and Immune-related Diseases",
      "pneumonia", "Infectious and Immune-related Diseases",
      "common variable immunodeficiency", "Infectious and Immune-related Diseases",
      "toxoplasmosis", "Infectious and Immune-related Diseases",
      "Plasmodium malariae malaria", "Infectious and Immune-related Diseases",
      "type 1 diabetes mellitus", "Infectious and Immune-related Diseases",
      "influenza", "Infectious and Immune-related Diseases",
      "chronic rhinitis", "Infectious and Immune-related Diseases",
      "periodontitis", "Infectious and Immune-related Diseases",
      "localized scleroderma", "Infectious and Immune-related Diseases",
      "lymphangioleiomyomatosis", "Infectious and Immune-related Diseases",
      "listeriosis", "Infectious and Immune-related Diseases",
      
      # Cancer (other than isolated cancers)
      "squamous cell lung carcinoma", "Cancer",
      "small cell lung carcinoma", "Cancer",
      "non-small cell lung carcinoma", "Cancer",
      "breast carcinoma", "Cancer",
      "breast cancer", "Cancer",
      "luminal B breast carcinoma", "Cancer",
      "luminal A breast carcinoma", "Cancer",
      "triple-negative breast carcinoma", "Cancer",
      "gastric cancer", "Cancer",
      "colorectal cancer", "Cancer",
      "colon sessile serrated adenoma/polyp", "Cancer",
      "follicular lymphoma", "Cancer",
      "B-cell acute lymphoblastic leukemia", "Cancer",
      "B-cell non-Hodgkin lymphoma", "Cancer",
      "acute myeloid leukemia", "Cancer",
      "acute promyelocytic leukemia", "Cancer",
      "plasma cell myeloma", "Cancer",
      "clear cell renal carcinoma", "Cancer",
      "nonpapillary renal cell carcinoma", "Cancer",
      "basal cell carcinoma", "Cancer",
      "colorectal neoplasm", "Cancer",
      "adenocarcinoma", "Cancer",
      "chromophobe renal cell carcinoma", "Cancer",
      "neuroendocrine carcinoma", "Cancer",
      "lung large cell carcinoma", "Cancer",
      "tongue cancer", "Cancer",
      "Wilms tumor", "Cancer",
      "pleomorphic carcinoma", "Cancer",
      "blastoma", "Cancer",
      
      # Neurodegenerative and Neurological Disorders
      "dementia", "Neurodegenerative and Neurological Disorders",
      "Alzheimer disease", "Neurodegenerative and Neurological Disorders",
      "Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
      "multiple sclerosis", "Neurodegenerative and Neurological Disorders",
      "Down syndrome", "Neurodegenerative and Neurological Disorders",
      "trisomy 18", "Neurodegenerative and Neurological Disorders",
      "frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
      "temporal lobe epilepsy", "Neurodegenerative and Neurological Disorders",
      "Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "amyotrophic lateral sclerosis 26 with or without frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
      
      # Respiratory Conditions
      "pulmonary fibrosis", "Respiratory Conditions",
      "respiratory system disorder", "Respiratory Conditions",
      "chronic obstructive pulmonary disease", "Respiratory Conditions",
      "cystic fibrosis", "Respiratory Conditions",
      "interstitial lung disease", "Respiratory Conditions",
      "hypersensitivity pneumonitis", "Respiratory Conditions",
      "non-specific interstitial pneumonia", "Respiratory Conditions",
      "aspiration pneumonia", "Respiratory Conditions",
      "pulmonary emphysema", "Respiratory Conditions",
      "pulmonary sarcoidosis", "Respiratory Conditions",
      
      # Cardiovascular Diseases
      "myocardial infarction", "Cardiovascular Diseases",
      "acute myocardial infarction", "Cardiovascular Diseases",
      "dilated cardiomyopathy", "Cardiovascular Diseases",
      "heart failure", "Cardiovascular Diseases",
      "arrhythmogenic right ventricular cardiomyopathy", "Cardiovascular Diseases",
      "congenital heart disease", "Cardiovascular Diseases",
      "non-compaction cardiomyopathy", "Cardiovascular Diseases",
      "cardiomyopathy", "Cardiovascular Diseases",
      "heart disorder", "Cardiovascular Diseases",
      
      # Metabolic and Other Disorders
      "type 2 diabetes mellitus", "Metabolic and Other Disorders",
      "chronic kidney disease", "Metabolic and Other Disorders",
      "digestive system disorder", "Metabolic and Other Disorders",
      "primary sclerosing cholangitis", "Metabolic and Other Disorders",
      "gastritis", "Metabolic and Other Disorders",
      "acute kidney failure", "Metabolic and Other Disorders",
      "tubular adenoma", "Metabolic and Other Disorders",
      "benign prostatic hyperplasia", "Metabolic and Other Disorders",
      "opiate dependence", "Metabolic and Other Disorders",
      "gingivitis", "Metabolic and Other Disorders",
      "hyperplastic polyp", "Metabolic and Other Disorders",
      "clonal hematopoiesis", "Metabolic and Other Disorders",
      "epilepsy", "Metabolic and Other Disorders",
      "age related macular degeneration 7", "Metabolic and Other Disorders",
      "kidney benign neoplasm", "Metabolic and Other Disorders",
      "malignant pancreatic neoplasm", "Metabolic and Other Disorders",
      "cataract", "Metabolic and Other Disorders",
      "macular degeneration", "Metabolic and Other Disorders",
      "hydrosalpinx", "Metabolic and Other Disorders",
      "tubulovillous adenoma", "Metabolic and Other Disorders",
      "gastric intestinal metaplasia", "Metabolic and Other Disorders",
      "Barrett esophagus", "Metabolic and Other Disorders",
      
      # Other Diseases
      "injury", "Other Diseases",
      "anencephaly", "Other Diseases",
      "primary biliary cholangitis", "Other Diseases",
      "keloid", "Other Diseases",
      "kidney oncocytoma", "Other Diseases",
      "respiratory failure", "Other Diseases",
      "pilocytic astrocytoma", "Other Diseases"
    )
    
    
    disease_data_grouped = 
      disease_data_grouped |> 
      select(-disease_groups) |> 
      left_join(disease_tbl) |> 
      mutate(disease_groups = if_else(disease_groups |> is.na(), "other", disease_groups))
    
    age_bin_table = 
      tbl |> 
      distinct(age_days, sex) |> 
      filter(!age_days |> is.na()) |> 
      mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |> 
      as_tibble() |> 
      mutate(age_bin = age_bin(age_days, sex))
    
    tbl |> 
      # TECH
      left_join(assay_data_grouped, copy=TRUE) |> 
      
      # DISEASE
      left_join(disease_data_grouped, copy=TRUE) |> 
      
      
      # TEMPORARY. de-group pancreas and liver
      mutate(tissue_groups = case_when(
        
        tissue %in% c("gallbladder") ~ "gallbladder",
        tissue %in% c("pancreas", "exocrine pancreas") ~ "pancreas",
        tissue %in% c("liver", "caudate lobe of liver", "hepatic cecum" ) ~ "liver",
        TRUE ~ tissue_groups
      )) |> 
      
      # SEX edit
      mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |> 
      
      # Age
      filter(age_days > 365) |> 
      left_join(age_bin_table, copy=TRUE) |> 
      
      # ETHNICITY
      left_join(ethnicity_grouped, copy=TRUE) |> 
      
      dplyr::select(cell_id, atlas_id, sample_id, donor_id, dataset_id, file_id_cellNexus_single_cell, title, collection_id, age_days, age_bin, sex, ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble, cell_type, disease_groups) |> 
      as_tibble() |> 
      
      # Set intercept
      mutate(
        ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
        assay_groups = fct_relevel(assay_groups, "10x Genomics 3"),
        disease_groups = fct_relevel(disease_groups, "Normal"),
        age_bin = fct_relevel(age_bin, "Adolescence")
      ) |>
      
      # Center based on adolescence
      mutate(age_days_scaled = age_days  |> scale(center = 15*365) |> as.numeric()) 
    
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
      sce,
      # load data and create pseudobulks
      {
        DelayedArray::setAutoBlockSize(size=1e9)
        
        get_metadata() |>
          filter(
            cell_type_unified_ensemble %in% c("epithelial", "plasma", "stromal"),
            feature_count > 600
            # ,
            # dataset_id %in% c(
            #   "7bb64315-9e5a-41b9-9235-59acf9642a3e",
            #   "e40c6272-af77-4a10-9385-62a398884f27",
            #   "d6dfdef1-406d-4efb-808c-3c5eddbfe0cb"
            # )
          ) |>
          edit_covariates(
            read_csv(glue("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/disease_data_grouped_further.csv"))
          ) |>
          filter(
            age_bin %in% c("Middle Age", "Senior"),
            disease_groups == "Normal",
            tissue_groups %in% c("large intestine", "stomach", "oesophagus"),
            
          ) |>
          as_tibble() |>
          droplevels() |>
          get_single_cell_experiment() |> 
          
          HDF5Array::saveHDF5SummarizedExperiment(
            "~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/plasma_sce_hdf5",
            as.sparse = TRUE,
            replace = TRUE,
            verbose = TRUE
          )
        
      } ,
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "tidySingleCellExperiment", "magrittr", "tibble", "glue", "cellNexus", "readr", "forcats", "scuttle", "BiocParallel", "purrr"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big")),
      memory = "persistent"
    ),
    
    tar_target(
      plasma_count,
        tar_read(input_relative, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_6_more_significant_figures/_targets") |>
        with_groups(sample_id, ~ .x |> mutate(plasma_prop = n/sum(n))) |> 
        dplyr::select(plasma_prop, sample_id, tissue_groups, cell_type_unified_ensemble) |> 
        filter(cell_type_unified_ensemble == "plasma") |> 
        dplyr::select(sample_id, tissue_groups, plasma_prop),
      packages = c("tidyverse", "targets"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big")),
      memory = "persistent"
    ),
    tar_target(
      pseudobulk_sample,
      {
        
        
        # system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/dharmesh_shared_mix/se_age_sex_tissues_epithelial.rds ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/")
        # sce = readRDS("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/se_age_sex_tissues_epithelial.rds")
        # 
        # 
        # system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/dharmesh_shared_mix/se_age_sex_tissues_epithelial_gut.rds ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/")
        # sce = readRDS("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/se_age_sex_tissues_epithelial_gut.rds")
      
        # colnames(sce)  = sce |> pull(sample_id)
        # 
        # sce = sce |> filter(ncells >=10)
        # 
        # sce = sce |> filter(tissue_groups %in% c("large intestine", "stomach", "oesophagus"))
        # 
        # se  = 
        #   sce |> 
        #   as("SummarizedExperiment")
        # 
        # # Filter for digestive tract
        # 
        # rownames(se) = rownames(sce) 
        # rownames(assay(se)) = rownames(sce) 
        
        # # Plasma prop
        # system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_estimates_1_0_6/estimates_age_bins___L3.rds ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/")
        # L3 = readRDS("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/estimates_age_bins___L3.rds")
        # plasma_prop = L3 |> sccomp_predict(~ age_bin*sex + (1 + age_bin*sex | tissue_groups))
        # plasma_prop = plasma_prop |> dplyr::filter(L3 == "plasma")
        # plasma_prop |> saveRDS("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/plasma_prop.rds")
        plasma_prop = readRDS("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/plasma_prop.rds")

        

        
         se = 
           sce |> 
           aggregateAcrossCells(
           ids = sce |>
             colData() |>
             as_tibble() |>
             dplyr::select(sample_id, cell_type_unified_ensemble, dataset_id, tissue_groups, age_bin, sex, ethnicity_groups, assay_groups) |>
             as("DataFrame"), 
           statistics = "sum", 
           BPPARAM = MulticoreParam(10)
          ) |>
           
           # Predicted proportions
           left_join(
             plasma_prop |> dplyr::rename(prop_plasma_predicted = proportion_mean) |> distinct(sample_id, prop_plasma_predicted)
           ) |>
           
           # Observed proportions
           left_join(plasma_count |> dplyr::rename(prop_plasma_observed = plasma_prop)) |> 
           mutate(prop_plasma_observed_arcsin_sqrt = asin(sqrt(prop_plasma_observed)) |> scale() |> as.numeric()) |> 
           
           
           # Some samples do not have prediction because they do not have immune cells within
           filter(!prop_plasma_observed |> is.na()) |> 
           as("SummarizedExperiment")
        
      rownames(se) = rownames(sce)
      se |> assay() |> rownames() = rownames(sce)
      #colnames(se)  = se |> tidybulk::pivot_sample() |> unite("pseudobulk_sample", sample_id, cell_type_unified_ensemble, tissue_groups) |> pull(pseudobulk_sample)
         
      se = se |>  filter(ncells > 10) 
      
        # Compute mean library size
        mean_library_size <- se |>
          assay("counts") |>
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
        
        # This technology is an outlier
        se = se |> 
          filter(assay_groups != "TruDrop") |> 
          filter(tissue_groups %in% c("large intestine", "stomach", "oesophagus")) |> 
          filter(age_bin %in% c("Middle Age", "Senior"))
        
        se = 
          se |> 
          
          # Get scaling factor 
          # DHARMESH PROB THE ONLY LIKE TO ADD
          identify_abundant(design = 
                                            se |>

                                            # Discretise the age for the following operation
                                            mutate(is_old_individual = age_days > 50*365) |>
                                            resolve_complete_confounders_of_non_interest(tissue_groups, cell_type_unified_ensemble, sex, is_old_individual) |>
                                            colData() |>
                                            droplevels() |>
                                            model.matrix(~ cell_type_unified_ensemble___altered + sex___altered + is_old_individual___altered, data = _  ),
                                          minimum_counts = 100
                            ) |> 
          scale_abundance(method = "TMMwsp", reference_sample = reference_sample) |> 
          mutate(offset = log(1/multiplier)) |> 
          
          # Resolve confounders
          group_split(cell_type_unified_ensemble) |> 
          map(tidybulk:::resolve_complete_confounders_of_non_interest,assay_groups, dataset_id) |> 
          do.call(cbind, args=_)
        
      se = se |> 
        mutate(
          prop_plasma_predicted_logit = prop_plasma_predicted |> boot::logit(),
          prop_plasma_predicted_logit_scaled = prop_plasma_predicted_logit |> scale() |> as.numeric(),
          ncells_scaled = ncells |> scale() |> as.numeric(),
          ncells_log_scaled = ncells |> log() |> scale() |> as.numeric()
        ) |> 
        
        # allow for interaction grouping
        unite("dataset_tissue_assay", dataset_id___altered, tissue_groups, assay_groups___altered, sep = "___", remove = FALSE)
        
        # As we have large sample size, we leave the estimation to gene-by-gene    
        # # Add dispersion
        # rowData(se)  = 
        #   rowData(se) |> 
        #   as_tibble(rownames = ".feature") |> 
        #   left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |> 
        #   data.frame(row.names = ".feature") |> DataFrame()
        

      
        se
        
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "tidySingleCellExperiment", "magrittr", "tibble", "glue", "cellNexus", "readr", "forcats", "scuttle", "BiocParallel", "purrr"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big_10_cores")),
      memory = "persistent"
    ),
    
    # Split in gene chunks
    tar_target(
      feature_df, 
      pseudobulk_sample |> 
        distinct(.feature, cell_type_unified_ensemble)|> 
        expand_grid(df = list(tribble(
          ~ analysis, ~formula, 
          "observed_proportion", c(
            "counts ~ 1 + offset(offset) + prop_plasma_observed_arcsin_sqrt + 
              (1 | dataset_tissue_assay)" ,
            "shape ~ 1 "
          ),
          "age_sex_interaction", c("counts ~ 1 + offset(offset) + age_bin*sex + assay_groups___altered +
              (1 | tissue_groups)" ,
            "shape ~ 1"
          )))
        ) |> 
        unnest(df) |> 
        group_by(.feature, cell_type_unified_ensemble, analysis, formula) |> 
        tar_group(), 
      iteration = "group",
      packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr", "brms"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      memory = "persistent"
    ),
    tar_target(
      se_df, 
      feature_df |> 
        mutate(
          se = map2(.feature, cell_type_unified_ensemble, ~ pseudobulk_sample[.x, , drop=FALSE] |> filter(cell_type_unified_ensemble == .y)
      )) |> 
        
        # Fit
        mutate(brms_fit = map2(se, formula, ~ {
          
          data = 
            .x |>
            as_tibble() |> 
            mutate(counts = counts |> as.integer()) |> 
            droplevels()
          
          n_NAs = data |> filter(counts |> is.na()) |> nrow()
          if(n_NAs > 0){
            warning(glue("You have {n_NAs} NAs in counts. They have been filtered out"))
            data = 
              data |> 
              filter(!counts |> is.na())
          }
          
          prior = c(
            prior(normal(i, 5), class = Intercept),
            prior(normal(0, 2), class = Intercept, dpar = shape),
            prior(normal(0, 5), class = b)
            # ,
            # prior(normal(0, 2), class = b, dpar = shape)
          ) |> 
            substitute(env = list(i = mean(log1p(data$counts / exp(data$offset))))) |> 
            eval()
          
          chains = 2
          inits <- list(Intercept = mean(log1p(data$counts / exp(data$offset))))
          inits <- replicate(chains, inits, simplify = FALSE)
          
          
          brm(
            formula =   bf(as.formula(.y[1]), as.formula(.y[2])),
            data = data,
            family = zero_inflated_negbinomial(),
            prior = prior,
            chains = chains,
            cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)), #, threads = 2,
            warmup = 600, 
            refresh = 10,
            backend = "cmdstanr", 
            #sparse = TRUE,
            #save_model = glue("{external_directory}~/temp.rds"),
            #algorithm = "pathfinder",
            init = inits,
            iter = 2000  # Increase iterations for better convergence
          )
          
        }))  , 
      pattern = map(feature_df),
      packages = 
        c( "brms", "glue") |> 
        c( "brms", "glue", "dplyr", "purrr", "SummarizedExperiment", "tidySummarizedExperiment", "tidySingleCellExperiment") |> 
        c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr") |> 
        c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(
        crew = tar_resources_crew("elastic")
        ),
      error = "continue"
      # format = "parquet"
    ),
    
    tar_target(
      summary_df, 
      se_df |> 
      
        
        # Summary
        # Random
        mutate(random_plus_fixed_effect_of_plasma_proportion = map(brms_fit, ~ .x |> coef(
          # Median instead and mad of mean and sd
          robust=TRUE
        ) |>
          _[[1]]  |> 
          # _[,,"prop_plasma_observed_arcsin_sqrt"] |> 
          as_tibble(rownames = "tissue_groups") |> 
          
          # Reshape + significance
          pivot_longer(-tissue_groups) |> 
          mutate(name = name |> str_replace("Q2\\.5", "Q2_5") |> 
                   str_replace("Q97\\.5", "Q97_5") |> 
                   str_replace("Est\\.Error", "Est_Error")
          ) |> 
          separate(name, c("stat", "parameter"), sep="\\.", extra="merge") |> 
          pivot_wider(names_from = "stat", values_from = "value") |> 
          mutate(significant = (Q2_5 * Q97_5) > 0)
        )) |>
        
        # Fixed
        mutate(fixed_effects = map(brms_fit, ~ .x |> fixef(
          # Median instead and mad of mean and sd
          robust=TRUE
        ) |> 
          as_tibble(rownames = "coefficient") |> 
 
          mutate(significant = (Q2.5 * Q97.5) > 0)  
        
        )) |> 
        
        mutate(Rhat = map_dbl(brms_fit, 
                              ~ summary(.x)$fixed |> 
                                as_tibble(rownames = "par") |> 
                                filter(par |> str_detect("plasma|sex|age")) |> 
                                pull(Rhat) |>
                                max()
        )) |> 
        
        # Adjust
        mutate(brms_fit_adjusted = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups=NA, ethnicity_groups = NA, ncells_log_scaled = NA), 
          robust = FALSE, correct_by_offset = FALSE,
          re_formula = ~ (1 + prop_plasma_logit_scaled  | tissue_groups) 
        ))) |> 
        select(-se, -brms_fit), 
      pattern = map(se_df),
      packages = 
        c( "brms", "glue") |> 
        c( "brms", "glue", "dplyr", "purrr", "SummarizedExperiment", "tidySummarizedExperiment", "tidySingleCellExperiment") |> 
        c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr") |> 
        c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(
        crew = tar_resources_crew("elastic")
      ),
      error = "continue"
      # format = "parquet"
    ),
    
    tar_target(
      adjusted_assay,
      {

        plan(callr, workers = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)))
        
        m = 
          summary_df |>
          select(1, 2, 3, 4, brms_fit_adjusted) |> 
          # filter(map_int(brms_fit_adjusted, nrow) == 4926 ) |> 
          mutate(brms_fit_adjusted = future_map2(brms_fit_adjusted, cell_type_unified_ensemble, ~ {
            library(tidySummarizedExperiment)
            .x |> 
                                            select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                                            mutate(sample_id = pseudobulk_sample |> dplyr::filter(cell_type_unified_ensemble == .y) |> colnames())
            }, 
                                          .progress = TRUE
          )) |>
          unnest(brms_fit_adjusted) |> 
          dplyr::filter(analysis == "observed_proportion") |> 
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
      },
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "future.callr", "furrr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
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

tar_read(estimates_chunk, store = glue::glue("{result_directory}/_targets"), branches = 1)

tar_meta(store = glue::glue("{result_directory}/_targets")) |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  dplyr::select(name, error)


tar_workspace(
  adjusted_assay, 
  script = glue::glue("{result_directory}/_targets.R"),
  store = glue::glue("{result_directory}/_targets")
)

# Effects
meta = tar_meta(starts_with("se_df_"), store = glue::glue("{result_directory}/_targets")) |> filter(!data |> is.na())

tar_read(summary_df, store = glue::glue("{result_directory}/_targets")    ) |> 

  filter(analysis == "observed_proportion") |> 
  select(-brms_fit_adjusted) |> 
  mutate(
    fixed_effects = map(fixed_effects, ~ .x |> mutate(significant = (Q2.5 * Q97.5) > 0)  ),
    random_plus_fixed_effect_of_plasma_proportion = map(random_plus_fixed_effect_of_plasma_proportion, ~ .x |> mutate(significant = (Q2.5.prop_plasma_observed_arcsin_sqrt * Q97.5.prop_plasma_observed_arcsin_sqrt) > 0)  )
    
  ) |> 
  saveRDS(glue::glue("{result_directory}/summary_epithelial_prop_plasma_DE_brms_digestive_tract_observed_proportions_all_tissues.rds"))

system(glue("~/bin/rclone copy {result_directory}/summary_epithelial_prop_plasma_DE_brms_digestive_tract_observed_proportions_all_tissues.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/"))


# Adjusted counts
pseudobulk_sample = tar_read(pseudobulk_sample, store = glue::glue("{result_directory}/_targets"))

m = 
  tar_read(
    adjusted_assay,
    store = glue::glue("{result_directory}/_targets")
  )  

pseudobulk_sample = pseudobulk_sample[rownames(m),, drop=FALSE ] 
assay(pseudobulk_sample, "counts_adjusted_plasma") = m

pseudobulk_sample |> saveRDS(glue::glue("{result_directory}/pseudobulk_epithelial_adjusted_digestive_tract_observed_proportions_all_tissues.rds"))

system(glue("~/bin/rclone copy {result_directory}/pseudobulk_epithelial_adjusted_digestive_tract_observed_proportions_all_tissues.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/"))


tar_read(  sce,  store = glue::glue("{result_directory}/_targets")) |> 
  # select(.cell, sex, age_days, dataset_id, observation_joinid , sample_id,  contains("cell"), disease, collection_id, title, contains("ethnicity"), assay, tissue, tissue_groups) |> 
  # dplyr::select(-run_from_cell_id) |> 
  # edit_covariates() |> 
  zellkonverter::writeH5AD(glue::glue("{result_directory}/single_cell_epithelial_adjusted_digestive_tract_observed_proportions_all_tissues.h5ad"),  verbose = TRUE)
system(glue("~/bin/rclone copy {result_directory}/single_cell_epithelial_adjusted_digestive_tract_observed_proportions_all_tissues.h5ad box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/"))






# rownames(seu) = rownames(seu) |> mapIds(org.Hs.eg.db, keys=_, column="SYMBOL", keytype="ENSEMBL", multiVals="first") 
# seu = seu[!rownames(seu) |> is.na(),]

job::job({
  
  library(scuttle)
  library(BiocParallel)
  
  cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
  bp <- MulticoreParam(workers = cores , progressbar = TRUE)  # Adjust the number of workers as needed
  
  pseudobulk =
    aggregateAcrossCells(
      sce,
      colData(sce)[,c("sample_id", "cell_type")],
      BPPARAM = bp
    )
})

pseudobulk |> assay() = pseudobulk |> assay() |> Matrix(sparse = TRUE)
pseudobulk |> assay() |> colnames() = pseudobulk |> colnames()
pseudobulk |> 
  saveRDS("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_for_dharmesh.rds", compress = "xz")
system("~/bin/rclone copy ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/pseudobulk_for_dharmesh.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")



seu = sce
seu = seu |> as.Seurat(data = NULL)




seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

seu_list = seu |> 
  nest(data = -sample_id) |> 
  filter(map_int(data, ~ .x |> ncol()) > 30) |>  
  mutate(data = map(data, SCTransform, verbose = FALSE, assay = "originalexp", .progress = TRUE))


HarmonyIntegration

var.features <- SelectIntegrationFeatures(object.list = seu_list |> pull(data), nfeatures = 2000)

seu = seu_list |> unnest(data) |> PrepSCTFindMarkers()

seu |> VariableFeatures() = var.features

seu <- seu |> RunPCA(verbose = FALSE, assay = "SCT")

seu |> DimPlot(group.by = "sample_id")

seu = seu |> ScaleData(assay = "originalexp") |> FindVariableFeatures(assay = "originalexp") |>  RunPCA(verbose = FALSE, assay = "originalexp")

seu = 
  seu |> 
  RunHarmony(  
    c("sample_id", "dataset_id"), 
    plot_convergence = TRUE
  )

seu <- seu |> RunUMAP(reduction = "harmony", dims = 1:17)

plot_age = seu |> mutate(stage = age_days |> divide_by(365) |> divide_by(10) |> round()) |>  DimPlot(group.by = "cell_type", split.by = c( "stage")) + ggtitle("Decades") + guides(color = "none")
plot_sex = seu |>  DimPlot(group.by = "cell_type", split.by = c( "sex")) + ggtitle("Sex")

library(patchwork)
plot_age | plot_sex


a =AnnotationHub()

infile = read.csv("data/gene_lists.csv")
data = infile[,"ENSEMBL"]

library(org.Hs.eg.db)



