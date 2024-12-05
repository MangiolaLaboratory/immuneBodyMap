

library(tidyverse)
library(sccomp)
library(magrittr)
library(glue)
library(forcats)
library(stringr)

library(arrow)
library(dplyr)
library(duckdb)


# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/dharmesh_shared_mix/drop_samples.csv /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/")




write_parquet_to_parquet = function(data_tbl, output_parquet, compression = "gzip") {
  
  # Establish connection to DuckDB in-memory database
  con_write <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # Register `data_tbl` within the DuckDB connection (this doesn't load it into memory)
  duckdb::duckdb_register(con_write, "data_tbl_view", data_tbl)
  
  # Use DuckDB's COPY command to write `data_tbl` directly to Parquet with compression
  copy_query <- paste0("
  COPY data_tbl_view TO '", output_parquet, "' (FORMAT PARQUET, COMPRESSION '", compression, "');
  ")
  
  # Execute the COPY command
  dbExecute(con_write, copy_query)
  
  # Unregister the temporary view
  duckdb::duckdb_unregister(con_write, "data_tbl_view")
  
  # Disconnect from the database
  dbDisconnect(con_write, shutdown = TRUE)
}

# dplyr::tbl(
#   DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   dbplyr::sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_1.parquet')")
# ) |> 
#   inner_join(input_relative |> select(-n), copy=TRUE) |> 
#   write_parquet_to_parquet("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/cell_metadata_1_0_1_sccomp_input.parquet")

# system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/cell_metadata_1_0_1_sccomp_input.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/caq_celltype_level_map.csv /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/")



library(targets)

tar_script({
  
  #-----------------------#
  # Input
  #-----------------------#
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(qs)
  library(crew)
  library(crew.cluster)
  
  #-----------------------#
  # Packages
  #-----------------------#
  tar_option_set(
    
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    #  debug = "dataset_id_sce_ce393fc1e85f2cbc", 
    # cue = tar_cue(mode = "never"), 
    workspace_on_error = TRUE,
    format = "qs",
    
    #-----------------------#
    # SLURM
    #-----------------------#
    controller = crew_controller_group(
      
      
      
      crew_controller_slurm(
        name = "slurm_1_80",
        script_lines = "#SBATCH --mem 80G",
        slurm_cpus_per_task = 30,
        workers = 20,
        tasks_max = 1,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_200",
        script_lines = "#SBATCH --mem 200G",
        slurm_cpus_per_task = 2,
        workers = 20,
        tasks_max = 1,
        verbose = T, 
        seconds_idle = 30
      )
    ),
     debug = "estimates",
    
    resources = tar_resources(crew = tar_resources_crew("slurm_1_80")) 
    #, # Set the target you want to debug.
    #
  )
  
  #-----------------------#
  # FUNCTIONS
  #-----------------------#
  
  #' Age Bin Assignment Based on Sex-Specific Life Stages
  #'
  #' This function assigns individuals to age bins (life stages) based on their
  #' age in days and biological sex. It provides sex-specific age categories
  #' reflecting typical immune system ageing phases in males and females. For
  #' unknown sex, it averages the age thresholds of the male and female bins.
  #'
  #' @param age_days Numeric vector. Ages in days for each individual.
  #' @param sex Character vector. Sex for each individual, with values
  #'   "male", "female", or "unknown". Must be of the same length as `age_days`.
  #'
  #' @return A character vector containing the corresponding life stage for each individual.
  #'   The life stages are defined as:
  #'   - For males: "Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior".
  #'   - For females: "Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior".
  #'   - For unknown sex: The function averages male and female thresholds to assign an appropriate life stage.
  #'
  #' @details
  #' The function determines life stages based on typical developmental, hormonal,
  #' and immune system changes at specific ages, with slight variations for males
  #' and females. For individuals with an unknown sex, the function assigns a life
  #' stage based on the average of male and female age thresholds.
  #'
  #' @examples
  #' age_days <- c(500, 4000, 15000, 25000, 30000)
  #' sex <- c("male", "female", "male", "female", "unknown")
  #' age_bins <- age_bin_sex_specific(age_days, sex)
  #' print(age_bins)
  #'
  #' @importFrom dplyr case_when
  #' 
  #' @export
  age_bin_sex_specific <- function(age_days, sex) {
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
  
  age_range_sex_specific <- function(age_bins, sex) {
    # Initialise an empty vector to store the age ranges in years
    age_ranges <- integer(length(age_bins))
    
    # Define age thresholds for male, female, and unknown
    male_thresholds <- c(3, 13, 21, 40, 55)
    female_thresholds <- c(3, 13, 19, 36, 50)
    unknown_thresholds <- c(3, 13, 20, 38, 52)
    
    # Map age groups to approximate age ranges for each sex category
    for (i in seq_along(age_bins)) {
      if (sex[i] == "male") {
        age_ranges[i] <- dplyr::case_when(
          age_bins[i] == "Infancy" ~ male_thresholds[1] / 2,
          age_bins[i] == "Childhood" ~ (male_thresholds[1] + male_thresholds[2]) / 2,
          age_bins[i] == "Adolescence" ~ (male_thresholds[2] + male_thresholds[3]) / 2,
          age_bins[i] == "Young Adulthood" ~ (male_thresholds[3] + male_thresholds[4]) / 2,
          age_bins[i] == "Middle Age" ~ (male_thresholds[4] + 55) / 2,
          age_bins[i] == "Senior" ~ 55,
          TRUE ~ NA_integer_
        )
      } else if (sex[i] == "female") {
        age_ranges[i] <- dplyr::case_when(
          age_bins[i] == "Infancy" ~ female_thresholds[1] / 2,
          age_bins[i] == "Childhood" ~ (female_thresholds[1] + female_thresholds[2]) / 2,
          age_bins[i] == "Adolescence" ~ (female_thresholds[2] + female_thresholds[3]) / 2,
          age_bins[i] == "Young Adulthood" ~ (female_thresholds[3] + female_thresholds[4]) / 2,
          age_bins[i] == "Middle Age" ~ (female_thresholds[4] + 50) / 2,
          age_bins[i] == "Senior" ~ 50,
          TRUE ~ NA_integer_
        )
      } else if (sex[i] == "unknown") {
        age_ranges[i] <- dplyr::case_when(
          age_bins[i] == "Infancy" ~ unknown_thresholds[1] / 2,
          age_bins[i] == "Childhood" ~ (unknown_thresholds[1] + unknown_thresholds[2]) / 2,
          age_bins[i] == "Adolescence" ~ (unknown_thresholds[2] + unknown_thresholds[3]) / 2,
          age_bins[i] == "Young Adulthood" ~ (unknown_thresholds[3] + unknown_thresholds[4]) / 2,
          age_bins[i] == "Middle Age" ~ (unknown_thresholds[4] + 52) / 2,
          age_bins[i] == "Senior" ~ 52,
          TRUE ~ NA_integer_
        )
      } else {
        stop("Each element of 'sex' must be either 'male', 'female', or 'unknown'.")
      }
    }
    
    return(age_ranges)
  }
  
  create_input_cell_counts = function(cellNexus_metadata, drop_sample_df, result_directory){
    
    
    ethnicity_grouped <- tribble(
      ~self_reported_ethnicity, ~ethnicity_groups,
      "unknown", "Other/Unknown",
      "European", "European",
      "Korean", "East Asian",
      "Asian", "East Asian",
      "Japanese", "East Asian",
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
    
    system(glue("~/bin/rclone copy box_adelaide:/minh_immune_map_disease/disease_data_grouped_further.csv {result_directory}/"))
    
    disease_data_grouped = 
      disease_data_grouped |> 
      select(-disease_groups) |> 
      left_join(read_csv(glue("{result_directory}/disease_data_grouped_further.csv"))) |> 
      mutate(disease_groups = if_else(disease_groups |> is.na(), "other", disease_groups))

      tbl(
        dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        sql(glue("SELECT * FROM read_parquet('{cellNexus_metadata}')"))
      ) |>
      
      # Filter empty droplets
      filter(!empty_droplet) |> 
      
      # SEX edit
      mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |> 
      
      # ETHNICITY
      left_join(ethnicity_grouped, copy=TRUE) |> 
      
      # TISSUE
      filter(!tissue_groups %in% c(
        "muscular system (skeletal muscles)",
        "ovary",
        "vasculature",
        "digestive tract junctions and connections",
        "peritoneal and abdominal cavity structures",
        "connective tissue",
        "miscellaneous glands"
      ), 
      !tissue_groups |> is.na()
      ) |> 
      
      # TECH
      left_join(assay_data_grouped, copy=TRUE) |> 
      
      # DISEASE
      left_join(disease_data_grouped, copy=TRUE) |> 
      
      # IMMUNE CELLS
      filter(is_immune) |> 
      filter(cell_type_unified_ensemble %in% c("cd8 naive", "cd16 mono", "cd4 tcm", "cd4 th17 em", "granulocyte", "cd4 th1/th17 em", "treg", "erythrocyte", "b memory", "b naive", "nk", "plasma", "cd4 th2 em", "mast", "cd4 th1 em", "cd8 tem", "mait", "tgd", "cdc", "cd4 fh em", "cd4 naive", "nkt", "macrophage", "cytotoxic", "cd8 tcm", "cd14 mono", "pdc", "ilc")) |> 
      
      # TEMPORARY. de-group pancreas and liver
      mutate(tissue_groups = case_when(
        
        tissue %in% c("gallbladder") ~ "gallbladder",
        tissue %in% c("pancreas", "exocrine pancreas") ~ "pancreas",
        tissue %in% c("liver", "caudate lobe of liver", "hepatic cecum" ) ~ "liver",
        TRUE ~ tissue_groups
      )) |> 
        
      dplyr::count(sample_id, donor_id, dataset_id, title, collection_id, age_days, sex, ethnicity_groups, tissue_groups, assay_groups, cell_type_unified_ensemble, disease_groups) |> 
      mutate(n = as.integer(n)) |> 
      as_tibble() |> 
      
      # Age
      filter(age_days > 365) |> 
      mutate(age_bin_sex_specific = age_bin_sex_specific(age_days, sex)) |> 
      
      # Here we drop those samples with a low cell type entropy. E.g. one cell type only.
      anti_join(drop_sample_df) |> 
      

    
    # Set intercept
      mutate(
        ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
        assay_groups = fct_relevel(assay_groups, "10x Genomics 3"),
        disease_groups = fct_relevel(disease_groups, "Normal"),
        age_bin_sex_specific = fct_relevel(age_bin_sex_specific, "Adolescence")
      ) |>
      
      # Center based on adolescence
      mutate(age_days_scaled = age_days  |> scale(center = 15*265) |> as.numeric()) |> 
      
      # Add hierarchy of cell types
      left_join(
        readr::read_csv("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/caq_celltype_level_map.csv"), 
        by = join_by(cell_type_unified_ensemble == cell_type_unified_harmonised)
      )
      
    
  }
  
  #-----------------------#
  # Pipeline
  #-----------------------#
  list(
    tar_target(
      result_directory,
      "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised", 
      deployment = "main"
    ),
    tar_target(
      drop_samples,
      {
        # evaluate result_directory for targets
        print(result_directory)
        
        system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/dharmesh_shared_mix/drop_samples.csv {result_directory}/"))
        
        read_csv(glue("{result_directory}/drop_samples.csv"))
        
      }, packages = c("glue", "readr")
    ),
    tar_target(
      input_relative, 
      create_input_cell_counts(
        "/vast/projects/cellxgene_curated/cellNexus/cell_metadata_cell_type_consensus_v1_0_1.parquet",
        drop_samples,
        result_directory
      ),
      packages = c(
        "dplyr",      # Data manipulation (mutate, filter, count, left_join)
        "tibble",     # Creating tibbles (tribble function)
        "duckdb",     # For connecting to DuckDB and reading parquet files
        "glue",       # For constructing SQL queries with glue syntax
        "readr",     # For reading CSV files
        "forcats"
      )
    ),
    tar_target(
      saved_input_relative,
      {
        # evaluate result_directory for targets
        print(result_directory)
        
        file_name = glue("{result_directory}/cell_metadata_1_0_1_sccomp_input_counts.rds")
        input_relative |> 
          saveRDS(file_name)
        
        system(glue("~/bin/rclone copy {file_name} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/"))
        
      }, 
      packages = "glue"
    ),
    tar_target(
      formula_df,
      tribble(
        ~ formula_composition, ~ formula_variability, ~ name,
        
        # continuous age
        "~ 1 + age_days_scaled*sex + disease_groups + ethnicity_groups + assay_groups + 
          (1 | dataset_id) + 
          (1 + age_days_scaled*sex + ethnicity_groups | tissue_groups)",  
        "~ age_days_scaled*sex + disease_groups",
        "estimates_continuous_age",
        
        # discrete
        "~ 1 + age_bin_sex_specific + disease_groups + sex + age_bin_sex_specific:sex + ethnicity_groups + assay_groups + 
          (1 | dataset_id) + 
          (1 + age_bin_sex_specific + sex + age_bin_sex_specific:sex + ethnicity_groups | tissue_groups)",  
        "~ age_bin_sex_specific + disease_groups",
        "estimates_age_bins", 
        
        # continuous + discrete
        "~ 1 + age_days_scaled + disease_groups + age_bin_sex_specific*sex + ethnicity_groups + assay_groups + 
          (1 | dataset_id) + 
          (1 + age_days_scaled + age_bin_sex_specific*sex + ethnicity_groups | tissue_groups)", 
        "~ age_days_scaled + disease_groups",
        "estimates_continuous_age_plus_age_bins",
        
        # disease tissue specific
        "~ 1 + age_bin_sex_specific + disease_groups + sex + age_bin_sex_specific:sex + ethnicity_groups + assay_groups + 
          (1 | dataset_id) + 
          (1 + age_bin_sex_specific + disease_groups + sex + age_bin_sex_specific:sex + ethnicity_groups | tissue_groups)",  
        "~ age_bin_sex_specific + disease_groups",
        "estimates_age_bins_disease", 
      ) |> 
        expand_grid(cell_type_level = glue("L{0:3}") |> as.character()) |> 
        group_by(name, cell_type_level) |> 
        tar_group(), 
      iteration = "group",
      packages = c("tibble", "glue", "targets", "dplyr", "tidyr")
    ),
    tar_target(
      estimates,
      input_relative |> 
        
        # I have to summarise further because I have counts already
        with_groups(
          c(sample_id,age_days_scaled,age_bin_sex_specific, sex, disease_groups,ethnicity_groups, assay_groups, dataset_id, tissue_groups, all_of(formula_df$cell_type_level)),
          ~ .x |> summarise(n = sum(n))
        ) |> 
        
        sccomp_estimate(
          formula_composition = formula_df$formula_composition |> as.formula(),
          formula_variability = formula_df$formula_variability |> as.formula(),        # Differential variability
          sample_column = "sample_id", 
          cell_group_column = formula_df$cell_type_level, # A level of the hierarchy
          abundance_column = "n",
          cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)),
          mcmc_seed = 42,
          verbose = T, 
          bimodal_mean_variability_association = TRUE,
          prior_overdispersion_mean_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40)),
          output_directory = "/vast/scratch/users/mangiola.s/my_draws", 
          max_sampling_iterations = 5000, 
          inference_method =   "hmc", 
          refresh = 1
        ),
      pattern = map(formula_df),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_80")), 
      packages = "sccomp"
    ),
    tar_target(
      saved_and_tranferred,
      {
        # evaluate result_directory for targets
        print(result_directory)
        
        local_file_name = glue("{result_directory}/{formula_df$name}___{formula_df$cell_type_level}.rds")
        local_file_name_FIT_FOR_PORTABILITY = glue("{result_directory}/{formula_df$name}___{formula_df$cell_type_level}_FIT_FOR_PORTABILITY.rds")
      
        # Save draws as monolythic
        attr(estimates, "fit")$save_object(file = local_file_name_FIT_FOR_PORTABILITY) 
        system(glue("~/bin/rclone copy {local_file_name_FIT_FOR_PORTABILITY} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_estimates_1_0_1/"))
        estimates |> attr("fit") = readRDS(local_file_name_FIT_FOR_PORTABILITY)

        # Save sccomp estimates
        estimates |>  sccomp_test() |> saveRDS(local_file_name)
        system(glue("~/bin/rclone copy {local_file_name} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_estimates_1_0_1/"))
        
      }, 
      pattern = map(formula_df, estimates), 
      resources = tar_resources(crew = tar_resources_crew("slurm_1_200")), 
      packages = c("glue", "sccomp", "magrittr")
    )
    
  )
}, 
ask = FALSE, 
script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets.R"
)


job::job({
  
  tar_make(
    script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets.R", 
    store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets", 
    reporter = "verbose" #, callr_function = NULL
  )
  
})

x = tar_read(input_relative, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/_targets")


tar_workspace(estimates_2e499e82b1dbf3a0, 
              script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets.R", 
              store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets"
              )

tar_meta(starts_with("estimates_"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets")

# Age proportion prediciton
estimates_age_bins |> 
  sccomp_predict(
    formula_composition = ~ 1 + age_bin_sex_specific*sex + (1 + age_bin_sex_specific*sex | tissue_groups), 
    number_of_draws = 100,
    summary_instead_of_draws = TRUE
  ) |> 
  mutate(age_bin_sex_specific = factor(
    age_bin_sex_specific, 
    c("Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior"),
    ordered = TRUE
  )) |> 
  mutate(age_bin_sex_specific_numeric = age_bin_sex_specific |> as.integer())  |> 
  saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/prediction_age_bins.rds")

system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/prediction_age_bins.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

tar_read(formula_df, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets")


# For Hong
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins.rds")
estimate_age_bins = estimate_age_bins |> dplyr::select(-count_data)
attr(estimate_age_bins, "fit") = NULL
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_age_bins_effect_tibble_only.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_age_bins_effect_tibble_only.rds box_adelaide:/immune_map_disease/")

# Save fit
library(magrittr)
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins.rds")
estimate_age_bins |> attr("fit") %$% save_object(file = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins_FIT_FOR_PORTABILITY.rds") 
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins_FIT_FOR_PORTABILITY.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
estimate_age_bins |> attr("fit") = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins_FIT_FOR_PORTABILITY.rds")
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_age_bins.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")


# estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_age_bins.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/21_11_2024_sccomp_archive_before_factor_ordering/estimates_age_bins.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

# Benchmark
tic()
estimate_age_bins |> sccomp_test(contrasts = c(  "respiratory system" = "sexmale + `sexmale___respiratory system`",
                                                 "blood" = "sexmale + `sexmale___blood`"))
toc()


