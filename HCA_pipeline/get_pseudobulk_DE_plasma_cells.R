
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

library(readr)
library(forcats)
library(glue)
edit_covariates = function(tbl, disease_tbl){
  
  
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
    
    dplyr::select(cell_id, sample_id, donor_id, dataset_id, file_id_cellNexus_single_cell, title, collection_id, age_days, age_bin, sex, ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble, cell_type, disease_groups) |> 
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

# result_directory = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1") 
system(glue("~/bin/rclone copy box_adelaide:/minh_immune_map_disease/disease_data_grouped_further.csv ./"))


# REMOVE OLD CACHE which is here get_default_cache_dir()
get_metadata() |> 
  edit_covariates(
    read_csv(glue("./disease_data_grouped_further.csv"))
    ) |>
     
  # filter( your filtering ) |>
  get_single_cell_experiment(atlas = "cellxgene")








sce = 
  # get_metadata() |> 
  
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata.1.0.5.parquet')")
  ) |> 

  filter(!file_id_cellNexus_single_cell |> is.na()) |> 
  filter(sample_id %in% c("324c08b06f9cc6e5f94c519861fece18", "d64e92e8c10522dc689b5eb2b87b00eb", "142a094ece1354293725ceb543fb831c", "f1f3f44f8a922b6bffeb3d1410ebefbe", "88ac6575907fa32b69eb5e1fba7944dd", "236bed4ef2516132fd27f0cc3da1e830", "815b4b736442e536213fdfaeba1c897a", "edd7ca8e111667dcdfc2d6ca49d91c63", "41b515ed2abc918fe3d83833068dc26e", "08e5b3d26728edd70c080e64f8bb4b0f", "f8699441d4ee0ba86560d0670a4d7f67", "8adeffa859e1938c8bd8771f98fd7b67", "c0e7a651d24619e3eecbf4f52055e0bc", "69fb03ba5c5b49b390669f459c8b1f5a", "007c3b188d1350067d50a4f17dd13aad", "d45b5d329e336b8864b29f6ba0ce4062", "251782cdf598068c0690ef9cd368c88f", "b1bd7bdff1e3d362ae4f2611e20f34cf", "60daa546af935991e11c176e8af3bb42", "1483aa5f3081631232668b2955b584c9", "365623a317dc279711f1749b6523e8f0", "437865f81c572118e382b12d74058e2c", "2b6d4dcc1c00988b9390cebd2c4939b7", "3b641516862bc2d1926b238fcd9ee4db___-1-HCATisStab7646029", "3b641516862bc2d1926b238fcd9ee4db___-1-HCATisStab7646030", "ac4969702dafceb269516b8c8c22c4d2", "0e7f5caa5854e50db93f7d49b3f155e6", "89e09195652f5954352b7472eb18fb26", "0a226bdefe11b725e6d761fd63ffcfbf", "c1b66d20e7d09f804b746fb5e0d2c895", "15086e9c78d3b6576e74a01d7e374721", "1644b16a62b0869717c3a590b9107810", "00bb4075752d2b2400de7cbcb177f51f", "8e1bc649ac42ba0ebfaeefb907a7335c", "4ce0a77ddae58f08c30987622e0010fc", "76cebfe97bba892b70034d2bb0407bb8", "33a7e5937e0b507546ec05cb3b197db6", "eb1fdac90e267aee4c056e0c600715fc", "35eea6755a192a7fb4c1247034bc3554", "152ccd2751fdbc4d507487c833f235d5", "a2d5028cf1f74a1165f5191d079fc45c", "196d178905a7c133ea11ea3c31c96d12", "25cfd68ff343955c7cf5c6d5c595fb24", "9ea072fe5cbd55a60feeb34024743aeb", "6d55502923b1be628a3411b005d1f212", "f64e09da25670f6ae9472e4ea6a19cc1___FCAImmP7352192-", "5ad20fced093387f483002cf7c4f1bef", "191e41d014d381a4c8b93ca2ceac9b10", "207c268915d4e7c5b8f5025a36ad680d", "850237d5f345b8c0833897a102962324", "66a15d552c5b4d43969d36c6632b9031", "d17a8ed86d8179b11f57424abe12447c", "c14a7613368276ba7004e29267312af9", "4aa38bebfeb6e1ebbdac21e8fe8cddb8", "fdf583ed81a6d0d0818aa3955191c60c", "f62031d662fd9a9861c0af8b784c7b72", "1ccd1c6b46022a0a7ef3dbe6dbaeb58e", "a72d5f0c962f7be89a3ecc4129f39fb2", "b5a691c284259293479c4a676bf69850", "c1276a0391770a26ff538d6af727ebd4", "7ee41772f59c60c972cad75673594aca___FCAImmP7862086-", "c017a423dd8ab5ca391de0f5da42f1fd", "207249aecaa55d00520026e59d02b598___FCAImmP7316893-", 
                          "f08e2be183fd75508b136b08a2924abc", "b8a99fc28db99109853fd4b6561190e6", "d87999e384720e16e6ff174dfc989b04", "e9f60341dfd36966a7372b319d17a30a", "9fb13b7f2af733bce02fe5eb92726da6___Liver_1.-1", "20a985795fedae625d91b3fbde4daa43", "3f3741d19b5ed262e13f234b0cba10cf", "cf50bd673c6b4be1c288bdd7b172df91", "7ea930d6bbe694419468eb9b96b7b698", "244a932a241bbcd2746b2cbb1916d06f", "aad0f57f4d2042b7cc809daba757a077", "f05e98ec0cb20bd5635caa79a544d8fe", "d9ac1987dacebd3876b50a13986374dd", "2a69f46dd2d306928ffe00893bb396b2", "ffd6706b0f19c066c30abe9fdd96ac60", "eadf351fa9fd61c0427342910e9a07b5", "5a214ea61c1ded896a05d3c76f529cf8", "e9132534f38d04681a1ae40faecf2f91", "4fa8e1d4c43011ffb37932c5ccace463", "2ddd147dac6deb0de640916f122f0927", "8dc3539ccecf875badffc8c562bdeb5f", "b9f478ed495cb7699fab8022f9922800", "3754e2a89e94b36ed78b0e3cffef122b", "07483d517f083cafcd70fe0e0e6c0107", "c3097e6e796076596ece1f2834a8b0cc", "c36c7215d1e083dac16ef93b33ac83cb", "a7fcaa665ab2706230c34cb622b6973a", "ea5bdbb220e6643f99d3ea5df4b080d3", "e4e6d6e1310ce4c283e80a20cf0a1acd", "b0a1999b5c66c8c27511fa03cc9ffd16", "c92a835eba5abe879b7028b62dd5189f", "ee578407f43242052fd3330e5f6e97e2", "2fb3bfad870b07b0abdaf599db0606a2", "c89736fe42238dac8e6b7709ff64b056", "a7ce5e530a45cc9d2cdb43da1cc5bc68", "8c9cffb50823a4b69acef91f7385acd9", "6578b5fb5a7257a30295eb629147c47e", "9d7f2e1c9a674cb339e44ca3905c917a", "522dae9725628ccc31b1c7517631b55b", "c22c9e1c36c37dba26d18ad33313a00f", "1c426357e9f440e553ffeaf9a66fd024", "ddf64bd854b3b637e3c2b0ca7fc71cb7", "f995cb92736f525d78b90d6ad69c4c3e___-1-HCATisStab7413622", "6be02231c40426b4139ca995b8a7882d", "5365a79e91bc8c59a35f5cdd904e9c96", "30d6f0a9b6f9d36d64c02aae25b92f1c", "18262b4b56c2ef22f4ee59d99c48d7dc", "12270d98b43c1c3ebac904257bceb9fc", "af94226061547db2b89f35e1a2b56dcd", "6e73b00d89f1ed1b3a556886a5ef6247", "a712a60b8bed0b684843b29a1c576d74", "517a30378a36e3cd2a64faded645225b", "99bff9753af35adb67f7495e8e4906fe", "0d80ed818e5cf5e014150922082df42f", "e04699bfa84d14c17e4550a151d5d1dc", "a4fdeda21ca66405c834c87546344e7e", "87482205199188c8b3624905e4b7bdd1", "d6f8a1e68abb75f1719685e73c2357b7", "72aff68d180c800c240a7e81f50858c3", "4b78e0159490c0d53b48a3994fee6624", "dc3d8b09a82b221bd401cdc74b7a43fe", "485f7847d8e44a4af8b7526bd6dc4b8c", "9146172bf36250f5ecf97f41d460a64d", "ff3f3215060b48696327c3dfe31b3a82", "14177f7850fecaef85fd6230e4c37bfe", "d7505e204bb90f0bf5aea8c86611ab01", 
                          "8c77ea423570dd448783211dc3aea653", "287a8b9da08ef5e00d7f1383eefe02f9", "411b0c2ad491a221c8025d0f933e2ff3", "a9e6b4485e689a947ab7cd878401b467", "4ebabb958c9f9c130df27d1da4edd0a4", "be626f20511df9bec5853c2f6fadf25e", "b4e27e3cc8ff8f3c03767225c6855f24", "38896afe3e2dbb90c4435aca4bb432fb", "deec7acb441cc0a9be537dea66635548", "a1e74af7c3d180400e96eef055d3ef96", "06a68b03b19d658f4b1975db23302081", "d4740d54b614098875712c25d997a8da", "44bc5416300d8745d5582af1592b7e79", "3e8f961a0bea4d6335ae12f76c5614dd", "22878e76754f89d2843585cd3e762008", "d4c0540de0007881d604a18556226bde", "34a8de0c539aeec51db3f2bc0cc5e160", "ad0947e467d0d718fb878742b1d78106", "ad72137b58a7fd9f8076022011cc4b2e", "8367c0051eb7f0bd4258f6634ee6531a", "1d56c6bcbb6f0cbc3e53f2365937b901", "9fb13b7f2af733bce02fe5eb92726da6___Liver_2.-1", "1e8759e4cc269b15f4076c77e3e2e18b", "281e3b5f19bcfc218f328e3ccc170d94", "11ed03fc1be9d3e3901d87129badc09a", "5b2827a3f88aefca7db4b32f06887c3d", 
                          "71fe8a45931887ef7c4044bf77bdd0f2", "1554596e035f9027550664f2f85f858c", "e07cf75e32a99afad9cf0b8711500e20", "9d951a4a20a20fcdb3943f2c591fda62", "8d65021b8890dbed51a3169ed6b66098", "12a2b06684d7bb7e9f5514b987180bc2", "2199875e6ec6ab3b2c775b6f25af72ba", "89ce6e6657707f73a6c2a5d0f73cb1bc", "53ecc73b6607077664d8e06716597702", "3d9599dada75108be69ce324c4c61f5c", "4462cc165eb3213a954bff284d263be6", "77a85321af9a2e9b26a2f7a8b832f2a2", "22fdaf0127511a76a74a2600a0892fb3", "1f5d7784b3deab3259690d59fa884982", "c186470cbf92d7b9be5c76470fcf3c28", "4978c7148f7809d01f038c27e1ea88e8", "1d5efe6cd4c3cc11c77c07e3e55ae908", "71d18b8276f3d4ad719710732a13ce14", "37a443d9158d75bf561cab3630931311", "0196ea6efea73f75461772670e2a94f4", "f2880c533a82828854c447e957e89974", "f499388d2b7d3f0ef972128979c16416", "ffb53d83fe9b540c758eddefdac928ca", "33dc9dc87e9f5d06b51563ceaa655bb9", "057d3011e2cc7573bd822c08c73a65ba", "2e4a6c83a3a9f809063a735c4de8da8e", "3b641516862bc2d1926b238fcd9ee4db___-1-HCATisStab7646031", "4be89e66d8479a43b30dd06e7cb386e8", "4690c1f8b45c64ab0f0b8c730382001c", "7542e4139772e4c348a55e23831393ed", "eec28f1c966b73c7c657c805e671bce4", "9db19e06abe119b72b4e7894b60cb31a", "56f151e5187ab1632a9c966da2600318", "467c99368ba6f147d6ea07dca3271d42", "c9e1e8fb4e67005b0886503a24d1f6af", "955e14dd72e096ce4de90aed7d9c901c", "cfef7899bfbc408d587601ec933a627b", "58d906cf7b7237250b92b15aed9ddab4", "395ce7190186c6d91c978b12e5e42dde", "e23dbb64d92e8d7deb7db2f145ae4932", "30a33e7eee245f84f2cbaec09c367774", "85a8d454188f29a3f5cd72cb76068578", "6fdfc72d95f5b4cd2d7dc15726946e4d", "50f9e82584cec7b05298ab39aa580e04", "c5ea775eaa01b3a6dfce646e4f4d25f6", "2841828d6cb4993db9e01a06bef7fd5e", "6e44fb534e997f1711b55936206a1cc6", "7ee41772f59c60c972cad75673594aca___FCAImmP7862085-", "7ad8848f3c4348b9e18ab6c8dc0e9a2a", "7d1c69cd68c6c1e8696db3195d52d6de", "8a55c92c3d79ef2bf7cf2f8d078fc515", "e96366f1d20e12c781b91947de5acd7b", "e4958abc555a7d10cc1a2919103f8ad8", "c9fc52515d26dc52d54de55359430737", "26a94c8572a37300c023b1ab35b0f04b", "8c5fc0b6896d0358e137ccbb42486028", "400d996bbf6aa503bfdf10649333eff5", "242586ad563959520c4b756db5bd448a", "45accbf8f171f27c1a67046aac58f428", "68a3d97da9939273826272bb1532511f", "59a8c618c06b1132fa0cf8707ac6fd85", "bbc95ff97c19ec8789e9d8c9cbdce40f", "4acb7d9a84b139bc62e5d5e706a36912", "f9999504ab3d8eeb1820b55b51068b86", "3ecc6621e202ef68fb2d1b43a3c3c214", 
                          "9aa1dcc4d24abbe237ccd01f9f7dc664", "2bb5e5b8039ae8b6c0335b21afa05840", "240ea57948d0f78736c327d207d14ca1", "f1094931956f4f696d09dd15f9330458", "16af564bc3cbc7ae6ac41cebefa08f2e", "f5601be8089fbdb0d346d07b3e637f6e", "b1f5e8fc34250b4b04eded62d59dbffe", "2a0473c3ce8e254cae5cb78a34f23ff8", "07c9502b0dc6946d095279fdf8ec814a", "86663d28fddfd350999c1d3410b27bb2", "a37b484d6a19e8b970e9be52d58bdac0", "6b1e92d20d31537f331e6c09c318d858", "5f3079d4303baf92547b6ec523cc44ab", "d249cb364eb422f2f28972cdb18f40b7", "0dfee248893a2ebb1d666f83d055cd74", "1a9746c7100232c1b4269c8fb98d530a", "c31d9cea4c2c615bfd6a1bbf5d809a29", "8833fe83aa513214d50d7bda69eaba03", "66e61fb40ab6c3faed02d25a9a421ac8", "90724f97073f3c0386d5beb66af097cf", "a4ed7eadb99c8d479f14db4d6e2045ac", "41d33e29a4b149d69a2d7c6206940f56", "4df0f90e7724939d6217336cf0e4a7bf", "0edfde5ae41acf53ec3f6ef7c2d49533", "2b921bb8d9f281cd4dbf73ebab1e7bc9", "7979a97107ee32eac6f812f03fd81657", "ebc46ad59dd2f712576ad11cc76cd65d", 
                          "207249aecaa55d00520026e59d02b598___FCAImmP7316892-", "50c242f6eac91db7e1be1a56dc9852ff", "51609f92cd2bb35cd7e7158c577046b6", "74bcf70ecca8af384c31564979f60d82", "fad6c855286af00fb0ad93bd7745f527", "f47f6d835ccfe9a5cb23423ad3bc3fbd", "892fc55f4621f8bda01b503e20a30605", "d68d07a31d49a664beb6a679c640376d", "c1beb872af9ee00494073a0a056dad6e", "989cfd6f5e379fb1fdfeffb47bf803d7", "b430f5406fbe85c83da3ee759dcc0bab", "4a6b3922c90c5636bbcaa61360f4de09", "7764e0477ef9b92fd7c09ea0c95a0171", "6a1c26d726f6281e7f17c8ea8e00553e", "170529ab25b23c2e27b1d190716e0e83", "feb4683498605f15e374394c53f55bb3", "a02d977aefa796a105d8595ed670ac01", "8d3341375aa1cb978793aa5bf3139aec", "3679c770df59b3232939c5788ffc66ac", "b87d85adbbf785cdc98cb7879d7483e0", "912fddeba29809d3c667d281d0e0d63d", "4c91aa3bdbb6a4c8dea678da82371230", "f9181d0892bbcc5d37c1473e541621f1", "a0916ab00e452c722e640300c88b4f66", "e029658622d33a00eda3c226e4b7aa62", "2ffe51d91b27e9e577fe9f6352e1455d", "ce005a713a2bad0c9e391414631efc0d", "33ee187e57053d9849944fc8373671dc", "9eab31abf604d8ff747a350fa8af7fcc", "1a588c0a6412c1042017911613b81f8f", "f368270d223d2f82257ea09ee3dd15f0", "161a6b7e9410ce9969f1682d671fbbd2", "abbd8481b02bf4894d2105fe7557c25c", "134fb29105f3f6bc51643cadc2c03213", "21a8b792422599811a25493fb9ea7ff9", "cb02a149a7253d4f3d69b1f4f9e0aa10", "78ea954aa461ea4cf3b6386af276bda0", "c932b02196ab9ca9782a973fe6b2f337", "8c81ea8471795391148b86d285b4928a", "326978c1fdeb1420b302a598b33eda73", "9c9b9a2b78e64bd6f59307ae817e7de0", "146ffe54a63ae5095020585d9d52dc4a", "c073a45f44c78f3c5dfd448180f81ff6", "f97562c1f9c9f6404ceb2726d26431b7", "0456df48da29597769fb85a3e21b9461", "7a83e38f165a0aa67fd9a99a25604ecc", "4cbc3ea0d2ff2be1f7d055159567b235", "a812e4cf975baae5eac574137979b0a0", "34b4b7216741c252891938da1b320f13", "a8b97f87d6bf54f597102b7d36c0262a", "5b1acfe55017b92594b59ac83a2eb555", "96771e61f24a5eb3d9e3b5f7178b1991", "39bf9632c69f5492424c5a8cf129c11f", "2734aa3ec894827a9326fcc04f0b93d5", "207249aecaa55d00520026e59d02b598___FCAImmP7316891-", "f75bb05a80d098d856a511ab715849b2", "80f1433ad44a1a79451af8540ea390f7", "5dc01d67bb9c755f265cb8f406f756f2", "39c9537ad5808b4497c62ab30101db21", "d83a910e0cd0bc85c833989c0d401973", "f64e09da25670f6ae9472e4ea6a19cc1___FCAImmP7352194-", "3d6bcd0671edb4f3cb829b3343cb2c8c", "200c8b51df098470167fbd8ad08f8ed0", "e1a07bc0e58356da0c138f3cad06ac97", "645cf62b6e7e15b2794b0a6038eba85d", 
                          "cd40c40180b897edef958d29ac964382", "9750f2034f9dc6d4c2231a29837f6fd4", "778ec1bb30108f94034d1bba71a616dc", "2f9dc166931b43d20ae3273a9a48ada8", "1ca526d2299a00f7f64cc1ec3f8eeb98", "c480fb85956eb1af12796eaeb4d7a949", "79cfa8a9b89002b71192e1472212686e", "207249aecaa55d00520026e59d02b598___FCAImmP7316889-", "7b729821af1a763fa6406548cf450b70", "2165c9031f205dbea37cebf0fbc75d4e", "27188e1154e94a932b55c765857ae805", "27e314bfd5a4668d5447db17c0380a70", "b95f52a3cad9d28b7f89b2208eca0812", "85a4420efd22792cb6588eb30d868e15", "cd1826c3733b84b12e118598db8c9c25", "05a3005571a76ed7362a164aeb7bceac", "df8f06ef6ed77c8b3a95e6800ef9acd1", "eea99222eded7e734cf5522cb86d45bf", "595196e6479375ab0796b58a45881e38", "81e1873ecb463dbd0374b1d23111ca42", "4acd223baae5e0dc84f0ece1d68a3593", "efc105ee8462651fd81c94a3dd9d1ae4", "a4d7fd96d086858e8a1b1b550ff3031d", "70dbd5805a57d4f4f2f413e576b7580f", "26d81135bc216cb599d7e9a224f50cce", "1056ba6b906eb2d8e56c312f2d4a754c", "b74b808a3f8751fb00986e334a69d114", "0389dc91de3e38a6b575524e2bc78509", "0389d4773cc583a12f2291d33c431e44", "c63a15c02ac26f7a2c6d4e7f4f9c15e2", "e2bd36f51110b103c89ab94d3a7ba53e", "40e4196a43875989f600af6fcee42293", "c0853c48dcfe80118e9678a230cd6f6f", "158be29ada5c79db87314963f4ebc331", "4a783360f238a3bc12feaff5597d8eac", "6ab23c16da634ae6b08b8d86aab40045", "a8501bedde05b707f60137f070ed44ca", "4d477616bf797ee3e3e1bdd6e991af81", "583de969439f5250046f667fc2952296", "61b75384744a590810cee684d2159dfe", "ca8af0355a6c4a57da5a6b5d8c1d0ead", "f4b9ce1995b7a1b4d955bdd7febc048d", "838c3a100b0ce51e700c9eaef44fd63f", "20ca081d77d005e4ad2aad949198b064", "b97b54780c6ef30096df345df36bc77c", "3166af597b3a0b5f3bf00cd690bf5b10", "070a982fda821e55fcfee00c733b4c85", "e3778beb7182dfac9189313c52e629f3", "9514cf5fb85a0712b85ed4c165ab7020", "1f75973fa1e1fb7959ee45d034bd15cb", "e6af28c962cfd30d3052f00e0706f97c", "4e704d948e667d11ec4a2c9e22451604", "258c93939a51d077e9d111373552e310", "5ba00f94f78fef64c615d880c40b0be9", "6441eb142387c56634ab8f5883641af5", "b96f2fe318f82fe8c8f4bfcb8872a657", "d71b03774a2f57be7c1550a2604e1f4a", "39c2cb03812e140da52f2fc74f578090", "55184eb405dcb30589220bc1171e9d1b", "5fe94665fe76cad4636694e169fd0049", "364fd67a3999c30d84c2b96cd32f6e2b", "b70eed3cb6044e3e464cd732c7ab5a7d", "60f83339aecc5ebd9232d36bca5e8ee9", "7f4842e85ff7a3f3ad8fb7e2939c45dd", "4094c930c52af53bbcada5e1a9622939", "e0699bf4b97408400c9756919af9a414", 
                          "d4eef1aa5434dcecda92fa1b156e3a39", "4506c01a3c290d746bf739adb424911e", "6276ead64c7b18e4b406827eb9d40f4f", "e3497ef5c01f7388829d39097a65d0e3", "fca78a33aee9feb950e10af01564ddef", "393c6e4d02ee04c77806020cb8f92057", "cf1cdfbf42dba4b852b42f741bd944c1", "a4a1975646c1a614214d294b42b96475", "9ca9c8688ed137c65c6d249d77077693", "207249aecaa55d00520026e59d02b598___FCAImmP7316890-", "70d8cea8e43c8eaca1128d48f330eb28", "ce872f4958be694ebf5cc8a0715132c0", "5cf422975e2d55c24152556d8ec2510c", "987b353300b874bc9da22c6b93e895f2", "97cab4e76d7951ac5b47e2de9430f0e4", "302c8b0edcca4a261da9e315351eb244", "82a20725d25730d5a3d15fc9ef336da1", "16d1107cf16a3ca3a65c2f66cae3224a", "8964d543c0e6e04ee2f6feba07b43d5d", "7172ff96c42f0adc802ebf4f77924f98", "3dd2a9b38849c75a9d022ac2aacac4a4", "657ae6baf28bbdc170d0707ca05f6caf", "cfcc374503bb8999887648680c5c3161", "399c8a38529244a07bbd271c037123bb", "e7114ee1672e840dd42ede26f15de43b", "8ff2b9dd42cfabb2b65c311ddf98ec69", "9295d3acfde84cf0d76beeca1eb12c6e", "ca7daa44f9215e2d5c72c5a1bfb8751e", "73298a675adda397316bc8e1dbbf55c8", "ba2acd872e3d606bbbaefc434a82576f", "2ec1d2a694f0370e9a250efe42204b86", "62223e873a37e9338b263043146435e4", "2fa7b79220dca1ebfd9a453c84580a73", "5aecce4af3842f6e266d3d89a8ce6e79", "3c6d63350ddd59954f44ca5c24d4e200", "5aa15aeb66a3fb4768b476560a360666", "bd19e0ba8d739f00bd49c5b5143cdee8", "a4e74b16d388b4d59d07f667b1582ee7", "14c4bde687a72b28aa96484b84843ae6", "b5e95ecb1453e57be00aacce57401b92", "e2434e71faa99e40d7bec6f7180dbca2", "08fde49dc006717786b79ccc4d9a2441", "a7dfedfb1905efed4f50aa5d28ac16cf", "349503c836d39477ce7717b2ad6c502c", "97def5a946d9c23a530fb392e277d8df", "3a3894257091828eb650b6b2111aad84", "90cbcbcb297c52a2600307e247f2a803", "80a07adb4fa177c07f623d335592a7e2", "1c42caab0613a39d82b45865f8a06aed", "4c92542fb18502973db5d51e70405054", "1fa12ff8dfbc926fdf195b4a9c3e026f", "539d0b265af0bd5dbaae107f1e3ac1d0", "6b856cd299752bf2696ff215673164c9", "d00f79ea2739eee72fa8d4d2e5df5307", "4c886aaa992a4f80309ef71b724e5a3d", "1c2941f873166a93acbb5f2bca7e10b7", "e56ea15ffb15f4e7a0878365ff17f774", "dfa15cd141844cf920392f197f5c67e0", "dfc19bf93a495e22cdd2bef76622259e", "d2582aa250d6c0403b0f7461c7b3ef33", "c24354cd630414430f81fced722ae8bf", "29ffb816b9082af76e105dc8bd8e0a27", "5a390b28994b4a058c339e86e59a968f", "c27200e3d0ee3de6542d8d6b3df6c7cd", "d4771f12dc4f58f2705fdcc4efc75a72", "5d6707bfdc33b51e2d75b5db856d5173",
                          "02539c451b9b43c13bf77c73e28a55f7", "8b268e44fe9b6da2aef755cfb2e1b509", "1eeda97bc1d852226d6d8e4af5bf12f9", "0c5b9c139029243d9baeea27848370ff", "d58edbcac74e93e226767d2b943d8257", "074a116aef9ea51cfa9aa7a71604d216", "241b6d9b5bcf1ae34eb7afcbf5631ad6", "14c498ed61ef2ebcd88eae2de3de0619", "7f7e90f721666cf2fa38cf2c352651ce", "b1ddacd6cc25ca7fe7a422e573a351c8", "38e950cfdce9a23f33d29a742e71d245", "b64c5b8b32064c7ce0182920c3906fb9", "8935da78a4021eb909d998af3769cfa5", "97cc420a6a1e882373325d6741b57de8", "dbba8fcc4ba497bbcd375918f309167a", "36ebd0e1f94296b5f97acf25a7643e43", "7469401c85e883ce428351b06b0906aa", "d9d3f786228d77b6957c9ac9ef7e769b", "c9df832723b0ec9e5c3c07d976854a21", "a14118234aff6dcfbf5ce0dd797532c5", "024e72d66737e43e40f79c51222f3840", "e3ee12d2c91651f69f78e800fba73978", "1714b3d8d53c9e8b840c317cb4a5d485", "fdac827cf31de1c092b1c435f4df295d", "e888b654b0126267aac935bd2db5b738", "0e4e054429a1737630632dd44110af49", "90461588050dd4181b32ae11de90e681", 
                          "1cf4d4c470f3a0a022f2bbb03fe65fc8", "240457c3dfb5ac5454702b2d019b8a09", "118f34abe3a4458c7e2ee51004265890", "9651ed4e85afe95a13017c1c6348e031", "1a17ec38cd2c8f234f93ec568a902c59", "3f3746053c23d09970832816c21ba6d7", "cfb9d03494c9c66586f72b690bc43ce3", "0e2e0ae153c49263b8466c8382a226bf", "9060e31c46cc60bb5ed7a61e6240b028", "f5201385e1aff4a957492af5b3ec91da", "25ab8fb89cf4cf2e624b20b36c10df68", "0a6d36c5a716888dd403709ce35fd393", "87a882e7854112c097254660f8d324ef", "1cc0b9f485e5a3d0de07767dfb7026bf", "428136d30292933c81f3f161bb5437a5", "f4a8cf082680783b1b61fa494ff52ebc", "8c117bedc3760c267e44c293e8d6ff5b", "16d253d43fd6324a850c802652e0462b", "273b5a1aa54730626a03e0eb2844354a", "0f1fe42caa6ac9d348b01a85ae0b01b9", "63eb22531e45d4a41e445bf0dcb6c719", "51ff7f64d367964417fed13f7673b14f", "be13a93d25f0302fbc4c498d8b98c566", "9ab879ce7803ef6ad68551e380b5f022", "7f3865e19996888502c5627b613c5ca6", "64ccf4f743468156e76462dd502da658", "d6c5f4cb42f14b590c1dbc37b161dcf5", "ce433a0f4bff8ea0b392c182f9edeec7", "9242cc6d7b31da1480a85226334c7cbe", "f995cb92736f525d78b90d6ad69c4c3e___-1-HCATisStab7413620", "f995cb92736f525d78b90d6ad69c4c3e___-1-HCATisStab7413619", "4e0c35c565ace058abd9027e5035601e", "9ff86a4038b23465961fb460c91c42c0", "422e647c4b2df501c536f966ecb81a4d", "d115ad590a3f7f54d2424b43c2830f30", "18e8a6d37a029f70ac1b5e0edf8ecbd8", "27026e759f151d0fe2eda6c3fe993680", "c5d1d69b6157f08e8a2b13380715117c", "b05abee40848521c700ce394c6dcb4d3", "7d55ceef1e549e71d0d39a79632520ff", "4358abc15916da212c15e56d0d39e2cd", "b2ed45002d900d2e94eb895e42333804", "b6cbeb466146e2b0b713a8383dfc8c0a", "183f045c5606ede7e3ce7fb1596589f4", "1a8a75ff190329bb5f7ec3fd87f4a6e8", "a398e7fb67563fd80fbfb879cfad4c1b", "27064ca4708350a9864dbf40db4d1ee8", "b49fefa280c0ad015077b92f914a7749", "d7ec7e2b8dc40f79ed42f90758282a18", "d3303cd670d68f4cc82ed47e04fb8f7a", "235423e2d31166069cfab31be978a81f", "cd964a98e934182e6f995d2f84468051", "1c722670b114e28b1742c52f7915f14e", "284ea367345cc44e551edac78ea56791", "7ee41772f59c60c972cad75673594aca___FCAImmP7862087-", "cee3b5dbcbb2b20795a9e93c98c98345", "40c8fc36041be47b3c9ab08e4fe16a7d", "96be639e4d4262bd09db7dfd20e04810", "45d7bbbbe2ebd6ad0abcd235015b559d", "8da4d6f30674f8d16a630dbc8559795f", "36f2d472fa391a370cda234b04b51e10", "1f98af634d34678682bc4fdbf221c248", "1db95c8096d0b798235351d95af61855", "2251ebbfea6de0fbd718d75fc6195940", 
                          "fcfff1d61536c3cffe2ca4b4cc9d5dea", "625062d55df8b90649f822ec811f671c", "19296970a7d7c3a85c54aa49f47fb433", "91fee01792f1382a788f840da903165a", "8d4c311729dc43b007616d7211832663", "22a0b71b77afd4ceaca8153b9ddf2d04", "8fd54d1a00fa1a49ffc55d3cff5aea4f", "71cc68ef783fb86b3bb8682fbe939c12", "3d26c1810cba809e69d5e091ed164759", "ffa16ac41fc1b7e236b0baa625ebf241", "958937baad440506ea1de66c7ed6b43d", "36a5b792bcc036a168ff2c28de410a45", "1f95ed2a0ff9de65e471d0c09cbbaf18", "1878b3379a1f9e6af8ea41f553f2937f", "e6b8b48918344da7509c98fcc2644221", "c3c89c737c1d033dba2bf88d7571d4a4", "e3710fee478740e7dc87b1d215f2525e", "12cfb09b51703ac587f3dce5355de7ae", "18d6ed0ee0adb037f7aeeb2b3535a3b2", "be8b45c6bb187fd223c3f7f401c27c72", "db4f8896d7e7498501359ec3213cbd1a", "b89d1738326fe6d5f53aeff629f9fb93", "9b22ce15958ca94a6f25a8eb20bb4692", "d2ea5c7b67d710fea492da47049bebf8", "5e81d71af0615e6c73b60dc45778671c", "62d34de86b676a71a2914b6cc4fa9095", "7afb05911ab53121ff0a6aa321c5b27e", "fe0af29a98f40737d30ab158105ffaff", "f72d2ada58562916020f420f22096a2f", "665542cfe33744f194337bf8c2d7c763", "d5317e77c9da1ec489d46bfcd2fe82f8", "4d17514827b9534bb1c84463d39ef13a", "0d8f0044c1f56071241e7675c7c64bae", 
                          "afb1d4d16fd292b3454f0f419b30a90a", "bfd9b7890f769b33a1b553763322cafb", "af5aad8154419de2396cced6e4b9648f", "4b19bc3e7421f9242551a74c72121269", "c239fbead2b3abc0318552f69f5f6e4a", "78114000e2efa31a4edc04a45910d356", "a02378d98ba7b6b92ba6e59b298ca5af", "f9637eaf20898cd7b198c94bfe48c46b", "bea4a1260653f8fa77dc9effdf72606c", "1bd5375fae02391f6009f3e8a70d4bcb", "8021d84613ff83b4b372e992d49e77e6", "7df3cc4d8d967c428ecad6ab70d8eea6", "7aa56d1bccedecbc81868cea7c7c0aee", "47be0795c53d7b7b71854b4da7b70d5d", "1ee4820bbeb18d3bd9830736dce01d0a", "7d51b1a3d2fe51f4cfcc2906f7ebc48d", "9cea5beacda424ef29980db1459a9836", "f995cb92736f525d78b90d6ad69c4c3e___-1-HCATisStab7413621", "44e503000646216b45cad4edff1946d9", "bb123b8c380533d7f9d0e6223406a682", "3b641516862bc2d1926b238fcd9ee4db___-1-HCATisStab7646028", "94d9ffec5b69557a36c45c0adca2a9d8", "106a4220990b91f9e70aac4861c42996", "79152a1ac421038cd43f0383ff656f61", "a42f0503c52b43eb380ccae05dc766dc", "993ac2a303b16e3cb8a66c2f42017b47", "0a71a027abc14dcfc90d14da8a65e627", "242dc223ab826598099e5a8d12d8d6a2", "b5842c33696d92923859e208f8321466", "7c68364a1a216f98cb016aafd9c9645e", "d6d0df06e3293400fe4d2a5ae9b31624", "e84ffb6e9b2cfeb4e66f067e0b16598c", "54317c514ecb4a7ef0c370ccd1013bd7", "2421643c2750a36fe414e04c7f167170", "c21d2ba633e6bdf7ee99c470569f2844", "01dd8d82ec19909b9180093aeb3640b0", "6a23d18d2fdbadc93564a487a1f4319a", "993ea3ecb6892d189a001f7431919bdd", "daec721547b6302663b4354d5c714cc2", "22c96aa29fd6fdbf900de96a5b70e8af", "00cecc59f936c8f3025d1993043fbca6", "8d2aeb5d5fbb79cd3736bd54d003e19f", "8b112af4c495476f4a268f732149a6a9", "21007bca7b083b6ec6207eb78815c373", "4cb02eaa3c8a00d62f47a7da2e502a4a", "ed9c05892053e540dbd9584a50cf621f", "a6389e77e7e661820dc385f2918841c0", "7c37aba93b528faf45debdb99c826c89", "ecb743062e2585d4b2c88b61f98db5e9", "3fe7fa1f2892e625c32aa4ed6389d8dc", "0aadf09cee9c3bbe744141e14c820a15", "e0a54684907153b66cb7ba89c9747776", "3fe8f9f419ecd023cffcd58621552b69", "a9c78a346b4124e5208de7a975b54e9c", "6e3892cdec9470897782e791c881e53e", "0b1b03b7f09796c932f34d2b115c378d", "e0e6cf4fd0a36c273dba3ed640c08bad", "e955e56b863ef995ff18ae09e22cd5dd", "99bd08be839586afc1358764b41218c2", "39eb5088369cfdd2bbf78cc46f50e1fa", 
                          "4c89ffa9f9d44ef3176d5adf19606cf0", "297302c83711753362d9018f18381965", "e0b90308fe5a1da07afbf82b72c20f88", "b3ef176e44ebf27fc90f16465f856a07", "53bf4ff6c378a0600a1913860f050cc8", "dc3362a3056fdb217ec272a14a78ab7c", "005a33802ba83df104f253afbd8035c0", "ecac2b8ac7dee46a30a55496b1b6189e", "c75bd3a02b0db21ae1caefc645e3b0ab", "9653fa26190ff57ff86f6e205bd8144c", "79596a0f7148dd5e01910384ed12e510", "7c534317e7b9cdd38a3af6bffee83a15", "03ec4c87400d89980a01ce49dbd0e4ad", "e3c5a00516ae7c9a7565587a9d0cf6b0", "d3495cf086a605dfd92154062f329323", "728a70a57ade816680da325e1aa14310", "91b87b9b21e7b34caae1da5e62447d0c", "5ed6fc53d24804bc032ddb33c6ff4934", "7ee41772f59c60c972cad75673594aca___FCAImmP7862084-", "51386fed8c93108c5499e66c66965ba7", "db665461e4eae719dc1c1063c62300e7", "750437b8c590e5ee58ee265b497d0d09")
  ) |> 
  filter(!sample_id %in% c("90461588050dd4181b32ae11de90e681", "c75bd3a02b0db21ae1caefc645e3b0ab", "e029658622d33a00eda3c226e4b7aa62")) |> 
  filter(cell_type_unified_ensemble %in% c("epithelial", "stromal", "secretory")) |>
  # add_count(sample_id, cell_type_unified_ensemble) |>
  # filter(n>30) |>
  
  # filter(file_id_cellNexus_single_cell=="004e0dd96de6f3091dac2cf8cc64ddc4___1.h5ad") |> 
  
  get_single_cell_experiment(assay = "counts", atlas_name="cellxgene", cache_directory ="/vast/projects/cellxgene_curated/cellxgene")
#get_seurat(assay = "X", atlas_name="cellxgene")

sce |> assays() |> names() = "counts"

sce = SingleCellExperiment(assay = sce |> assays(), colData = colData(sce) |> apply(2, as.character))

sce |> 
 # select(.cell, sex, age_days, dataset_id, observation_joinid , sample_id,  contains("cell"), disease, collection_id, title, contains("ethnicity"), assay, tissue, tissue_groups) |> 
  select(-run_from_cell_id) |> 
  edit_covariates() |> 
  zellkonverter::writeH5AD("~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/single_cell_for_dharmesh.h5ad",  verbose = TRUE)
system("~/bin/rclone copy ~/PostDoc/immuneHealthyBodyMap/HCA_pipeline/single_cell_for_dharmesh.h5ad box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")




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



