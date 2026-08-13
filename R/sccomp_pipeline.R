
library(targets)

tar_script({
    
    #-----------------------#
    # Input
    #-----------------------#
    library(tidyverse)
    library(targets)
    library(glue)
    
    #-----------------------#
    # Packages
    #-----------------------#
    tar_option_set(
      memory = "transient", 
      garbage_collection = 100, 
      workspace_on_error = TRUE
    )
    
    #-----------------------#
    # FUNCTIONS
    #-----------------------#
    
    edit_covariates = function(tbl){
      
      tissue_grouped = list(
        
        # Respiratory System
        "respiratory system" = c(
          "lung", "lung parenchyma", "alveolus of lung",  "bronchus",
          "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
          "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
          "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
        ),
        
        "trachea" = c( "epithelium of trachea", "trachea"),
        
        # Cardiovascular System
        "cardiovascular system" = c(
          "heart", "heart left ventricle", "heart right ventricle", "cardiac ventricle",
          "cardiac atrium", "right cardiac atrium", "left cardiac atrium", "apex of heart",
          "aorta", "coronary artery", 
          "venous blood", "anterior wall of left ventricle", "myocardium", "interventricular septum", "ventricular tissue", "basal zone of heart"
        ),
        
        "vasculature" = c("kidney blood vessel", "artery", "vein", "vasculature", "mesenteric artery"),
        
        # Umbilical Cord Blood
        "umbilical cord blood" = "umbilical cord blood",
        
        # Oesophagus
        "oesophagus" = c(
          "esophagus", "lower esophagus", "esophagus muscularis mucosa",
          "submucosal esophageal gland",
          
          # Epithelium
          "epithelium of esophagus"
        ),
        
        # Stomach
        "stomach" = c(
          "stomach", "body of stomach", "cardia of stomach"
        ),
        
        # Small Intestine
        "small intestine" = c(
          "small intestine", "duodenum", "jejunum", "ileum",
          
          # Epithelium
          "epithelium of small intestine", "jejunal epithelium", "ileal epithelium",
          "submucosa of ileum", "lamina propria of small intestine"
        ),
        
        # Large Intestine
        "large intestine" = c(
          "large intestine", "colon", "left colon", "right colon",
          "sigmoid colon", "descending colon", "transverse colon",
          "ascending colon", "hepatic flexure of colon", "caecum",
          "rectum", "appendix", "vermiform appendix",
          
          # epithelium
          "colonic epithelium", "submucosa of ascending colon", "lamina propria of large intestine",
          "mucosa of colon", "lamina propria of mucosa of colon", "caecum epithelium"
        ),
        
        # Digestive System (General)
        "digestive system (general)" = c(
          "intestine", "hindgut", "lamina propria", "mucosa"
        ),
        
        # Nasal, Oral, and Pharyngeal Regions
        "nasal, oral, and pharyngeal regions" = c(
          "nasal cavity", "nasopharynx", "oral mucosa", "tongue", "anterior part of tongue",
          "posterior part of tongue", "gingiva", "nose", "saliva"
        ),
        
        # Cerebral Lobes and Cortical Areas
        "cerebral lobes and cortical areas" = c(
          "frontal lobe", "left frontal lobe", "right frontal lobe", "primary motor cortex",
          "dorsolateral prefrontal cortex", "superior frontal gyrus", "orbitofrontal cortex",
          "medial orbital frontal cortex", "Broca's area", "prefrontal cortex",
          "temporal lobe", "left temporal lobe", "right temporal lobe", 
          "angular gyrus", "entorhinal cortex",
          "parietal lobe", "left parietal lobe", "right parietal lobe", "primary somatosensory cortex",
          "occipital lobe", "right occipital lobe", "primary visual cortex",
          "occipital cortex", "insular cortex", "parietal cortex", "temporal cortex",
          "frontal cortex", "Brodmann (1909) area 4", "temporoparietal junction",
          "middle temporal gyrus", "cingulate cortex", "brain", "brain white matter", "cerebral cortex", "cerebral nuclei"
        ),
        
        # Limbic and Basal Systems
        "limbic and basal systems" = c(
          "anterior cingulate cortex", "anterior cingulate gyrus", "hippocampal formation",
          "hypothalamus", "thalamic complex", "dentate nucleus", "basal ganglion",
          "caudate nucleus", "putamen", "substantia nigra pars compacta",
          "lateral ganglionic eminence", "medial ganglionic eminence",
          "caudal ganglionic eminence", "ganglionic eminence"
        ),
        
        # Brainstem and Cerebellar Structures
        "brainstem and cerebellar structures" = c(
          "pons", "midbrain", "myelencephalon", "telencephalon", "forebrain",
          "cerebellum", "cerebellum vermis lobule", "cerebellar cortex",
          "hemisphere part of cerebellar posterior lobe", "white matter of cerebellum"
        ),
        
        # General Brain and Major Structures
        "general brain and major structures" = c(
          "spinal cord", "neural tube", "cervical spinal cord white matter"
        ),
        
        # Muscular System (Skeletal Muscles)
        "muscular system (skeletal muscles)" = c(
          "rectus abdominis muscle", "gastrocnemius", "muscle of abdomen", "muscle organ",
          "muscle tissue", "pelvic diaphragm muscle", "skeletal muscle tissue", "muscle of pelvic diaphragm"
        ),
        
        # Connective Tissue
        "connective tissue" = c(
          "connective tissue", "tendon of semitendinosus", "vault of skull", "bone spine",
          "rib"
        ),
        
        # Adipose Tissue
        "adipose tissue" = c(
          "adipose tissue", "subcutaneous adipose tissue", "visceral abdominal adipose tissue",
          "perirenal fat", "omental fat pad", "subcutaneous abdominal adipose tissue",
          "abdominal adipose tissue"
        ),
        
        # Endocrine System
        "endocrine system" = c(
          "thyroid gland", "adrenal tissue", "adrenal gland", "islet of Langerhans",
          "endocrine pancreas", "pineal gland"
        ),
        
        # Lymphatic System
        "lymphatic system" = c(
          "lymph node", "mesenteric lymph node", "thoracic lymph node",
          "cervical lymph node", "bronchopulmonary lymph node", "tonsil", "inguinal lymph node"
        ),
        
        # Integumentary System (Skin)
        "integumentary system (skin)" = c(
          "skin of abdomen", "skin of forearm", "skin of scalp", "skin of face", "skin of leg",
          "skin of chest", "skin of back", "skin of hip", "skin of body", "skin of cheek",
          "skin of temple", "skin of shoulder", "skin of external ear", "skin of trunk",
          "skin of prepuce of penis", "skin epidermis", "arm skin", "lower leg skin",
          "hindlimb skin", "zone of skin", "dermis", "skin of nose", "skin of forehead",
          "skin of pes", "axilla"
        ),
        
        # Gastrointestinal Accessory Organs
        "gallbladder" =  "gallbladder",
        
        # Gastrointestinal Accessory Organs
        "pancreas" = c( "pancreas", "exocrine pancreas" ),
        
        # Gastrointestinal Accessory Organs
        "liver" = c( "liver", "caudate lobe of liver", "hepatic cecum" ),
        
        # Spleen
        "spleen" = "spleen",
        
        # Thymus
        "thymus" = "thymus",
        
        # Blood
        "blood" = "blood",
        
        # Bone Marrow
        "bone marrow" = "bone marrow",
        
        # Female Reproductive System
        "female reproductive system" = c(
          "uterus", "myometrium", "fallopian tube", "ampulla of uterine tube",
          "fimbria of uterine tube", "uterine cervix", "endometrium",
          "decidua", "decidua basalis", "placenta", "yolk sac", "isthmus of fallopian tube"
        ),
        "ovary" = "ovary", 
        
        # Male Reproductive System
        "male reproductive system (other)" = c(
          "testis", "gonad"
        ),
        
        # Prostate
        "prostate" = c(
          "prostate gland", "transition zone of prostate", "peripheral zone of prostate"
        ),
        
        # Renal System
        "renal system" = c(
          "kidney", "cortex of kidney", "renal medulla", "renal papilla",
          "renal pelvis", "ureter", "bladder organ"
        ),
        
        # Miscellaneous Glands
        "miscellaneous glands" = c(
          "parotid gland", "lacrimal gland", "sublingual gland", "mammary gland",
          "chorionic villus"
        ),
        
        # Eye and Visual-Related Structures
        "sensory-related structures" = c(
          "retina",
          "retinal neural layer",
          "macula lutea",
          "macula lutea proper",
          "sclera",
          "trabecular meshwork",
          "conjunctiva",
          "pigment epithelium of eye",
          "cornea",
          "iris",
          "ciliary body",
          "peripheral region of retina",
          "eye trabecular meshwork",
          "perifoveal part of retina",
          "choroid plexus",
          "lens of camera-type eye",
          "corneo-scleral junction",
          "fovea centralis",
          "eye",
          "inner ear",
          "vestibular system",
          "primary auditory cortex"
        ),
        
        # Digestive Tract Junctions and Connections
        "digestive tract junctions and connections" = c(
          "esophagogastric junction", "duodeno-jejunal junction", "hepatopancreatic ampulla",
          "hepatopancreatic duct", "pyloric antrum"
        ),
        
        # Peritoneal and Abdominal Cavity Structures
        "peritoneal and abdominal cavity structures" = c(
          "peritoneum", "omentum", "retroperitoneum", "mesentery"
        ),
        
        # Breast
        "breast" = c(
          "breast", "upper outer quadrant of breast"
        )
      ) |> 
        enframe(name ="tissue_groups") |> 
        distinct() |> 
        unnest(value) |> 
        rename(tissue = value) 
      
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
      
      temp_path = tempdir()
      system(glue("~/bin/rclone copy box_adelaide:/minh_immune_map_disease/disease_data_grouped_further.csv {temp_path}/"))
      
      disease_data_grouped = 
        disease_data_grouped |> 
        left_join(
          read_csv(glue("{temp_path}/disease_data_grouped_further.csv")) |> 
            rename(disease_groups_further = disease_groups)
        ) |> 
        mutate(disease_groups = if_else(!disease_groups_further |> is.na(), disease_groups_further, disease_groups)) |> 
        select(disease,  disease_groups)
      
      tbl |> 
        
        # TISSUE
        select(-any_of("tissue_groups")) |> 
        left_join(tissue_grouped, copy=TRUE) |> 
        
        # TECH
        left_join(assay_data_grouped, copy=TRUE) |> 
        
        # DISEASE
        left_join(disease_data_grouped, copy=TRUE) |> 
        
        # make disease tissue specific, omit Normal
        mutate(disease_groups = paste(disease_groups, tissue_groups, sep = "_")) |> 
        mutate(disease_groups = if_else(disease_groups |> str_detect("Normal_.+"), "Normal", disease_groups))  |> 
        
        
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
        mutate(age_years = age_days / 365) |> 
        mutate(age_bin = dplyr::case_when(
          age_years < 3 ~ "Infancy",
          age_years < 12 ~ "Childhood",
          age_years < 20 ~ "Adolescence",
          age_years < 40 ~ "Young Adulthood",
          age_years < 50 ~ "Middle Age",
          age_years < 60 ~ "Senior_50",
          age_years < 70 ~ "Senior_60",
          age_years >= 70 ~ "Senior_70",
          TRUE ~ NA_character_
        )) |> 
        mutate(age_decade = ceiling(age_years/10) |> as.integer() |> as.character()) |> 
        
        # ETHNICITY
        left_join(ethnicity_grouped, copy=TRUE) |> 
        
        dplyr::select(
          cell_id, sample_id, donor_id, dataset_id, title, collection_id, age_days, age_bin, age_decade, sex, 
          ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble,
          cell_type, disease_groups, is_immune
        ) |> 
        as_tibble() |> 
        
        # Center based on adolescence
        mutate(age_days_scaled = age_days  |> scale(center = 50*365) |> as.numeric()) 
      
    }     
    
    join_census_from_cache = function(metadata){
      census_path = path.expand("~/Library/Caches/org.R-project.R/R/cellNexus/census_cell_metadata.2.3.0.parquet")
      if (!file.exists(census_path)) {
        stop("Census metadata not found at ", census_path, 
             ". Install/cache cellNexus census metadata, or restore join via cellxgene.census.")
      }
      
      con = dbplyr::remote_con(metadata)
      census = duckdb::tbl_file(con, census_path)
      
      metadata |> 
        left_join(
          census |> 
            select(
              observation_joinid, dataset_id, tissue, self_reported_ethnicity,
              cell_type, assay, disease, sex, title, collection_id
            ),
          by = join_by(observation_joinid, dataset_id)
        )
    }
    
    create_input_cell_counts = function(
    drop_sample_df, 
    caq_celltype_level_map,
    ethnicity_imputed
    ){
      
      my_tbl = 
        get_metadata() |> 
        join_census_from_cache() |> 
        
        # Filter low quality cells
        dplyr::filter(empty_droplet == FALSE) |>
        
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
        
        # Filter embrios
        filter(age_days > 365) |> 
        
        # NON immune cells collapsed; keep immune cell types of interest
        mutate(cell_type_unified_ensemble = if_else(is_immune, cell_type_unified_ensemble, "non_immune")) |> 
        filter(cell_type_unified_ensemble %in% c("non_immune", "cd8 naive", "cd16 mono", "cd4 tcm", "cd4 th17 em", "granulocyte", "cd4 th1/th17 em", "treg", "b memory", "b naive", "nk", "plasma", "cd4 th2 em", "mast", "cd4 th1 em", "cd8 tem", "mait", "tgd", "cdc", "cd4 fh em", "cd4 naive", "nkt", "macrophage", "cd8 tcm", "cd14 mono", "pdc", "ilc")) |> 
        
        edit_covariates() |> 
        
        # Drop samples with low cell type entropy
        anti_join(drop_sample_df, copy = TRUE) |> 
        
        # Donor-level sample_id (dataset/tech/disease kept to avoid collisions)
        mutate(sample_id = paste0(donor_id, dataset_id, tissue_groups, age_days, assay_groups, disease_groups)) 
      
      my_tbl |> 
        dplyr::count(
          sample_id, donor_id, dataset_id, title, collection_id, age_days, age_bin, age_days_scaled, age_decade,
          sex, ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble, is_immune,
          disease_groups) |> 
        mutate(n = as.integer(n)) |> 
        as_tibble() |> 
        
        # Add hierarchy of cell types L1, L2, L3
        left_join(
          caq_celltype_level_map,
          by = join_by(cell_type_unified_ensemble == cell_type_unified_harmonised)
        ) |> 
        mutate(
          across(
            matches("^L[0-9]"),
            ~ if_else(is_immune, as.character(.), "non_immune")
          )
        ) |> 
        
        # Add imputed ethnicities
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
        ) 
      
    }
    
    check_rclone_installation = function(){
      rclone_path <- path.expand("~/bin/rclone")
      if (!file.exists(rclone_path)) {
        stop("rclone was not found in the expected location '~/bin/rclone'.")
      }
    }
    
    #-----------------------#
    # Pipeline
    # Single model: age_decade × sex × tissue, donor grouping, L3, immune-only, Normal only
    # Output: estimates_age_decade_interaction_donor_grouping___L3___disease_FALSE___immune_only_TRUE.rds
    #-----------------------#
    list(
      tar_target(
        result_directory,
        {
          dir <- "sccomp_on_cellNexus_1_0_14_add_smooths"
          dir.create(dir, showWarnings = FALSE, recursive = TRUE)
          normalizePath(dir)
        },
        deployment = "main"
      ),
      tar_target(
        drop_samples,
        {
          print(result_directory)
          check_rclone_installation()
          system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/dharmesh_shared_mix/drop_samples.csv {result_directory}/"))
          read_csv(glue("{result_directory}/drop_samples.csv"))
        }, packages = c("glue", "readr")
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
        caq_celltype_level_map,
        {
          check_rclone_installation()
          temp_path = tempdir()
          system(glue("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/caq_celltype_level_map.csv {temp_path}/"))
          read_csv(glue("{temp_path}/caq_celltype_level_map.csv"))
          
        }, packages = c("glue", "readr")
      ),
      tar_target(
        input_relative_donor, 
        create_input_cell_counts(
          drop_samples,
          caq_celltype_level_map,
          ethnicity_imputed
        ),
        packages = c("dplyr", "tidyr", "tibble", "duckdb", "dbplyr", "glue", "readr", "forcats", "tidybulk", "cellNexus")
      ),
      tar_target(
        estimates,
        {
          counts = input_relative_donor |> 
            filter(is_immune) |> 
            filter(disease_groups___altered == "Normal") |> 
            with_groups(
              c(sample_id, age_days_scaled, age_bin, age_decade, sex, disease_groups___altered,
                ethnicity_groups_imputed, assay_groups___altered, dataset_id___altered, tissue_groups, L3),
              ~ .x |> summarise(n = sum(n))
            )
          
          sccomp_estimate(
            .data = counts,
            formula_composition = ~ 1 + age_decade * sex * tissue_groups + ethnicity_groups_imputed + assay_groups___altered +
              (1 | dataset_id___altered),
            formula_variability = ~ age_decade,
            sample = "sample_id", 
            cell_group = "L3",
            abundance = "n",
            cores = parallelly::availableCores(),
            mcmc_seed = 42,
            verbose = TRUE, 
            bimodal_mean_variability_association = TRUE,
            max_sampling_iterations = 5000,
            inference_method = "hmc",
            refresh = 1
          )
        },
        packages = c("sccomp", "dplyr", "tidyr")
      ),
      tar_target(
        saved_and_tranferred,
        {
          # Touch estimates so targets registers the dependency
          est <- estimates
          print(result_directory)
          
          local_file_name = glue("{result_directory}/estimates_age_decade_interaction_donor_grouping___L3___disease_FALSE___immune_only_TRUE.rds")
          local_file_name_FIT_FOR_PORTABILITY = glue("{result_directory}/estimates_age_decade_interaction_donor_grouping___L3___disease_FALSE___immune_only_TRUE_FIT_FOR_PORTABILITY.rds")
          
          attr(est, "fit")$save_object(file = local_file_name_FIT_FOR_PORTABILITY) 
          
          check_rclone_installation()
          system(glue("~/bin/rclone copy {local_file_name_FIT_FOR_PORTABILITY} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_estimates_1_0_13/"))
          attr(est, "fit") <- readRDS(local_file_name_FIT_FOR_PORTABILITY)
          
          est |> sccomp_test() |> saveRDS(local_file_name)
          
          check_rclone_installation()
          system(glue("~/bin/rclone copy {local_file_name} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_14_add_smooths/"))
          
          local_file_name
        }, 
        packages = c("glue", "sccomp", "magrittr")
      )
      
    )
}, 
ask = FALSE, 
script = "_targets.R"
)

tar_make(
  script = "_targets.R", 
  store = "_targets", 
  reporter = "verbose"
)

