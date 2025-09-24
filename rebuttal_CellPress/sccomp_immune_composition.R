

job::job({
  
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
      workspace_on_error = TRUE, 
      # workspaces = "estimates",
      # format = "qs", 
      
      
      #-----------------------#
      # SLURM
      #-----------------------#
      controller = crew_controller_group(
        
        
        
        crew_controller_slurm(
          name = "slurm_1_80",
          tasks_max = 1,
          workers = 40,
          seconds_idle = 30,
          crashes_max = 7,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = 80, 
            cpus_per_task = 30, 
            verbose = T
          )
        ),
        crew_controller_slurm(
          name = "slurm_1_200",
          tasks_max = 1,
          workers = 20,
          seconds_idle = 30,
          crashes_max = 7,
          options_cluster = crew_options_slurm(
            memory_gigabytes_required = 200, 
            cpus_per_task = 2, 
            verbose = T
          )
        )
      ),
      debug = "input_relative_sample",
      
      resources = tar_resources(crew = tar_resources_crew("slurm_1_80")) 
      #, # Set the target you want to debug.
      #
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
        
        # left_join(age_bin_table, copy=TRUE) |> 
        
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
    
    create_input_cell_counts = function(
    drop_sample_df, 
    caq_celltype_level_map,
    ethnicity_imputed, count_by = c("sample_id", "donor_id"), 
    sce_scored_delayed_plasma,
    result_directory
    ){
      
      count_by <- match.arg(count_by)
      
      my_tbl = 
        get_metadata() |> 
        # tbl(
        #   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
        #   sql(glue("SELECT * FROM read_parquet('{cellNexus_metadata}')"))
        # ) |>
        
        # Filter low quality cells
        dplyr::filter(
          empty_droplet == FALSE,
          alive == TRUE,
          scDblFinder.class != "doublet"
        ) |>
        
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
        
        # NON immune cells
        mutate(cell_type_unified_ensemble = if_else(is_immune, cell_type_unified_ensemble, "non_immune")) |> 
        
        # IMMUNE CELLS
        # filter(is_immune) |> 
        filter(cell_type_unified_ensemble %in% c("non_immune", "cd8 naive", "cd16 mono", "cd4 tcm", "cd4 th17 em", "granulocyte", "cd4 th1/th17 em", "treg", "b memory", "b naive", "nk", "plasma", "cd4 th2 em", "mast", "cd4 th1 em", "cd8 tem", "mait", "tgd", "cdc", "cd4 fh em", "cd4 naive", "nkt", "macrophage", "cd8 tcm", "cd14 mono", "pdc", "ilc")) |> 
        
        edit_covariates() |> 
        
        # Here we drop those samples with a low cell type entropy. E.g. one cell type only.
        anti_join(drop_sample_df, copy = TRUE) |> 
        
        
        # Attach plasma subtype
        left_join(
          sce_scored_delayed_plasma |> 
            readRDS() |> 
            as_tibble()
        ) |> 
        
        # calculate median per dataset
        mutate(median_plasma_score = median(TotalScore, na.rm = TRUE), .by = dataset_id) |> 
        mutate(cell_type_unified_ensemble = 
                 case_when(
                   cell_type_unified_ensemble == "plasma" & !TotalScore |> is.na() & TotalScore > median_plasma_score ~ "plasma_long_lived",
                   cell_type_unified_ensemble == "plasma" & !TotalScore |> is.na() & TotalScore <= median_plasma_score ~ "plasma_short_lived",
                   TRUE ~ cell_type_unified_ensemble
                 ))
      
      
      if(count_by=="donor_id")
        my_tbl =  
        my_tbl |> 
        
        # I add dataset_id because there are duplicated donor_id
        # I add technology because there are 71 donors with multiple technologies. It is acceptable in the grand scheme of things of 4K+ donors
        # Same for disease
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
          caq_celltype_level_map |> 
            
            # Plasma hierarchy
            mutate(L4=L3) |>
            bind_rows(tibble(
              cell_type_unified_harmonised = c("plasma_long_lived", "plasma_short_lived"), 
              L0 = "b",   L1 ="plasma",       L2   ="plasma",     L3  = "plasma",    
              L4 =c("plasma_long_lived", "plasma_short_lived")
            )) |> 
            mutate(L4 = case_when(
              L4 |> is.na() ~ NA,
              cell_type_unified_harmonised %in% c("plasma_long_lived", "plasma_short_lived") ~ L4,
              TRUE ~ "non_plasma"
            )),
          
          by = join_by(cell_type_unified_ensemble == cell_type_unified_harmonised)
        ) |> 
        mutate(
          across(
            matches("^L[0-9]"),        # select columns whose names start L0, L1, … L9
            ~ if_else(is_immune, as.character(.), "non_immune")
          )
        ) |> 
        
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
    #-----------------------#
    list(
      tar_target(
        result_directory,
        "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12", 
        deployment = "main"
      ),
      tar_target(
        drop_samples,
        {
          # evaluate result_directory for targets
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
        input_relative_sample, 
        create_input_cell_counts(
          drop_samples,
          caq_celltype_level_map,
          ethnicity_imputed,
          sce_scored_delayed_plasma,
          result_directory,
          count_by = "sample_id"
        ),
        packages = c( "dplyr", "tidyr", "tibble",   "duckdb",   "glue",  "readr", "forcats","tidybulk", "cellNexus" )
      ),
      tar_target(sce_scored_delayed_plasma, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sce_scored_delayed_plasma.rds", format = "file"),
      tar_target(
        input_relative_donor, 
        create_input_cell_counts(
          drop_samples,
          caq_celltype_level_map,
          ethnicity_imputed,
          sce_scored_delayed_plasma,
          result_directory, 
          count_by = "donor_id"
        ),
        packages = c( "dplyr", "tidyr", "tibble",   "duckdb",   "glue",  "readr", "forcats","tidybulk", "cellNexus"  )
      ),
      tar_target(
        saved_input_relative,
        {
          # evaluate result_directory for targets
          print(result_directory)
          
          file_name = glue("{result_directory}/cell_metadata_1_0_10_sccomp_input_counts.rds")
          input_relative_sample |> 
            saveRDS(file_name)
          
          check_rclone_installation()
          system(glue("~/bin/rclone copy {file_name} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/"))
          
        }, 
        packages = "glue"
      ),
      tar_target(
        formula_df,
        tribble(
          ~ name, ~ formula_composition, ~ formula_variability, ~ counts,
          
          # continuous age
          "estimates_continuous_age",
          "~ 1 + age_days_scaled*sex + disease_groups___altered + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_days_scaled*sex + ethnicity_groups_imputed | tissue_groups)",  
          "~ age_days_scaled*sex + disease_groups___altered",
          input_relative_sample,
          
          # discrete
          "estimates_age_bins", 
          "~ 1 + age_bin + disease_groups___altered + sex + age_bin:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed | tissue_groups)",  
          "~ age_bin + disease_groups___altered",
          input_relative_sample,
          
          # discrete decade
          "estimates_age_decade", 
          "~ 1 + age_decade + disease_groups___altered + sex + age_decade:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed | tissue_groups)",  
          "~ age_decade + disease_groups___altered",
          input_relative_sample,
          
          # discrete decade
          "estimates_age_decade_plasma", 
          "~ 1 + age_decade + disease_groups___altered + sex + age_decade:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed | tissue_groups)",  
          "~ age_decade + disease_groups___altered",
          input_relative_sample ,
          
          # interaction decade, this model is to test that the signal really comes from each tissue and not learned
          # at the body level, because of lack of data
          "estimates_age_decade_interaction", 
          "~ 1 + age_decade*sex*tissue_groups + disease_groups___altered*tissue_groups + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered)", 
          "~ age_decade", # + disease_groups___altered",
          input_relative_sample,
          
          # discrete decade DONOR grouping
          "estimates_age_decade_donor_grouping", 
          "~ 1 + age_decade + disease_groups___altered + sex + age_decade:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed | tissue_groups)",  
          "~ age_decade + disease_groups___altered",
          input_relative_donor,
          
          # interaction decade DONOR grouping
          "estimates_age_decade_interaction_donor_grouping", 
          "~ 1 + age_decade*sex*tissue_groups + disease_groups___altered*tissue_groups + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered)", 
          "~ age_decade", # + disease_groups___altered",
          input_relative_donor,
          
          # discrete + interaction ethnicity sex
          "estimates_age_bins_sex_ethnicity_interaction", 
          "~ 1 + age_bin + disease_groups___altered + sex + age_bin:sex + ethnicity_groups_imputed * sex + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed * sex | tissue_groups)",  
          "~ age_bin + disease_groups___altered",
          input_relative_sample,
          
          # continuous + discrete
          "estimates_continuous_age_plus_age_bins",
          "~ 1 + age_days_scaled + disease_groups___altered + age_bin*sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_days_scaled + age_bin*sex + ethnicity_groups_imputed | tissue_groups)", 
          "~ age_days_scaled + disease_groups___altered",
          input_relative_sample,
          
          # disease tissue specific
          "estimates_age_bins_disease", 
          "~ 1 + age_bin + disease_groups___altered + sex + age_bin:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_bin + disease_groups___altered + sex + age_bin:sex + ethnicity_groups_imputed | tissue_groups)",  
          "~ age_bin + disease_groups___altered",
          input_relative_sample
        ) |> 
          expand_grid(
            cell_type_level = glue("L{0:4}") |> as.character(), 
            drop_disease = c(TRUE, FALSE),
            immune_only = c(TRUE,FALSE)
          ) |> 
          
          # Keep NON immune for one model only
          filter(immune_only | name == "estimates_age_bins") |> 
          
          # Filter immune if needed
          mutate(counts = map2(counts, immune_only,  ~ { 
            if(.y) .x |> filter(is_immune) 
            else .x |> 
              nest(data = -c(sample_id, tissue_groups)) |> 
              filter(map_lgl(data, ~ .x |> filter(cell_type_unified_ensemble=="non_immune") |> nrow() > 0)) |>  
              mutate(non_immune_count = map_int(data, ~ .x |> filter(cell_type_unified_ensemble=="non_immune") |> pull(n))) |> 
              mutate(total_count = map_int(data, ~ .x |> pull(n) |> sum())) |> 
              mutate(proportion_non_immune = non_immune_count/total_count) |> 
              filter(proportion_non_immune > 2/3) |> 
              filter(!tissue_groups %in% c(
                "bone marrow",
                "lymphatic system",
                "spleen",
                "thymus", 
                "blood"
              )) |> 
              unnest(data)
          } )) |> 
          
          # Drop Disease if needed
          mutate(
            counts = map2(counts, drop_disease, 
                          ~ {
                            if(.y) .x |> filter(disease_groups___altered == "Normal")
                            else .x 
                          }),
            formula_composition = if_else(drop_disease, formula_composition |> str_remove_all("\\+ disease_groups___altered"), formula_composition),
            formula_variability = if_else(drop_disease, formula_variability |> str_remove_all("\\+ disease_groups___altered"), formula_variability)
          ) |> 
          mutate(local_file_name = glue("{name}___{cell_type_level}___disease_{!drop_disease}___immune_only_{immune_only}")) |> 
          
          # TEMPORARY
          # FILTER FOR JUST ONE MODEL
          filter(name %in% c("estimates_age_decade", "estimates_age_decade_plasma", "estimates_age_decade_interaction", "estimates_age_decade_donor_grouping", "estimates_age_decade_interaction_donor_grouping")) |> 
          filter(cell_type_level %in% c("L0", "L3") | (cell_type_level == "L4" & name == "estimates_age_decade_plasma")) |> 
          
          group_by(local_file_name) |> 
          tar_group(), 
        iteration = "group",
        packages = c("tibble", "glue", "targets", "dplyr", "tidyr", "purrr", "stringr")
      ),
      tar_target(
        estimates,
        
        formula_df$counts[[1]] |> 
          
          # With L0 I have to summarise further because I have counts already
          with_groups(
            c(sample_id,age_days_scaled,age_bin, age_decade, sex, disease_groups___altered,ethnicity_groups_imputed, assay_groups___altered, dataset_id___altered, tissue_groups, all_of(formula_df$cell_type_level)),
            ~ .x |> summarise(n = sum(n))
          ) |>
          # 
          sccomp_estimate(
            formula_composition = formula_df$formula_composition |> as.formula(),
            formula_variability = formula_df$formula_variability |> as.formula(),        # Differential variability
            sample = "sample_id", 
            cell_group = formula_df$cell_type_level, # A level of the hierarchy
            abundance = "n",
            cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1)),
            mcmc_seed = 42,
            verbose = T, 
            bimodal_mean_variability_association = TRUE,
            prior_mean = list(intercept = c(0, 0.8), coefficients = c(0, 3)),
            prior_overdispersion_mean_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40)),
            output_directory = "/vast/scratch/users/mangiola.s/my_draws", 
            
            max_sampling_iterations = 5000,  
            
            # # TEMPORARY DEBUG
            # max_sampling_iterations = 1,
            # warmup_samples =1,  
            
            inference_method =   "hmc",
            refresh = 1
          ),
        pattern = map(formula_df),
        resources = tar_resources(crew = tar_resources_crew("slurm_1_80")),
        error = "continue", 
        packages = "sccomp"
        
        # TEMPORARY
        , cue = tar_cue(mode = "never")
      ),
      tar_target(
        saved_and_tranferred,
        {
          # evaluate result_directory for targets
          print(result_directory)
          
          local_file_name = glue("{result_directory}/{formula_df$local_file_name}.rds")
          local_file_name_FIT_FOR_PORTABILITY = glue("{result_directory}/{formula_df$local_file_name}_FIT_FOR_PORTABILITY.rds")
          
          # Save draws as monolythic
          attr(estimates, "fit")$save_object(file = local_file_name_FIT_FOR_PORTABILITY) 
          
          check_rclone_installation()
          system(glue("~/bin/rclone copy {local_file_name_FIT_FOR_PORTABILITY} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_estimates_1_0_10/"))
          estimates |> attr("fit") = readRDS(local_file_name_FIT_FOR_PORTABILITY)
          
          # Save sccomp estimates
          estimates |>  sccomp_test() |> saveRDS(local_file_name)
          
          check_rclone_installation()
          system(glue("~/bin/rclone copy {local_file_name} box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_12/"))
          
        }, 
        pattern = map(formula_df, estimates), 
        resources = tar_resources(crew = tar_resources_crew("slurm_1_200")),
        error = "continue", 
        packages = c("glue", "sccomp", "magrittr")
      )
      
    )
  }, 
  ask = FALSE, 
  script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets.R"
  )
  
  tar_make(
    # callr_function = NULL,
    script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets.R", 
    store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets", 
    reporter = "verbose" #, callr_function = NULL
  )
  
})


system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins___L3___disease_TRUE___immune_only_TRUE.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_12/")


library(targets)
tar_meta(store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets") |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  dplyr::select(name, error)






library(targets)

x = tar_read(caq_celltype_level_map, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets")


tar_workspace(create_input_cell_counts, 
              script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets.R", 
              store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets"
)

library(tidyverse)
library(sccomp)
library(magrittr)
library(glue)
library(forcats)
library(stringr)
library(arrow)
library(dplyr)
library(duckdb)

tar_meta(starts_with("estimates_"), store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets")

# Age proportion prediciton
estimates_age_bins |> 
  sccomp_predict(
    formula_composition = ~ 1 + age_bin*sex + (1 + age_bin*sex | tissue_groups), 
    number_of_draws = 100,
    summary_instead_of_draws = TRUE
  ) |> 
  mutate(age_bin = factor(
    age_bin, 
    c("Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior"),
    ordered = TRUE
  )) |> 
  mutate(age_bin_numeric = age_bin |> as.integer())  |> 
  saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/prediction_age_bins.rds")

system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/prediction_age_bins.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

tar_read(formula_df, store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/_targets")


# For Hong
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins___L3.rds")
estimate_age_bins = estimate_age_bins |> dplyr::select(-count_data)
attr(estimate_age_bins, "fit") = NULL
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins_effect_tibble_only.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins_effect_tibble_only.rds box_adelaide:/immune_map_disease/")

# Save fit
library(magrittr)
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins.rds")
estimate_age_bins |> attr("fit") %$% save_object(file = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins_FIT_FOR_PORTABILITY.rds") 
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins_FIT_FOR_PORTABILITY.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
estimate_age_bins |> attr("fit") = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins_FIT_FOR_PORTABILITY.rds")
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")


# estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/estimates_age_bins.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_12/21_11_2024_sccomp_archive_before_factor_ordering/estimates_age_bins.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

# Benchmark
tic()
estimate_age_bins |> sccomp_test(contrasts = c(  "respiratory system" = "sexmale + `sexmale___respiratory system`",
                                                 "blood" = "sexmale + `sexmale___blood`"))
toc()


estimate_age_bins |> 
  sccomp_test()



# annotate plasma
out_dir <- file.path(getwd(), "./", "plasma_h5se")


build_plasma_sce <- function(cache_dir,
                             subset_n = 100,
                             min_features = 8000,
                             seed = NULL) {
  metadata <- cellNexus::get_metadata()
  
  plasma_metadata <- metadata |>
    dplyr::filter(cell_type_unified_ensemble == "plasma") |>
    dplyr::filter(
      empty_droplet == FALSE,
      alive == TRUE,
      scDblFinder.class != "doublet",
      feature_count >= min_features
    ) |>
    dplyr::collect()
  
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(subset_n) && nrow(plasma_metadata) > subset_n) {
    plasma_metadata <- plasma_metadata |>
      head(subset_n)
  }
  
  sce <- cellNexus::get_single_cell_experiment(
    plasma_metadata,
    cache_directory = cache_dir
  )
  
  sce
}


compute_singscore_delayed <- function(
    sce, assay_name = NULL, upSet = singscore::tgfb_gs_up, 
    downSet = singscore::tgfb_gs_dn, workers = 1, output_file = NULL
) {
  assay_names <- SummarizedExperiment::assayNames(sce)
  if (is.null(assay_name)) {
    assay_name <- if ("logcounts" %in% assay_names) "logcounts" else assay_names[[1]]
  }
  message("Ranking genes using assay: ", assay_name)
  expr_mat <- SummarizedExperiment::assay(sce, assay_name)
  feature_ids <- rownames(expr_mat)
  
  # Debug information
  message("Expression matrix dimensions: ", nrow(expr_mat), " x ", ncol(expr_mat))
  message("Expression matrix rownames length: ", length(rownames(expr_mat)))
  message("Expression matrix class: ", class(expr_mat))
  # Choose ID space (prefer best overlap with matrix rownames)
  up_sym <- upSet; down_sym <- downSet
  # Map symbols -> ENSG
  up_ens <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = up_sym, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
  down_ens <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = down_sym, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
  
  # Use delayed ranking with optional output file
  ranked <- singscore::rankGenes(
    expr_mat, 
    workers = workers 
    
    # TEMPORARY
    # , stableGenes = up_ens    
  )
  
  # Convert character vectors to GeneSet objects
  library(GSEABase)
  up_geneset <- GeneSet(up_ens, setName = "upSet")
  down_geneset <- GeneSet(down_ens, setName = "downSet")
  
  # Use regular simpleScore for now since the delayed version has issues
  # Convert DelayedMatrix to matrix for simpleScore
  # Set the stable attribute properly
  scoredf_delayed <- singscore::simpleScore(
    ranked, 
    upSet = up_geneset, 
    downSet = down_geneset
  )
  
  # Add original cell
  scoredf_delayed$cell_id = sce |> colData() |> _[,"original_cell_"]
  scoredf_delayed$dataset_id = sce |> colData() |> _[,"dataset_id"]
  
  # Return the delayed results and the ranked data
  return(scoredf_delayed)
}


# Build from cache (default)
sce <- build_plasma_sce(cache_dir = "/vast/projects/cellxgene_curated/cellNexus", subset_n = Inf, min_features = 8000) 

library(tidySummarizedExperiment)
# BiocManager::install("stemangiola/singscore@delayedArray")
library(singscore)
# Rename features using

# from sce to singscore execution
long_lived_genes <- c("SDC1", "CD38", "TNFRSF17", "TNFRSF13B", "CD28", "CXCR4", "CD69", "ITGA4", "ITGB1", "ITGB7", "ITGAE", "PRDM1", "XBP1", "MZB1", "DERL3", "HSPA5", "BCL2", "MCL1")
short_lived_genes <- c("CD19", "PTPRC", "HLA-DRA", "HLA-DRB1", "CD27", "CXCR3", "S1PR1", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "PAX5", "BCL6", "BACH2", "AICDA", "CCR9", "CCR10")

long_lived_genes <- c("SDC1", "TNFRSF17", "CXCR4", "PRDM1", "XBP1", "MZB1", "MCL1")
short_lived_genes <- c("CD19", "PTPRC", "HLA-DRA", "MKI67", "PAX5", "BCL6", "BACH2")

sce_scored_delayed <- compute_singscore_delayed(
  sce, 
  assay_name = "counts", 
  upSet =long_lived_genes  , 
  downSet = short_lived_genes, 
  workers = 20,
  output_file = NULL  # Disable output file for now to avoid HDF5 conflicts
)

sce_scored_delayed |> saveRDS("sce_scored_delayed_plasma.rds")

sce_scored_delayed |>
  as_tibble(rownames="cell") |> 
  left_join(colData(sce) |> as_tibble(rownames="cell")) |>
  ggplot(aes(x=TotalScore)) +
  geom_density(aes(color = dataset_id)) +
  guides(color = "none") +
  ylim(c(0, 30)) + 
  theme_minimal() 

plot(sce_scored_delayed$TotalScore, sce_scored_delayed_down$TotalScore)

# sce_scored |>
#colData() |>
#saveRDS(file = "dev/plasma_h5se_scored.rds")

