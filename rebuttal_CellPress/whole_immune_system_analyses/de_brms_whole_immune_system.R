library(targets)

# cell type:
# [1] "b naive"         "cd14 mono"       "cd16 mono"       "cd4 naive"      
# [5] "cd4 tcm"         "cd4 th17 em"     "cd8 naive"       "cd8 tcm"        
# [9] "cd8 tem"         "t cd4"           "t cd8"           "tgd"            
# [13] "treg"            "monocytic"       "b memory"        "cd4 th1/th17 em"
# [17] "cytotoxic"       "cd4 th2 em"      "progenitor"      "mait"           
# [21] "nk"              "t"               "macrophage"      "cd4 fh em"      
# [25] "cd4 tem"         "cd4 th1 em"      "cdc"             "b"              
# [29] "dc"              "plasma"          "granulocyte"     "erythrocyte"    
# [33] "pdc"             "ilc"             "neuron"          "glial"          
# [37] "pericyte"        "endothelial"     "immune"          "blood"          
# [41] "muscle"          "stromal"         "mesothelial"     "epithelial"     
# [45] "liver"           "mast"            "nkt"             "renal"          
# [49] "endocrine"       "reproductive"    "secretory"       "fat"            
# [53] "pneumocyte"      "myoepithelial"   "sensory"         "lens"           
# [57] "epidermal"       "cartilage"       "bone" 

# SET Script ------

tar_script({
  
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(qs)
  library(crew)
  library(crew.cluster)
  
  # Set file path -----
  ## Phoenix HPC setting -----
  # hdf5_path = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseudobulk_sample_is_immune"
  # metadata_path = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/cell_metadata_1_0_6_sccomp_input_counts.rds"
  
  ## Pawsey setting -----
  hdf5_path = "/scratch/pawsey1192/zhanchen/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseudobulk_sample_cell_type"
  target_cell_type = "cd4 fh em"

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
        crashes_max = 7,
        options_cluster = crew_options_slurm(
          # script_lines = '#SBATCH -A saigencir003',
          script_lines = '#SBATCH --account=pawsey1192 \n#SBATCH --time=1-00:00:00 \nsource /software/projects/pawsey1192/zhanchen/miniconda3/bin/activate R_443',
          memory_gigabytes_required = c(5, 10, 20, 40, 80, 160), 
          cpus_per_task =16, 
          partition = 'work', # for pawsey
          # time_minutes = c(60*4, 60*4, 60*4, 60*4, 60*24, 60*24), 
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
          # script_lines = '#SBATCH -A saigencir003',
          script_lines = '#SBATCH --account=pawsey1192 \n#SBATCH --time=1-00:00:00 \nsource /software/projects/pawsey1192/zhanchen/miniconda3/bin/activate R_443',
          memory_gigabytes_required = c(80, 160), 
          cpus_per_task = 2, 
          partition = 'work', # for pawsey
          # time_minutes = c(60*24, 60*24),
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
          # script_lines = '#SBATCH -A saigencir003',
          script_lines = '#SBATCH --account=pawsey1192 \n#SBATCH --time=1-00:00:00 \nsource /software/projects/pawsey1192/zhanchen/miniconda3/bin/activate R_443',
          memory_gigabytes_required = c(160), 
          cpus_per_task = 30, 
          partition = 'work', # for pawsey
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
  
  get_adjusted_matrix = function(effect_removed_df, column_adjusted){
    
    column_adjusted = enquo(column_adjusted)
    
    m = 
      effect_removed_df |>
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
      rename(tissue = value) |>
      mutate()

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
      left_join(
        readr::read_csv("/home/zhanchen/From_scratch/Mangiola_ImmuneAtlas/disease_data_grouped_further.csv") |>
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
      filter(age_days > 365) |>
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

      # left_join(age_bin_table, copy=TRUE) |>

      # ETHNICITY
      left_join(ethnicity_grouped, copy=TRUE) |>

      dplyr::select(
        sample_id, donor_id, dataset_id, title, collection_id, age_days, age_bin, sex, ethnicity_groups, 
        tissue_groups, tissue, assay_groups, cell_type_unified_ensemble, cell_type,
        disease_groups
      ) |>
      as_tibble() |>

      # Center based on adolescence
      mutate(age_days_scaled = age_days  |> scale(center = 50*365) |> as.numeric())

  }

  offset_calcuation = function(se, method = 'TMMwsp', reference_sample){
    
    # Check if package is installed, otherwise install
    tidybulk:::check_and_install_packages("edgeR")
    
    # Drop genes with NAs, as edgeR::calcNormFactors does not accept them
    my_counts_filtered = se %>% assays() %>% as.list() %>% .[[1]] %>% na.omit()
    # Calcuate library size
    library_size_filtered = my_counts_filtered %>% colSums(na.rm  = TRUE)
    
    # Calculate TMM
    nf <-
      edgeR::calcNormFactors(
        my_counts_filtered,
        refColumn = reference_sample,
        method = method
      )
    
    # Calculate multiplier
    multiplier = library_size_filtered[reference_sample] * nf[reference_sample] %>% divide_by(library_size_filtered * nf)
    
    # Calcuate offset
    offset = log(1/multiplier)
    
    # Add to sample info
    colData(se)$normalisation_factor = nf
    colData(se)$multiplier = multiplier
    colData(se)$offset = offset
    
    return(se)
    
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
    
        # reference_sample ------
    # calculate reference_sample for scalig gene counts
    # this ensures calculation only done once
    tar_target(
      reference_sample,
      {
        metadata <- 
          get_metadata(cache_directory="/home/zhanchen/From_scratch/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseduobulk/")
        
        metadata <- metadata %>% 
          select(
            sample_id, donor_id, dataset_id, title, collection_id, 
            age_days, sex, self_reported_ethnicity, 
            tissue, assay, cell_type_unified_ensemble, cell_type,
            disease
          ) %>% as_tibble() %>% distinct() %>% 
          edit_covariates
        
        #---------------------------------#
        # Edit or add more filters here for analyses
        #---------------------------------#
        # filter(is_gene_shared) |> 
        # filter(is_immune & do_analyse) 
        
        se <- 
          loadHDF5SummarizedExperiment(hdf5_path) %>% 
          filter(cell_type_unified_ensemble == target_cell_type) %>% 
          filter(do_analyse) %>% 
          filter(is_gene_shared) |>
          select(-c(
            age_days, sex, ethnicity_groups, tissue_groups,
            assay_groups, disease_groups, age_days_scaled
          )) 
        
        se <- se %>% 
          left_join(
            metadata %>% 
              distinct(
                sample_id, 
                age_days, age_bin, age_days_scaled,
                sex, ethnicity_groups,
                tissue, tissue_groups,
                assay_groups,
                disease_groups
              ),
            by = 'sample_id', 
            copy = T
          ) 
        
        # TEMPORARY BECAUSE I FORGOT TO INTEGRATE AGE BINS
        # se = se |>
        #   left_join(
        #     readRDS(metadata_path) |>
        #       filter(age_days > 365) |> 
        #       mutate(age_years = age_days / 365) |> 
        #       mutate(age_bin = dplyr::case_when(
        #         age_years < 3 ~ "Infancy",
        #         age_years < 12 ~ "Childhood",
        #         age_years < 20 ~ "Adolescence",
        #         age_years < 40 ~ "Young Adulthood",
        #         age_years < 50 ~ "Middle Age",
        #         age_years < 60 ~ "Senior_50",
        #         age_years < 70 ~ "Senior_60",
        #         age_years >= 70 ~ "Senior_70",
        #         TRUE ~ NA_character_
        #       )) %>% 
        #       distinct(sample_id,  age_days, age_bin)
        #   )
        
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
        
        reference_sample
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble", "forcats", "readr", 'stringr', "cellNexus"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big_30_cores")),
      memory = "persistent",
      error = "stop"
    ),
    
    # This target loads and processes the pseudobulk sample data. It imports a HDF5 SummarizedExperiment, 
    # applies filters to retain shared genes, immune cells, and samples marked for analysis, integrates age metadata,
    # filters for common genes and samples with an appropriate number of detected genes, computes the mean library size, 
    # selects a reference sample, and performs normalisation and scaling.
    # pseudobulk_sample ------
    tar_target(
      pseudobulk_sample,
      {
        metadata <- 
          get_metadata(cache_directory="/home/zhanchen/From_scratch/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseduobulk/")
        
        metadata <- metadata %>% 
          select(
            sample_id, donor_id, dataset_id, title, collection_id, 
            age_days, sex, self_reported_ethnicity, 
            tissue, assay, cell_type_unified_ensemble, cell_type,
            disease
          ) %>% as_tibble() %>% distinct() %>% 
          edit_covariates
        
        #---------------------------------#
        # Edit or add more filters here for analyses
        #---------------------------------#
        # filter(is_gene_shared) |> 
        # filter(is_immune & do_analyse) 
        
        se <- 
          loadHDF5SummarizedExperiment(hdf5_path) %>% 
          filter(cell_type_unified_ensemble == target_cell_type) %>% 
          filter(do_analyse) %>% 
          filter(is_gene_shared) |>
          select(-c(
            age_days, sex, ethnicity_groups, tissue_groups,
            assay_groups, disease_groups, age_days_scaled
          )) 
        
        se <- se %>% 
          left_join(
            metadata %>% 
              distinct(
                sample_id, 
                age_days, age_bin, age_days_scaled,
                sex, ethnicity_groups,
                tissue, tissue_groups,
                assay_groups,
                disease_groups
              ),
            by = 'sample_id', 
            copy = T
          ) 
        
        # TEMPORARY BECAUSE I FORGOT TO INTEGRATE AGE BINS
        # se = se |>
        #   left_join(
        #     readRDS(metadata_path) |>
        #       filter(age_days > 365) |> 
        #       mutate(age_years = age_days / 365) |> 
        #       mutate(age_bin = dplyr::case_when(
        #         age_years < 3 ~ "Infancy",
        #         age_years < 12 ~ "Childhood",
        #         age_years < 20 ~ "Adolescence",
        #         age_years < 40 ~ "Young Adulthood",
        #         age_years < 50 ~ "Middle Age",
        #         age_years < 60 ~ "Senior_50",
        #         age_years < 70 ~ "Senior_60",
        #         age_years >= 70 ~ "Senior_70",
        #         TRUE ~ NA_character_
        #       )) %>% 
        #       distinct(sample_id,  age_days, age_bin)
        #   )
        
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
        # mean_library_size <- se |>
        #   assay("counts") |>
        #   _[nrow(se) |> seq_len() |> sample(size = 2000), ] |> 
        #   colSums() |>
        #   mean()
        # 
        # # Optional: retrieve the sample name (column name in the SummarizedExperiment)
        # reference_sample <- colnames(se)[
        #   se |>
        #     assay("counts") |>
        #     colSums() |>
        #     {\(x) abs(x - mean_library_size)}() |>  # Calculate absolute difference from the mean
        #     which.min()                             # Identify the smallest difference
        # ]
        
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
          # scale_abundance(method = "TMMwsp", reference_sample = reference_sample) |> 
          offset_calcuation(method = "TMMwsp", reference_sample = reference_sample) |>

          # Drop sex unknown as causes problem during fit
          mutate(
            sex = if_else(sex |> is.na(), "unknown", sex),
            ethnicity_groups = if_else(ethnicity_groups |> is.na(), "Other/Unknown", ethnicity_groups)
          ) |> 
          filter(sex != "unknown") |> 
          filter(!age_bin |> is.na()) |> 
          
          # Eliminate complete confounders
          tidybulk:::resolve_complete_confounders_of_non_interest(assay_groups, dataset_id, disease_groups) |> 
          
          # library size factor is the reciproque of the multiplier (correction factor)
          mutate(offset = log(1/multiplier)) |> 
          
          # Set intercept
          mutate(
            ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
            assay_groups___altered = fct_relevel(assay_groups___altered, "10x Genomics 3"),
            disease_groups___altered = fct_relevel(disease_groups___altered, "Normal"),
            # age_bin = fct_relevel(age_bin, "Adolescence")
            age_bin = fct_relevel(age_bin, "Senior_50")
          ) 
        
        # # Add dispersion
        # rowData(se)  = 
        #   rowData(se) |> 
        #   as_tibble(rownames = ".feature") |> 
        #   left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |> 
        #   data.frame(row.names = ".feature") |> DataFrame()
        
        se
        
        # load process data to save time when testing
        # loadHDF5SummarizedExperiment('/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseduobulk_sample_tar_load_altered/')
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble", "forcats", "readr", 'stringr', "cellNexus"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big_30_cores")),
      memory = "persistent", 
      error = "stop"
    ),
    
    # pseudobulk_sample_id ------
    # This target extracts unique sample ids from the pseudobulk sample  
    tar_target(
      pseudobulk_sample_id,
      pseudobulk_sample |> colnames(),
      packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # feature_df ------
    # This target extracts unique features from the pseudobulk sample and groups them into 
    # chunks for parallel processing.
    tar_target(
      feature_df, 
      pseudobulk_sample |> 
        distinct(.feature)|>
        # testing genes that ran for long time
        # filter(.feature %in% readRDS('/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/ning_data/ethnicity_umap_selected_genes.rds')) |>
        # slice_sample(n=1500) %>% 
        group_by(.feature) |> 
        tar_group(), 
      iteration = "group",
      packages = c( "tidySummarizedExperiment", "targets", "purrr", "dplyr"),
      resources = tar_resources(crew = tar_resources_crew("elastic"))
    ),
    
    # se_df -----
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
    
    # estimates_chunk ------
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
        
        # Manually revise data colnames to suit brms bug
        colnames(data) = colnames(data) |> stringr::str_replace_all("_+", "_")
        
        # # Check if dispersion estimation has failed
        # if(data |> pull(dispersion) |> unique() |> is.na()){
        #   warning("The dispersion calculation has failed. 1 is given as default prior.")
        #   data = data |> mutate(dispersion = 1)
        # }
        
        # Define the model formula
        formula <- bf(
          
          # Formula for counts
          counts ~ 1 + offset(offset) + age_bin*sex + disease_groups_altered + ethnicity_groups + assay_groups_altered + 
            (1 | dataset_id_altered) + 
            (1 + age_bin*sex + ethnicity_groups | tissue_groups),
          
          # Formula for dispersion
          shape ~ 1 + disease_groups_altered + assay_groups_altered + ethnicity_groups + (1 | tissue_groups)  # Model 'shape' as a function of scaled 'disp'
          
          # Using the externally, eBayes inferred overdispersion
          # shape ~ 1 + offset(log(1/dispersion))
        )
        
        # prior Version 0
        prior = c(
          prior(student_t(3, i, 1.5), class = Intercept),
          prior(student_t(3, 0, 1), class = Intercept, dpar = shape),
          prior(student_t(3, 0, 5), class = b),
          prior(student_t(3, 0, 2), class = b, dpar = shape)
          # prior(beta(0.5381488, 10.3577433), class = "zi", lb = 0, ub = 1) # addition zi from V2
        ) |>
          substitute(env = list(i = mean(log1p(data$counts / exp(data$offset))))) |>
          eval()
        # 
        # chains = 2
        # inits <- list(Intercept = mean(log1p(data$counts / exp(data$offset))))
        # inits <- replicate(chains, inits, simplify = FALSE)
        
        # test:
        # prior = prior(normal(-0.0002056948, 0.07690437))
        
        # HPC pipeline: param V1: learned from version 0
        # prior = c(
        #   prior(student_t(4.45496, 0.008599254, 1.143344), class = "b"),
        #   prior(student_t(18.16242, 0.07952513, 0.9926044), class = "b", dpar = "shape"),
        #   prior(normal(5.441626, 2.25460683), class = "Intercept"),
        #   prior(normal(0.1459487, 0.8347875), class = "Intercept", dpar = "shape"),
        #   prior(student_t(4.7655009	, 0.887529, 0.8684176), class = "sd", lb = 0),
        #   prior(student_t(53.08894, 0.9080073, 0.2870678), class = "sd", dpar = "shape", lb = 0),
        #   prior(beta(0.5541155, 9.337894), class = "zi")
        # ) 
        
        # HPC pipeline: param V2: updated from v1 and set inits
        # prior = c(
        #   prior(student_t(6.153327, 0.06161134, 0.9263627), class = "b"),
        #   prior(student_t(40.51669, 0.07603337, 0.8252114), class = "b", dpar = "shape"),
        #   prior(normal(6.057503, 2.438534), class = "Intercept"),
        #   prior(normal(0.4260793, 1.470536), class = "Intercept", dpar = "shape"),
        #   prior(student_t(52.19541	, 0.5703259, 0.4147664), class = "sd", lb = 0),
        #   prior(normal(0.8670409, 0.1779553), class = "sd", dpar = "shape", lb = 0),
        #   prior(beta(0.5381488, 10.3577433), class = "zi", lb = 0, ub = 1)
        # )
        # 
        # chains = 2
        # 
        # Kc <- 39
        # Kc_shape <- 28
        # M_1 <- 1; N_1 <- 105
        # M_2 <- 19; N_2 <- 26
        # M_3 <- 1; N_3 <- 26
        # 
        # inits <- lapply(1:chains, function(i) {
        #   list(
        #     # Fixed effects for count part
        #     b = 0.06161134 + 0.9263627 * rt(Kc, 6.153327),
        #     Intercept = rnorm(1, 6.057503, 2.438534),
        #     
        #     # Fixed effects for shape submodel
        #     b_shape = 0.07603337 + 0.8252114 * rt(Kc_shape, 40.51669),
        #     Intercept_shape = rnorm(1, 0.4260793, 1.470536),
        #     
        #     # Zero-inflation probability
        #     zi = rbeta(1, 0.5381488, 10.3577433),
        #     
        #     # Group-level standard deviations and effects
        #     sd_1 = abs(0.5703259 + 0.4147664 * rt(M_1, 52.19541)),      # count
        #     z_1 = replicate(M_1, rnorm(N_1, mean = 0 , sd = 0.08547970), simplify = FALSE),
        #     
        #     sd_2 = abs(0.5703259 + 0.4147664 * rt(M_2, 52.19541)),      # zi
        #     z_2 = matrix(rnorm(M_2 * N_2, mean = 0 , sd = 0.08547970), nrow = M_2, ncol = N_2),
        #     L_2 = diag(M_2),                                            # no correlation (identity)
        #     
        #     sd_3 = abs(rnorm(M_3, 0.8670409, 0.1779553)),               # shape
        #     z_3 = replicate(M_3, rnorm(N_3, mean = 0 , sd = 0.08547970), simplify = FALSE)
        #   )
        # })
        
        # HPC pipeline: param V3: cap df in v2 and dynamically set mu of intercept
        # prior = c(
        #   # prior(student_t(3, 0.06161134, 0.9263627), class = "b"),
        #   # prior(student_t(3, 0.07603337, 0.8252114), class = "b", dpar = "shape"),
        #   prior(normal(i, 2.438534), class = "Intercept"),
        #   prior(normal(0.4260793, 1.470536), class = "Intercept", dpar = "shape"),
        #   prior(normal(0, 5), class = b),
        #   prior(normal(0, 2), class = b, dpar = shape),
        #   # prior(student_t(3	, 0.5703259, 0.4147664), class = "sd", lb = 0),
        #   # prior(normal(0.8670409, 0.1779553), class = "sd", dpar = "shape", lb = 0),
        #   prior(beta(0.5381488, 10.3577433), class = "zi", lb = 0, ub = 1)
        #   ) |>
        #   substitute(env = list(i = mean(log1p(data$counts / exp(data$offset))))) |>
        #   eval()
        
        chains = 2

        # dynamically extract param from stan data
        # code used from brm
        bterms <- brmsterms(
          formula = brms:::validate_formula(
            formula, data = data, family = zero_inflated_negbinomial(),
            autocor = NULL, sparse = NULL, cov_ranef = NULL
          )
        )
        bframe <- brms:::brmsframe(bterms, data)
        sdata <- brms:::.standata(
          bframe, data = data, prior = prior,
          data2 = NULL, stanvars = NULL, threads = NULL
        )

        Kc <- sdata$Kc
        Kc_shape <- sdata$Kc_shape
        M_1 <- sdata$M_1; N_1 <- sdata$N_1
        M_2 <- sdata$M_2; N_2 <- sdata$N_2
        M_3 <- sdata$M_3; N_3 <- sdata$N_3

        inits <- lapply(1:chains, function(i) {
          list(

            #### revert v0 prior
            # Fixed effects for count part
            b = rnorm(Kc, 0, 5),
            # dynamically set mu for intercept
            Intercept = rnorm(1, mean(log1p(data$counts / exp(data$offset))), 1.5),

            # Fixed effects for shape submodel
            b_shape = rnorm(Kc_shape, 0, 2),
            Intercept_shape = rnorm(1, 0, 1)


            # # Fixed effects for count part
            # b = 0.06161134 + 0.9263627 * rt(Kc, 3),
            # # dynamically set mu for intercept
            # Intercept = rnorm(1, mean(log1p(data$counts / exp(data$offset))), 2.438534),
            #
            # # Fixed effects for shape submodel
            # b_shape = 0.07603337 + 0.8252114 * rt(Kc_shape, 3),
            # Intercept_shape = rnorm(1, 0.4260793, 1.470536)

            # Zero-inflation probability
            # zi = rbeta(1, 0.5381488, 10.3577433)
            #
            # # Group-level standard deviations and effects
            # sd_1 = abs(0.5703259 + 0.4147664 * rt(M_1, 3)),      # count
            # z_1 = replicate(M_1, rnorm(N_1, mean = 0 , sd = 0.08547970), simplify = FALSE),
            #
            # sd_2 = abs(0.5703259 + 0.4147664 * rt(M_2, 3)),      # zi
            # z_2 = matrix(rnorm(M_2 * N_2, mean = 0 , sd = 0.08547970), nrow = M_2, ncol = N_2),
            # L_2 = diag(M_2),                                            # no correlation (identity)
            #
            # sd_3 = abs(rnorm(M_3, 0.8670409, 0.1779553)),               # shape
            # z_3 = replicate(M_3, rnorm(N_3, mean = 0 , sd = 0.08547970), simplify = FALSE)
          )
        })
        
        brm(
          formula = formula,
          data = data,
          family = zero_inflated_negbinomial(),
          prior = prior,
          chains = chains,
          cores = pmin(as.numeric(parallelly::availableCores()), chains), 
          threads = threading(threads = (as.numeric(parallelly::availableCores()) / chains) |> floor()),
          warmup = 400, 
          refresh = 10,
          backend = "cmdstanr", 
          #sparse = TRUE,
          #save_model = glue("{external_directory}~/temp.rds"),
          #algorithm = "pathfinder",
          # sample_prior = TRUE, 
          init = inits,
          iter = 600  # Increase iterations for better convergence
        )
        
      })) |> 
        
        # Drop data because it is withn the brms object
        select(-se), 
      pattern = map(se_df),
      packages = c( "brms", "glue", "stringr", "dplyr", "purrr", "SummarizedExperiment", "tidySummarizedExperiment"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      cue = tar_cue(mode = "never"),
      error = "null"
    ),
    
    ## summary ----- 
    # This target summarises the fitted Bayesian models by performing hypothesis tests for ethnicity contrasts 
    # and extracting convergence diagnostics (Rhat) for the ethnicity parameters.
    tar_target(
      summary,
      estimates_chunk |>
    #     mutate(summary_ethnicity = map(brms_fit, ~ .x |> hypothesis(
    #       c(
    #         "Europeans" = "(ethnicity_groupsAfrican
    # + ethnicity_groupsEastAsian
    # + ethnicity_groupsHispanicDLatinAmerican
    # + ethnicity_groupsSouthAsian
    # + `ethnicity_groupsJapanese`) / 5 = 0",
    #         "EastAsian" = "(
    #    ethnicity_groupsAfrican
    #  + ethnicity_groupsHispanicDLatinAmerican
    #  + ethnicity_groupsSouthAsian
    #  + `ethnicity_groupsJapanese`
    #  - 5 * ethnicity_groupsEastAsian
    #  ) / 5 = 0",
    #         "SouthAsian" = "(
    #    ethnicity_groupsAfrican
    #  + ethnicity_groupsHispanicDLatinAmerican
    #  + ethnicity_groupsEastAsian
    #  + `ethnicity_groupsJapanese`
    #  - 5 * ethnicity_groupsSouthAsian
    #  ) / 5 = 0",
    #         "African" = "(
    #    ethnicity_groupsEastAsian
    #  + ethnicity_groupsHispanicDLatinAmerican
    #  + ethnicity_groupsSouthAsian
    #  + `ethnicity_groupsJapanese`
    #  - 5 * ethnicity_groupsAfrican
    #  ) / 5 = 0",
    #         "HispanicDLatinAmerican" = "(
    #    ethnicity_groupsAfrican
    #  + ethnicity_groupsEastAsian
    #  + ethnicity_groupsSouthAsian
    #  + `ethnicity_groupsJapanese`
    #  - 5 * ethnicity_groupsHispanicDLatinAmerican
    #  ) / 5 = 0",
            
    #         "Japanese" = "(
    #    ethnicity_groupsAfrican
    #  + ethnicity_groupsHispanicDLatinAmerican
    #  + ethnicity_groupsSouthAsian
    #  + ethnicity_groupsEastAsian
    #  - 5 * `ethnicity_groupsJapanese`
    #  ) / 5 = 0"
    #       ),
    #       #        c(
    #       #          "African" = "(ethnicity_groupsEuropean
    #       #  + ethnicity_groupsEastAsian
    #       #  + ethnicity_groupsHispanicDLatinAmerican
    #       #  + ethnicity_groupsSouthAsian
    #       #  + `ethnicity_groupsJapanese`) / 5 = 0",
    #       #          
    #       #          "Europeans" = "(
    #       #   ethnicity_groupsEastAsian
    #       # + ethnicity_groupsHispanicDLatinAmerican
    #       # + ethnicity_groupsSouthAsian
    #       # + `ethnicity_groupsJapanese`
    #       # - 5 * ethnicity_groupsEuropean
    #       # ) / 5 = 0",
    #       #          
    #       #          "EastAsian" = "(
    #       #   ethnicity_groupsEuropean
    #       # + ethnicity_groupsHispanicDLatinAmerican
    #       # + ethnicity_groupsSouthAsian
    #       # + `ethnicity_groupsJapanese`
    #       # - 5 * ethnicity_groupsEastAsian
    #       # ) / 5 = 0",
    #       #          
    #       #          "SouthAsian" = "(
    #       #   ethnicity_groupsEuropean
    #       # + ethnicity_groupsHispanicDLatinAmerican
    #       # + ethnicity_groupsEastAsian
    #       # + `ethnicity_groupsJapanese`
    #       # - 5 * ethnicity_groupsSouthAsian
    #       # ) / 5 = 0",
    #       #          
    #       #          "HispanicDLatinAmerican" = "(
    #       #   ethnicity_groupsEuropean
    #       # + ethnicity_groupsEastAsian
    #       # + ethnicity_groupsSouthAsian
    #       # + `ethnicity_groupsJapanese`
    #       # - 5 * ethnicity_groupsHispanicDLatinAmerican
    #       # ) / 5 = 0",
    #       #          
    #       #          "Japanese" = "(
    #       #   ethnicity_groupsEuropean
    #       # + ethnicity_groupsHispanicDLatinAmerican
    #       # + ethnicity_groupsSouthAsian
    #       # + ethnicity_groupsEastAsian
    #       # - 5 * `ethnicity_groupsJapanese`
    #       # ) / 5 = 0"
    #       #        ),
          
    #       # Median instead and mad of mean and sd
    #       robust=TRUE,
    #       alpha = 0.1
    #     )
    #     )) |>

         mutate(

          summary_ethnicity = map(
            
            brms_fit,  function(x) {
              
              params = x$fit %>% summary() |> _[[1]] |> rownames()
              params = params[grepl("^b_ethnicity_groups", params)] %>% sub("^b_", "", .) %>% setdiff(c("ethnicity_groupsOtherDUnknown", "ethnicity_groupsNativeAmerican&PacificIslander")) %>% paste0("`", . , "`")
              ethnicity_groups_names <- sub("`ethnicity_groups(.*)`", "\\1", params) 
              
              equations <- sapply(seq_along(params), function(i) {
                this_ethnicity <- ethnicity_groups_names[i]
                this_param <- params[i]
                other_params <- params[-i]
                avg_expr <- paste0("(", paste(other_params, collapse = " + "), ")/", length(other_params) + 1)
                eq <- paste0(this_param, " - ", avg_expr, " = 0")
                eq
              })
              names(equations) <- ethnicity_groups_names
              equations = append(
                equations, 
                c('Europeans' = paste0("(", paste(params, collapse = " + "), ")/", length(params), ' = 0'))
              )
              
              return(
                x |> hypothesis(equations, robust=TRUE, alpha = 0.1)
              )
              
            }
            
          )
          
        ) %>% 

        mutate(
          
          summary_tissue = map(
            brms_fit,  function(x) {
              
              params = x$fit %>% summary() |> _[[1]] |> rownames()
              params = params[grepl("^r_tissue_groups\\[.*?,Intercept\\]$", params)] %>% sub("^r_", "", .) %>% paste0("`", . , "`")
              tissue_names <- sub("`tissue_groups\\[(.*),Intercept\\]`", "\\1", params)
              
              equations <- sapply(seq_along(params), function(i) {
                this_tissue <- tissue_names[i]
                this_param <- params[i]
                other_params <- params[-i]
                avg_expr <- paste0("(", paste(other_params, collapse = " + "), ")/", length(other_params))
                eq <- paste0(this_param, " - ", avg_expr, " = 0")
                eq
              })
              names(equations) <- tissue_names
              
              return(
                x |> hypothesis(equations, class = "r", robust=TRUE, alpha = 0.1)
              )
              
            }
          )
          
        ) %>% 
        
        mutate(Rhat_ethnicity = map_dbl(brms_fit, 
                                        ~ summary(.x)$fixed |> 
                                          as_tibble(rownames = "par") |> 
                                          filter(par |> str_detect("ethnicity")) |> 
                                          pull(Rhat) |>
                                          max()
        )) |> 
        
        mutate(Rhat_tissue = map_dbl(brms_fit, 
                                     ~ summary(.x)$random$tissue_groups |> 
                                       as_tibble() |> 
                                       pull(Rhat) |>
                                       max()
        )) %>%  
        
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      error = "null"
    ),
    
    ## effect_removed -----
    # This target generates adjusted model estimates by removing unwanted effects from the fitted Bayesian models,
    # thereby isolating the effects of interest. Here, nuisance covariates are set to NA and removed from the predictions.
    # This target produces adjusted estimates from the Bayesian models, removing unwanted effects while retaining 
    # the tissue group random effect, thus preserving variability associated with tissue-specific factors.
    tar_target(
      effect_removed, 
      estimates_chunk |> 
        mutate(brms_fit_adjusted_ethnicity = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = TRUE, 
          re_formula = ~ 0
        ))) |> 
        mutate(brms_fit_adjusted_ethnicity_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          newdata = .x$data |> mutate(assay_groups_altered=NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = FALSE, correct_by_offset = FALSE,
          re_formula = ~ 0
        ))) |> 
        mutate(brms_fit_adjusted_tissue = map(brms_fit, ~ .x |> remove_unwanted_effect(
          newdata = .x$data |> mutate(assay_groups=NA, sex = NA, age_bin = NA, disease_groups = NA, ethnicity_groups = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = TRUE, 
          re_formula = ~ (1 | tissue_groups)
        ))) |> 
        mutate(brms_fit_adjusted_tissue_new = map(brms_fit, ~ .x |> remove_unwanted_effect_new(
          newdata = .x$data |> mutate(assay_groups_altered=NA, ethnicity_groups = NA, sex = NA, age_bin = NA, disease_groups_altered = NA, dataset_id_altered = NA), # age_bin*sex + disease_groups + ethnicity_groups + assay_groups
          robust = FALSE, correct_by_offset = FALSE,
          re_formula = ~ (1 | tissue_groups)
        ))) |> 
        mutate(brms_fit_adjusted_ethnicity_estimate = map(brms_fit_adjusted_ethnicity, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |> 
        mutate(brms_fit_adjusted_ethnicity_new_estimate = map(brms_fit_adjusted_ethnicity_new, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |> 
        mutate(brms_fit_adjusted_tissue_estimate = map(brms_fit_adjusted_tissue, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |>
        mutate(brms_fit_adjusted_tissue_new_estimate = map(brms_fit_adjusted_tissue_new, ~ {
          
          df = .x |> as_tibble()
          if (nrow(df) == length(pseudobulk_sample_id)){
            return(
              df |>
                select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |> 
                mutate(sample_id = pseudobulk_sample_id)
            )
          }else{
            return(NULL)
          }
          
        })) |>
        select(-brms_fit),
      
      pattern = map(estimates_chunk),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      error = "null"
    ),
    
    ## param -----
    tar_target(
      param,
      estimates_chunk %>%
        mutate(
          param = map(
            brms_fit, 
            ~ summary(.x$fit) |> as.data.frame()
          ) 
        )%>% 
        select(-brms_fit),
      pattern = map(estimates_chunk),
      packages = c( "brms", "dplyr", "purrr", "rstan"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      error = "null"
    ),
    
    # adjusted_matrix -----
    tar_target(
      adjusted_assay_ethnicity,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    ),
    
    tar_target(
      adjusted_assay_ethnicity_new,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_ethnicity_new_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    ),
    
    tar_target(
      adjusted_assay_tissue,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    ),
    
    tar_target(
      adjusted_assay_tissue_new,
      get_adjusted_matrix(effect_removed, brms_fit_adjusted_tissue_new_estimate),
      packages = c( "brms", "glue", "dplyr", "purrr", "rstan", "magrittr", "stringr", "tidySummarizedExperiment") ,
      resources = tar_resources(
        crew = tar_resources_crew("elastic_big_30_cores")
      )
    )
    
  ) # end ) of all target list
  
  
}, ask = FALSE, script = glue::glue(script_path))

