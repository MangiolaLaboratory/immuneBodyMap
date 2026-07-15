library(targets)

# pick up the current cell type from driver script
if (!exists("cur_ct", envir = .GlobalEnv)) {
  stop("`cur_ct` must be defined in the global environment before sourcing this file")
}

# ---------------------------------------------------------------------------
# Portable path resolution (publication package)
# Absolute HPC paths below are replaced with env / injected paths.
# Required:
#   PSEUDOBULK_H5AD           - harmonised pseudobulk .h5ad (Zenodo primary input)
#   IMMUNE_HEALTHY_BODY_MAP_ROOT or TAR_PIPELINE_ROOT - locate this pipeline tree
# Optional:
#   EDIT_COVARIATES_PATH      - override edit_covariates.R
#   DISEASE_GROUPING_CSV      - override disease grouping CSV
#   TARGETS_USE_CREW          - "true" (default) use Slurm crew; "false" for local
#   TARGETS_SLURM_ACCOUNT     - Slurm account (default saigencir003)
# ---------------------------------------------------------------------------
.resolve_pipeline_root <- function() {
  env <- Sys.getenv("TAR_PIPELINE_ROOT", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  root <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
  if (nzchar(root)) {
    return(normalizePath(
      file.path(root, "vignettes", "targets_pipeline"),
      winslash = "/",
      mustWork = FALSE
    ))
  }
  # When sourced by path, try relative to this file
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in seq_len(sys.nframe())) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) {
        ofile <- f
        break
      }
    }
  }
  if (!is.null(ofile)) {
    return(normalizePath(dirname(ofile), winslash = "/", mustWork = FALSE))
  }
  normalizePath(".", winslash = "/", mustWork = FALSE)
}

.pipeline_root <- .resolve_pipeline_root()
.fun_dir <- normalizePath(file.path(.pipeline_root, "functions"), winslash = "/", mustWork = FALSE)

.resolve_hdf5 <- function() {
  env <- Sys.getenv("PSEUDOBULK_H5AD", unset = "")
  if (nzchar(env) && file.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  data_root <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_DATA", unset = "")
  if (nzchar(data_root)) {
    cand <- file.path(data_root, "zenodo_release", "pseudobulk", "pseudobulk_se.h5ad")
    if (file.exists(cand)) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Set PSEUDOBULK_H5AD to the harmonised pseudobulk_se.h5ad ",
    "(Zenodo primary input), or place it under ",
    "IMMUNE_HEALTHY_BODY_MAP_DATA/zenodo_release/pseudobulk/.",
    call. = FALSE
  )
}

.resolve_edit_covariates <- function() {
  candidates <- c(
    Sys.getenv("EDIT_COVARIATES_PATH", unset = ""),
    file.path(.fun_dir, "edit_covariates.R"),
    {
      root <- Sys.getenv("IMMUNE_HEALTHY_BODY_MAP_ROOT", unset = "")
      if (nzchar(root)) {
        file.path(root, "rebuttal_CellPress", "edit_covariates.R")
      } else {
        ""
      }
    }
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Could not find edit_covariates.R. Set EDIT_COVARIATES_PATH or ",
    "IMMUNE_HEALTHY_BODY_MAP_ROOT.",
    call. = FALSE
  )
}

.hdf5_abs <- .resolve_hdf5()
.edit_cov_abs <- .resolve_edit_covariates()
.slurm_account <- Sys.getenv("TARGETS_SLURM_ACCOUNT", unset = "saigencir003")
.use_crew <- tolower(Sys.getenv("TARGETS_USE_CREW", unset = "true")) %in%
  c("1", "true", "yes", "y")

cli::cli_alert_info("targets pipeline root: {.path {.pipeline_root}}")
cli::cli_alert_info("pseudobulk h5ad: {.path {.hdf5_abs}}")
cli::cli_alert_info("edit_covariates: {.path {.edit_cov_abs}}")
cli::cli_alert_info("crew/slurm: {(.use_crew)}")


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

tar_script_expr <-
  substitute({
    tar_script({
      library(tidyverse)
      library(targets)
      library(tarchetypes)
      library(glue)
      # library(qs)
      library(crew)
      library(crew.cluster)
      
      # Set file path -----
      ## Phoenix HPC setting -----
      # hdf5_path = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseudobulk_sample_is_immune"
      # metadata_path = "/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/cell_metadata_1_0_6_sccomp_input_counts.rds"
      
      ## Input paths (injected absolute paths at script-write time) -----
      hdf5_path = HDF5_ABS
      target_cell_type = CT
      
      tar_option_set(
        memory = "transient",
        garbage_collection = 100,
        storage = "worker",
        retrieval = "worker",
        error = "continue",
        
        cue = tar_cue(mode = "never"),
        
        workspace_on_error = TRUE,
        format = "qs",
        
        # debug = "estimates_chunk",  # forces estimates_chunk to cue=always (full re-fit). Comment out for production runs.
        
        controller = if (isTRUE(USE_CREW)) crew_controller_group(
          crew_controller_slurm(
            name = "elastic_mini",
            workers = 500,
            tasks_max = 200,
            seconds_idle = 600,
            seconds_interval = 10,
            crashes_max = 7,
            options_cluster = crew_options_slurm(
              script_lines = paste0("#SBATCH -A ", SLURM_ACCOUNT),
              memory_gigabytes_required = c(20, 40, 80, 160),
              cpus_per_task = 1,
              time_minutes = c(60 * 4, 60 * 4, 60 * 4, 60 * 4, 60 * 24, 60 * 24),
              verbose = T
            )
          ),
          crew_controller_slurm(
            name = "elastic_multi_cores",
            workers = 256,
            tasks_max = 25,
            seconds_idle = 300,
            crashes_max = 7,
            options_cluster = crew_options_slurm(
              script_lines = paste0("#SBATCH -A ", SLURM_ACCOUNT),
              memory_gigabytes_required = c(20, 40, 80, 160),
              cpus_per_task = 8,
              time_minutes = c(60 * 4, 60 * 4, 60 * 4, 60 * 4, 60 * 24, 60 * 24),
              verbose = T
            )
          ),
          crew_controller_slurm(
            name = "elastic_big_mem",
            workers = 150,
            tasks_max = 20,
            seconds_idle = 30,
            crashes_max = 5,
            options_cluster = crew_options_slurm(
              script_lines = paste0("#SBATCH -A ", SLURM_ACCOUNT),
              memory_gigabytes_required = c(80, 160),
              cpus_per_task = 2,
              time_minutes = c(60 * 4, 60 * 4, 60 * 4, 60 * 4, 60 * 24, 60 * 24),
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
              script_lines = paste0("#SBATCH -A ", SLURM_ACCOUNT),
              memory_gigabytes_required = c(160, 200, 400),
              cpus_per_task = 30,
              time_minutes = c(60 * 4, 60 * 4, 60 * 4, 60 * 4, 60 * 24, 60 * 24),
              verbose = T
            )
          )
        ) else NULL
      )
      
      
      #-----------------------#
      # Functions -----
      #-----------------------#
      ### edit_covariates ------
      source(EDIT_COV_ABS)
      ### remove_unwanted_effect ------
      source(file.path(FUN_DIR, "remove_unwanted_effect.R"))
      ### get_adjusted_matrix ------
      source(file.path(FUN_DIR, "get_adjusted_matrix.R"))
      ### fit_to_age_monotonic_changes ------
      source(file.path(FUN_DIR, "fit_to_age_monotonic_changes.R"))
      ### fit_to_sex_differ_at_decades ------
      source(file.path(FUN_DIR, "fit_to_sex_differ_at_decades.R"))
      
      #-----------------------#
      # Pipeline --------
      #-----------------------#
      
      list(
        # This target loads and processes the pseudobulk sample data. It imports a HDF5 SummarizedExperiment,
        # applies filters to retain shared genes, immune cells, and samples marked for analysis, integrates age metadata,
        # filters for common genes and samples with an appropriate number of detected genes, computes the mean library size,
        # selects a reference sample, and performs normalisation and scaling.
        # pseudobulk_sample ------
        tar_target(
          pseudobulk_sample,
          {
            cli::cli_alert_info("\n Loading data...\n")
            
            se = zellkonverter::readH5AD(hdf5_path, use_hdf5 = TRUE, reader = "R") %>% filter(cell_type_unified_ensemble == target_cell_type)
            
            # Clean row names early to avoid indexing issues
            rownames(se) = se %>% rownames() %>% str_remove('_X$')
            
            tbl = se %>% colData() %>% as_tibble() %>%
              mutate(cell_type = cell_type_unified_ensemble) %>%
              edit_covariates()
            
            cli::cli_alert_info("\n Edit metadata...\n")
            
            se <- se %>%
              select(-(
                intersect(se %>% colData() %>% colnames(), tbl %>% colnames()) %>%
                  setdiff('sample_id')
              )) %>%
              inner_join(tbl, copy = T)
            
            cli::cli_alert_info(
              se %>% colData() %>% as_tibble() %>% is.na() %>% colSums() %>% {
                .[. > 0]
              } %>% names() %>% paste0(collapse  = ' ; ') %>% paste0('\nNA columns: ', .)
            )
            
            cli::cli_alert_info("\n Filter samples that have enough genes > 0 but not too many ...\n")
            # Filter common genes
            # se = se[((assay(se, "gene_presence") > 0) |> rowSums() > (ncol(se) * 0.95)),,drop=FALSE ]
            
            # Filter samples that have enough genes > 0 but not too many
            samples_with_right_number_of_detected_genes =
              (se |> assay() > 0) |>
              colSums() |>
              divide_by(nrow(se)) |>
              dplyr::between(0.1, 1)
            
            se = se[, samples_with_right_number_of_detected_genes]
            
            cli::cli_alert_info("\nCalculating reference sample for scaling gene counts...\n")
            
            # Compute mean library size
            mean_library_size <- se |>
              assay("counts") |>
              _[nrow(se) |> seq_len() |> sample(size = 2000), ] |>
              colSums() |>
              mean()
            
            # Optional: retrieve the sample name (column name in the SummarizedExperiment)
            reference_index <-
              se |>
              assay("counts") |>
              colSums() |>
              {
                \(x) abs(x - mean_library_size)
              }() |>                          # Calculate absolute difference from the mean
              which.min()                     # Identify the smallest difference
            reference_sample <- colnames(se)[reference_index]
            
            reference_path = glue::glue(
              "{targets::tar_config_get('store')}/reference_sample_{target_cell_type |> make.names()}.rds"
            )
            saveRDS(reference_sample, file = reference_path)
            cli::cli_alert_info("\nReference sample saved: {file.exists(reference_path)}")
            
            design =
              se |>
              
              # Discretise the age for the following operation
              # mutate(is_old_individual = age_days > 50 * 365) |>
              
              # This is to resolve some confounders to preserve the genes.
              # In this case we care about data variability, not the actual meaning of the variables
              resolve_complete_confounders_of_non_interest(sex, age_decade) |>
              colData() |>
              droplevels() |>
              model.matrix( ~ sex___altered + age_decade___altered, data = _)
            
            h <- hat(design)
            MinSampleSize <- 1 / max(h)
            
            counts_matrix <- assay(se, "counts")
            lib.size <- colSums(counts_matrix)
            CPM <- edgeR::cpm(counts_matrix, lib.size = lib.size)
            
            quantile_cpm = CPM %>%
              apply(
                1, quantile, 
                probs =  1 - MinSampleSize %>%
                  ceiling() %>%
                  {. / ncol(se)}
              )
            
            mini_cpm_threshold = quantile_cpm %>% quantile(0.5) %>% unname()
            cli::cli_alert_info("\nMini_cpm_threshold = {mini_cpm_threshold}")
            
            se =
              se |>
              keep_abundant(
                design = design,
                minimum_count_per_million = mini_cpm_threshold
              ) |>
              
              # Get scaling factor
              scale_abundance(method = "TMMwsp", reference_sample = reference_sample) |>
              # offset_calcuation(method = "TMMwsp", reference_sample = reference_sample) |>
              
              # Drop sex unknown as causes problem during fit
              mutate(
                sex = if_else(sex |> is.na(), "unknown", sex),
                ethnicity_groups = if_else(
                  ethnicity_groups |> is.na(),
                  "Other/Unknown",
                  ethnicity_groups
                )
              ) |>
              filter(sex != "unknown") |>
              filter(!age_bin |> is.na()) |>
              
              # Eliminate complete confounders
              tidybulk:::resolve_complete_confounders_of_non_interest(assay_groups, dataset_id, disease_groups) |>
              
              # library size factor is the reciproque of the multiplier (correction factor)
              mutate(offset = log(1 / multiplier)) |>
              
              # Set intercept
              mutate(
                ethnicity_groups = fct_relevel(ethnicity_groups, "European"),
                assay_groups___altered = fct_relevel(assay_groups___altered, "10x Genomics 3"),
                disease_groups___altered = fct_relevel(disease_groups___altered, "Normal"),
                # age_bin = fct_relevel(age_bin, "Adolescence")
                age_bin = fct_relevel(age_bin, "Senior_50"),
                age_decade = fct_relevel(age_decade, "5")
              )
            
            # # Add dispersion
            # rowData(se)  =
            #   rowData(se) |>
            #   as_tibble(rownames = ".feature") |>
            #   left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |>
            #   data.frame(row.names = ".feature") |> DataFrame()
            
            cli::cli_alert_info("\n se dimensions: {nrow(se)} x {ncol(se)}")
            
            se
            
            # load process data to save time when testing
            # loadHDF5SummarizedExperiment('/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/taskforce_shared_folder/pseduobulk_sample_tar_load_altered/')
          },
          packages = c(
            "tidybulk",
            "HDF5Array",
            "tidySummarizedExperiment",
            "magrittr",
            "tibble",
            "forcats",
            "readr",
            'stringr',
            'dplyr'
          ),
          resources = tar_resources(crew = tar_resources_crew("elastic_big_mem")),
          memory = "persistent",
          error = "stop"
        ),
        
        # sample_id ------
        # This target extracts unique sample ids from the pseudobulk sample
        tar_target(
          pseudobulk_sample_id,
          pseudobulk_sample |> colnames(),
          packages = c("tidySummarizedExperiment", "targets", "purrr", "dplyr"),
          # resources = tar_resources(crew = tar_resources_crew("elastic"))
          deployment = "main"
        ),
        
        # feature_df ------
        # This target extracts unique features from the pseudobulk sample and groups them into
        # chunks for parallel processing.
        tar_target(
          feature_df,
          pseudobulk_sample |>
            distinct(.feature) |>
            # testing genes that ran for long time
            # filter(.feature %in% readRDS('/hpcfs/groups/phoenix-hpc-mangiola_laboratory/Mangiola_ImmuneAtlas/ning_data/ethnicity_umap_selected_genes.rds')) |>
            # slice_sample(n=10) %>%
            group_by(.feature) |>
            tar_group(),
          iteration = "group",
          packages = c("tidySummarizedExperiment", "targets", "purrr", "dplyr"),
          # resources = tar_resources(crew = tar_resources_crew("elastic_mini")),
          deployment = "main"
        ),
        
        # se_df 
        # This target creates a list-column of SummarizedExperiment objects,
        # with each object corresponding to a distinct feature.
        # tar_target(
        #   se_df,
        #   feature_df |> mutate(se = map(.feature, ~
        #                                   pseudobulk_sample[.x, , drop = FALSE]))  ,
        #   pattern = map(feature_df),
        #   packages = c("brms", "glue"),
        #   # resources = tar_resources(crew = tar_resources_crew("elastic_mini"))
        #   deployment = "main"
        # ),
        
        # estimates_chunk ------
        # This target fits Bayesian models on chunks of the data. It processes each feature's data, handles missing values,
        # defines the model specification with priors, and runs the Bayesian inference using the brm function.
        tar_target(
          estimates_chunk,
          {
            warning(glue::glue("***** Gene:___{feature_df$.feature}___*****"))
            
            feature_df |> 
              mutate(se = map(.feature, ~ pseudobulk_sample[.x, , drop = FALSE])) |> 
              mutate(brms_fit = map(se, ~ {
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
                if (n_NAs > 0) {
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
                
                ### Define the model formula -----
                formula <- bf(
                  # Formula for counts
                  counts ~ 1 + offset(offset) + age_decade * sex + disease_groups_altered + ethnicity_groups + assay_groups_altered +
                    (1 | dataset_id_altered) +
                    (1 + age_decade * sex + ethnicity_groups |
                       tissue_groups),
                  
                  # Formula for dispersion
                  shape ~ 1 + disease_groups_altered + assay_groups_altered + ethnicity_groups + (1 |
                                                                                                    tissue_groups)  # Model 'shape' as a function of scaled 'disp'
                  
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
                  substitute(env = list(i = mean(log1p(
                    data$counts / exp(data$offset)
                  )))) |>
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
                    formula,
                    data = data,
                    family = zero_inflated_negbinomial(),
                    autocor = NULL,
                    sparse = NULL,
                    cov_ranef = NULL
                  )
                )
                bframe <- brms:::brmsframe(bterms, data)
                sdata <- brms:::.standata(
                  bframe,
                  data = data,
                  prior = prior,
                  data2 = NULL,
                  stanvars = NULL,
                  threads = NULL
                )
                
                Kc <- sdata$Kc
                Kc_shape <- sdata$Kc_shape
                M_1 <- sdata$M_1
                N_1 <- sdata$N_1
                M_2 <- sdata$M_2
                N_2 <- sdata$N_2
                M_3 <- sdata$M_3
                N_3 <- sdata$N_3
                
                inits <- lapply(1:chains, function(i) {
                  list(
                    #### revert v0 prior
                    # Fixed effects for count part
                    b = rnorm(Kc, 0, 5),
                    # dynamically set mu for intercept
                    Intercept = rnorm(1, mean(log1p(
                      data$counts / exp(data$offset)
                    )), 1.5),
                    
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
                
                res <- brm(
                  formula = formula,
                  data = data,
                  family = zero_inflated_negbinomial(),
                  prior = prior,
                  chains = chains,
                  cores = pmin(as.numeric(parallelly::availableCores()), chains),
                  threads = threading(
                    threads = {
                      tpc <- floor(as.numeric(parallelly::availableCores()) / chains)
                      if (tpc <= 1L) NULL else tpc
                    }
                  ),
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
                
                if (nrow(res$data) != nrow(data)) {
                  warning(
                    glue(
                      "The number of rows in the data and the number of rows in the brms object are different. The data has been filtered out."
                    )
                  )
                }
                
                return(res)
                
              })) |> 
              
              # Drop data because it is withn the brms object
              select(-se)
          },
          pattern = map(feature_df),
          packages = c(
            "brms",
            "glue",
            "stringr",
            "dplyr",
            "purrr",
            "SummarizedExperiment",
            "tidySummarizedExperiment"
          ),
          resources = tar_resources(crew = tar_resources_crew("elastic_multi_cores")),
          cue = tar_cue(mode = "never")
        ),
        
        # hypothesis_age_monotonic -----
        # tar_target(
        #   hypothesis_age_monotonic_and_adjust_tissue,
        #   estimates_chunk |>
        #     mutate(
        #       hypothesis_age_monotonic = map(brms_fit, fit_to_age_monotonic_changes)
        #     ) |>
        #     select(-brms_fit),
        #   pattern = map(estimates_chunk),
        #   packages = c("brms", "glue", "dplyr", "purrr", "rstan", "tibble", "purrr"),
        #   resources = tar_resources(crew = tar_resources_crew("elastic_mini"))
        # ),
        
        # Age and Tissue ------
        # hypothesis_age_monotonic AND Tissue_adjustment  
        
        tar_target(
          hypothesis_age_monotonic_and_adjust_tissue,
          estimates_chunk |>
            ## hypothesis_age_monotonic -----
            mutate(
              hypothesis_age_monotonic = map(brms_fit, fit_to_age_monotonic_changes)
            ) %>% 
          ## hypothesis testing ----
          # This target summarises the fitted Bayesian models by performing hypothesis tests for ethnicity/tissue contrasts
          # and extracting convergence diagnostics (Rhat) for the ethnicity parameters.
          
            mutate(summary_tissue = map(brms_fit, function(x) {
              params = x$fit %>% summary() |> _[[1]] |> rownames()
              params = params[grepl("^r_tissue_groups\\[.*?,Intercept\\]$", params)] %>% sub("^r_", "", .) %>% paste0("`", . , "`")
              tissue_names <- sub("`tissue_groups\\[(.*),Intercept\\]`", "\\1", params)
              
              equations <- sapply(seq_along(params), function(i) {
                this_tissue <- tissue_names[i]
                this_param <- params[i]
                other_params <- params[-i]
                avg_expr <- paste0("(",
                                   paste(other_params, collapse = " + "),
                                   ")/",
                                   length(other_params))
                eq <- paste0(this_param, " - ", avg_expr, " = 0")
                eq
              })
              names(equations) <- tissue_names
              
              return(x |> hypothesis(
                equations,
                class = "r",
                robust = TRUE,
                alpha = 0.1
              ))
              
            })) %>%
            
            mutate(
              Rhat_tissue = map_dbl(
                brms_fit,
                ~ summary(.x)$random$tissue_groups |>
                  as_tibble() |>
                  pull(Rhat) |>
                  max()
              )
            ) %>%
            
          ## Tissue effect_removed -----
          # This target generates adjusted model estimates by removing unwanted effects from the fitted Bayesian models,
          # thereby isolating the effects of interest. Here, nuisance covariates are set to NA and removed from the predictions.
          # This target produces adjusted estimates from the Bayesian models, removing unwanted effects while retaining
          # the tissue group random effect, thus preserving variability associated with tissue-specific factors.
          
            mutate(
              brms_fit_adjusted_tissue = map(
                brms_fit,
                ~ .x |> remove_unwanted_effect(
                  newdata = .x$data |> mutate(
                    offset = 0,
                    sex = NA,
                    age_decade = NA,
                    disease_groups_altered = NA,
                    ethnicity_groups = NA,
                    assay_groups_altered = NA
                  ),
                  # offset(offset) + age_decade * sex + disease_groups_altered + ethnicity_groups + assay_groups_altered
                  robust = TRUE,
                  re_formula = ~ (1 | tissue_groups)
                )
              )
            ) |>
            
            mutate(
              brms_fit_adjusted_tissue_estimate = map(brms_fit_adjusted_tissue, ~ {
                df = .x |> as_tibble()
                if (nrow(df) == length(pseudobulk_sample_id)) {
                  return(df |>
                           select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |>
                           mutate(sample_id = pseudobulk_sample_id))
                } else{
                  return(NULL)
                }
                
              })
            ) |>
            
            select(-brms_fit),
          
          pattern = map(estimates_chunk),
          packages = c(
            "brms",
            "glue",
            "dplyr",
            "purrr",
            "rstan",
            "magrittr",
            "stringr"
          ),
          resources = tar_resources(crew = tar_resources_crew("elastic_mini")),
          error = "null"
        ),
        
        tar_target(
          hypothesis_sex_differ_at_decades,
          estimates_chunk |>
            # hypothesis_sex_test -----
          mutate(
            hypothesis_sex_differ_at_decades = map(brms_fit, fit_to_sex_differ_at_decades)
          ) %>%
            
            select(-brms_fit),
          
          pattern = map(estimates_chunk),
          packages = c(
            "brms",
            "glue",
            "dplyr",
            "purrr",
            "rstan",
            "magrittr",
            "stringr"
          ),
          resources = tar_resources(crew = tar_resources_crew("elastic_mini")),
          error = "null"
        ),
        
        
        # adjust_age  -----
        tar_target(
          adjust_age,
          estimates_chunk |>
            
            mutate(
              brms_fit_adjusted_age = map(
                brms_fit,
                ~ .x |> remove_unwanted_effect(
                  newdata = .x$data |> mutate(
                    offset = 0,
                    sex = NA,
                    tissue_groups = NA,
                    disease_groups_altered = NA,
                    ethnicity_groups = NA,
                    assay_groups_altered = NA
                  ),
                  # offset(offset) + age_decade * sex + disease_groups_altered + ethnicity_groups + assay_groups_altered
                  robust = TRUE,
                  re_formula = NA
                )
              )
            ) |>
            
            mutate(
              brms_fit_adjusted_age_estimate = map(brms_fit_adjusted_age, ~ {
                df = .x |> as_tibble()
                if (nrow(df) == length(pseudobulk_sample_id)) {
                  return(df |>
                           select(adjusted___Estimate) |> #, adjusted___Q2.5, adjusted___Q97.5) |>
                           mutate(sample_id = pseudobulk_sample_id))
                } else{
                  return(NULL)
                }
                
              })
            ) |>
            
            select(-brms_fit),
          
          pattern = map(estimates_chunk),
          packages = c(
            "brms",
            "glue",
            "dplyr",
            "purrr",
            "rstan",
            "magrittr",
            "stringr"
          ),
          resources = tar_resources(crew = tar_resources_crew("elastic_mini")),
          error = "null"
        )
  
        
      ) # end ) of all target list
      
    }, ask = FALSE) # end of tar_script
    
  }, list(
    CT = cur_ct,
    HDF5_ABS = .hdf5_abs,
    EDIT_COV_ABS = .edit_cov_abs,
    FUN_DIR = .fun_dir,
    USE_CREW = .use_crew,
    SLURM_ACCOUNT = .slurm_account
  )) # cur_ct is the cell type; paths resolved above

eval(tar_script_expr)
cli::cli_alert_success("Write tar_script successfully!")
