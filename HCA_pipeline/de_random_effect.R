
library(tidySummarizedExperiment)
library(HPCell)
library(magrittr)
library(tibble)
devtools::load_all("~/PostDoc/tidybulk/")



# Dispersion 2 days calculation
job::job({
  library(tictoc)
  tic()
  glmGamPoi_subsample = 
    pseudobulk_sample |> 
    glmGamPoi::glm_gp(
      on_disk = T,
      subsample = TRUE,
      design = 
        pseudobulk_sample |> 
        tidybulk::resolve_complete_confounders_of_non_interest(dataset_id, assay_groups, tissue_groups, disease_groups) |> 
        colData() |> 
        droplevels() |> 
        model.matrix(~ dataset_id + assay_groups + disease_groups, data = _  ) , 
      verbose = TRUE,
      use_assay = "counts"
    ) 
  toc()
})


# overdispersion with NO subsampling ~ 6K samples
glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_3.rds")$overdispersions
glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])


pseudobulk_sample_small_mat = pseudobulk_sample_small
as = assay(pseudobulk_sample_small_mat) |> as.matrix() |> apply(2, as.integer)
rownames(as) = rownames(pseudobulk_sample_small_mat)
assay(pseudobulk_sample_small_mat) = as

pseudobulk_sample_small_mat <- 
  pseudobulk_sample_small_mat |> 
  DESeqDataSet(design = ~ sex) |> 
  estimateSizeFactors() |> 
  estimateDispersions()

design = model.matrix(~ sex, data = pseudobulk_sample_small |> colData() |> droplevels())

x = glmGamPoi::glm_gp(pseudobulk_sample_small, on_disk = T, design = design)$overdispersion_shrinkage_list$ql_disp_estimate
y = rowData(pseudobulk_sample_small_mat)$dispGeneEst
z = pseudobulk_sample_small |> 
  assay() |> 
  as.matrix() |> 
  edgeR::estimateDisp(tagwise = TRUE, robust = TRUE, design = design) %$% 
  tagwise.dispersion

design = model.matrix(~ condition, data = dds |> colData() |> droplevels())

library(tictoc)
library(DESeq2)
library(edgeR)
library(glmGamPoi)
pseudobulk_sample_small_mat <- 
  tidybulk::breast_tcga_mini_SE |> 
  select(-count_scaled) |> 
  DESeqDataSet(design = ~ Call) |> 
  estimateSizeFactors() |> 
  estimateDispersions(fitType="glmGamPoi")

design = model.matrix(~ Call, data = tidybulk::breast_tcga_mini_SE |> colData() |> droplevels())

microbenchmark::microbenchmark(list = list(
  glmGamPoi = 
    tidybulk::breast_tcga_mini_SE |> 
    glmGamPoi::glm_gp(design = design, use_assay = "count", ) %$% 
    overdispersion_shrinkage_list %$% 
    ql_disp_estimate,
  glmGamPoi_disk = 
    tidybulk::breast_tcga_mini_SE |> 
    glmGamPoi::glm_gp(on_disk = T, design = design, use_assay = "count") %$% 
    overdispersion_shrinkage_list %$% 
    ql_disp_estimate,
  DESeq2 = rowData(pseudobulk_sample_small_mat)$dispGeneEst,
  egdeR = tidybulk::breast_tcga_mini_SE |> 
    mutate(counts = count) |> 
    SE2DGEList() |> 
    calcNormFactors() |> 
    edgeR::estimateDisp(tagwise = TRUE,   design = design) %$% 
    tagwise.dispersion
))

glmGamPoi = 
  tidybulk::breast_tcga_mini_SE |> 
  glmGamPoi::glm_gp(design = design, use_assay = "count") %$% overdispersions
glmGamPoi_disk = 
  tidybulk::breast_tcga_mini_SE |> 
  glmGamPoi::glm_gp(on_disk = T, design = design, use_assay = "count") %$% overdispersions
DESeq2 = rowData(pseudobulk_sample_small_mat)$dispGeneEst
egdeR = tidybulk::breast_tcga_mini_SE |> 
  mutate(counts = count) |> 
  SE2DGEList() |> 
  calcNormFactors() |> 
  edgeR::estimateDisp(tagwise = TRUE,   design = design) %$% 
  tagwise.dispersion

pairs(data.frame(glmGamPoi_disk, glmGamPoi, DESeq2, egdeR))



result_directory = "/vast/scratch/users/mangiola.s/DE_pseudobulk_sample_cellNexus_1_0_3"
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
          memory_gigabytes_required = c(10, 20, 40, 80, 160), 
          cpus_per_task = 2, 
          time_minutes = c(60*4, 60*4, 60*4, 60*24, 60*24),
          verbose = T
        )
      ),
      crew_controller_slurm(
        name = "elastic_big",
        workers = 50,
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
  # Pipeline
  #-----------------------#
  list(
    
    
    # tar_target(
    #   result_directory,
    #   "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/de_ethnicity_pseudobulk_sample"
    # ),
    tar_target(
      glmGamPoi_overdispersions,
      {
        glmGamPoi_overdispersions  = readRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/glmGamPoi_all_samples_no_subsampling_cellNexus_1_0_3.rds")$overdispersions
        glmGamPoi_overdispersions[glmGamPoi_overdispersions>1e5] = max(glmGamPoi_overdispersions[glmGamPoi_overdispersions<1e5])
        glmGamPoi_overdispersions
      }, 
      deployment = "main"
      
    ),
    
    tar_target(
      pseudobulk_sample,
      {
        se = 
          loadHDF5SummarizedExperiment("/vast/projects/cellxgene_curated/cellNexus/pseudobulk_sample_is_immune") |> 
          filter(is_gene_shared) |> 
          filter(is_immune & do_analyse) 
        
        samples_with_right_number_of_detected_genes = 
          (se |> assay() > 0) |> 
          colSums() |> 
          divide_by(nrow(se)) |> 
          between(0.3, 0.7)
        
        se = se[,samples_with_right_number_of_detected_genes] 
        
        se = 
          se |> 
          keep_abundant(design = 
                          se |> 
                          mutate(is_old_individual = age_days > 50*365) |> 
                          resolve_complete_confounders_of_non_interest(tissue_groups, sex, ethnicity_groups, is_old_individual) |> 
                          colData() |> 
                          droplevels() |> 
                          model.matrix(~ tissue_groups + sex + ethnicity_groups + is_old_individual, data = _  ), 
                        minimum_counts = 100
          ) |> 
          
          # Get scaling factor
          scale_abundance(method = "TMMwsp", reference_sample = "0edf00b9d5cd39b046f90be198fb07db___1") |> 
          
          # Drop sex unknown as causes problem during fit
          filter(sex != "unknown") |> 
          
          # Eliminate complete confounders
          tidybulk:::resolve_complete_confounders_of_non_interest(assay_groups, dataset_id, disease_groups) |> 
          
          mutate(offset = log(multiplier))
        
        # Add dispersion
        rowData(se)  = 
          rowData(se) |> 
          as_tibble(rownames = ".feature") |> 
          left_join(glmGamPoi_overdispersions |> enframe(name = ".feature", value = "dispersion")) |> 
          mutate(dispersion = 1 / dispersion) |> 
          data.frame(row.names = ".feature") |> DataFrame()
          
        se
        
      }, 
      packages = c("tidybulk", "HDF5Array", "tidySummarizedExperiment", "magrittr", "tibble"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big"))
      
    ),
    
    # Split in gene chunks
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
    tar_target(
      se_df, 
        feature_df |> mutate(se = map(.feature, ~ 
          pseudobulk_sample[.x, , drop=FALSE]
        ))  , 
      pattern = map(feature_df),
      packages = c( "brms", "glue"),
      resources = tar_resources(crew = tar_resources_crew("elastic_big"))
      
    ),
    tar_target(
      estimates_chunk, 
        
        se_df |> mutate(brms_fit = map(se, ~ {
          
          data = 
            .x |>
            as_tibble() |> 
            mutate(counts = counts |> as.integer()) |> 
            droplevels()
          
          # Define the model formula
          formula <- bf(
            counts ~ 1 + offset(offset) + age_bin_sex_specific*sex + disease_groups + ethnicity_groups + assay_groups + 
              (1 | dataset_id) + 
              (1 + age_bin_sex_specific*sex + ethnicity_groups | tissue_groups),
            shape ~ 1 + dispersion  # Model 'shape' as a function of scaled 'disp'
          )
          
          prior = c(
            prior(normal(i, 2), class = Intercept),
            prior(normal(0, 2), class = Intercept, dpar = shape),
            prior(normal(0, 2), class = b, dpar = shape)
          ) |> 
            substitute(env = list(i = mean(log1p(data$counts * exp(data$offset))))) |> 
            eval()
          
          chains = 2
          inits <- list(Intercept = mean(log1p(data$counts * exp(data$offset))))
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
          
          
          
          }))  , 
      pattern = map(se_df),
      packages = c( "brms", "glue", "dplyr", "purrr", "SummarizedExperiment", "tidySummarizedExperiment"),
      resources = tar_resources(crew = tar_resources_crew("elastic")),
      cue = tar_cue(mode = "never")
      
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

x = tar_read_raw("estimates_chunk_062dd2621e3cb7e1", store = glue::glue("{result_directory}/_targets"), branches = 1)
x |> saveRDS("/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_chunk_062dd2621e3cb7e1.rds")
system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/estimates_chunk_062dd2621e3cb7e1.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/removal_unwanted_effects/")


tar_workspace(
  estimates_chunk_f3b204e35ac2d218, 
  store = glue::glue("{result_directory}/_targets"),
  script = glue::glue("{result_directory}/_targets.R")
)



tar_meta(store = glue::glue("{result_directory}/_targets")) |> 
  arrange(desc(time)) |>
  filter(!error |> is.na()) |> 
  select(name, error)
