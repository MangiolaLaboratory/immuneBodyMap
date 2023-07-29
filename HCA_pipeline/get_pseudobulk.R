library(tidyverse)
library(targets)
library(glue)



# Get input



result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5"

tar_script({
	
	#-----------------------#
	# Input
	#-----------------------#
	library(tidyverse)
	library(targets)
	library(tarchetypes)
	library(tidyseurat)
	library(glue)
	library(CuratedAtlasQueryR)

	#-----------------------#
	# Packages
	#-----------------------#
	tar_option_set( 
		packages = c("CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Seurat", "tidyseurat", "glue"), 
		storage = "worker", 
		retrieval = "worker", 
		error = "continue", 		
		format = "qs"
	)
	
	#-----------------------#
	# Future SLURM
	#-----------------------#
	
	library(future)
	library("future.batchtools")
	slurm <- 
		`batchtools_slurm` |>
		future::tweak( template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
									 resources=list(
									 	ncpus = 2,
									 	memory = 40000,
									 	walltime = 172800
									 )
		)
	plan(slurm)
	
	
	root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap"
	
	split_metadata = function(metadata_DB_path){
		metadata = readRDS(metadata_DB_path) |> mutate(age_days = age_days_original) 
		
		get_metadata() |>
			filter(
				sample_ %in% (
					!!metadata |>
						distinct(sample_) |> 
						tidyr::extract(sample_, "sample_", "([a-zA-Z0-9]+)_.+") |>
						pull(sample_)
				)) |>
			
			# Attach lineage
			left_join(
				read_csv("~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv") |> 
					replace_na(list(lineage_1 = "other_non_immune")) |>
					mutate(is_immune = as.character(lineage_1 == "immune")), 
				by = join_by(cell_type),
				copy = TRUE
			) |>
			filter(is_immune == "FALSE" & !is.na(lineage_1)) |>
			mutate(cell_type_harmonised = "non_immune") |>
			as_tibble() |>
			
			# Bind immune cells
			bind_rows(metadata) |>
			
			# Format covatriates
			mutate(assay = assay |> str_replace_all(" ", "_") |> str_replace_all("-", "_")  |> str_remove_all("'")) |>
			mutate(
				ethnicity = case_when(
					ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
					ethnicity |> str_detect("African") ~ "African",
					TRUE ~ ethnicity
				)
			) |>
			
			# Mutate days
			filter(development_stage!="unknown") |> 
			
			# Establish the baseline for simplified ethnicity. European as it is the most represented
			# This is so I have a tight intercept term for data simulation
			mutate(ethnicity_simplified = case_when(
				ethnicity %in% c("European", "Chinese", "African", "Hispanic or Latin American") ~ ethnicity,
				TRUE ~ "Other"
			)) |> 
			mutate(
				ethnicity_simplified = 
					ethnicity_simplified |> 
					fct_relevel(c("European", "Chinese", "African", "Hispanic or Latin American", "Other")
					)) |> 
			
			# Establish the baseline for simplified assay
			# Summarise assays to get more stable data simulations 
			# 10x as baseline
			mutate(assay_simplified = if_else(assay |> str_detect("10x"), "10x", assay)) |> 
			mutate(assay_simplified = factor(assay_simplified)) |>
			
			# Establish the baseline for disease
			mutate(disease = if_else(disease == "normal", "aaa_normal", disease)) |>
			
			# group
			group_split(tissue_harmonised, cell_type_harmonised) 
		
	}
	
	get_pseudobulk = 	function(metadata) {
			
				metadata |>
				get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated") |>
				tidySingleCellExperiment::aggregate_cells( .sample = sample_	) 
		
		}
	
	#-----------------------#
	# Pipeline
	#-----------------------#
	
	list(
		tar_target(
			metadata, 
			split_metadata(glue("{root_directory}/sccomp_on_HCA_0.2.3.4/input_relative.rds")), 
			deployment = "main",
			iteration = "list"
		),
		tar_target(
			pseudobulk, 
			get_pseudobulk(metadata), 
			pattern = map(metadata), 
			iteration = "list"
		)
	)
		

}, ask = FALSE, script = glue("{result_directory}/_targets__pseudobulk.R"))


tar_make_future(
	script = glue("{result_directory}/_targets__pseudobulk.R"),
	store = glue("{result_directory}/_targets__pseudobulk"), 
	workers = 200, 
	garbage_collection = TRUE
)


# # CREATE MAKEFILE
# output_directory = glue("{root_directory}/pseudobulk_0.2.3.4")
# output_directory |> dir.create( showWarnings = FALSE, recursive = TRUE)
# script_directory = glue("{root_directory}/HCA_pipeline")
# tab = "\t"
# metadata |>
#   distinct( tissue_harmonised, cell_type_harmonised, is_immune) |>
#   as_tibble() |>
#   mutate(
#     tissue_harmonised = tissue_harmonised |> str_replace_all(" ", "__"),
#     cell_type_harmonised = cell_type_harmonised |> str_replace_all(" ", "__")
#   ) |>
#   mutate(
#     file_name = glue("{tissue_harmonised}____{cell_type_harmonised}____{is_immune}.rds") |>
#       as.character()|>
#       str_replace_all("/", "__"),
#     output_file_path =
#       glue("{output_directory}/{file_name}") |>
#       as.character()
#   ) |>
# 	mutate(commands = pmap(
# 	  list(output_file_path, tissue_harmonised, cell_type_harmonised, is_immune),
# 		~
# 			c(
# 				glue("CATEGORY=pseudobulk\nMEMORY=50000\nCORES=1\nWALL_TIME=60000"),
# 				glue(
# 					"{..1}:{metadata_DB}\n{tab}Rscript {script_directory}/get_pseudobulk.R {..2} {..3} {..4} {..1}"
# 				)
# 			)
# 	))  |>
# 	pull(commands) |>
# 	unlist() |>
# 	write_lines(glue("{output_directory}/get_pseudobulk.makeflow"))
# 
# 
# # Read arguments
# args = commandArgs(trailingOnly = TRUE)
# tissue_harmonised = args[[1]] |> str_replace_all("__", " ")
# cell_type_harmonised = args[[2]] |> str_replace_all("__", " ")
# is_immune = args[[3]] 
# output_file = args[[4]] 
# 
# seu = 				
#   metadata |>
#   filter(
#       tissue_harmonised == !!tissue_harmonised & 
#       cell_type_harmonised == !!cell_type_harmonised &
#       is_immune == !!is_immune
#   ) |>
#   get_single_cell_experiment(
#     cache_directory = "/vast/projects/cellxgene_curated"
#   ) |>
# 	as.Seurat(data = NULL) 
# 
# print("calculating stats")
# 
# stats = 
# 	seu |>
# 	select(sample_, nFeature_originalexp) |>
# 	with_groups(sample_, ~ .x |>
# 	summarise(
# 		mean_nFeature = mean(nFeature_originalexp), 
# 		median_nFeature = median(nFeature_originalexp)
# 	)) 
# 
# print("aggregating")
# 
# 
# seu |> 
# 	
# 	# Calculate summary stats
# 	left_join(stats, by = "sample_") |>
# 	mutate(sample_ = glue("{sample_}___{cell_type_harmonised}")) |>
#   tidyseurat::aggregate_cells( .sample = sample_,
#   	.metadata_columns_to_keep = c('cell_type_harmonised', 'sample_', 'file_id', 'assay', 'age_days', 'development_stage', 'sex', 'ethnicity', 'tissue_harmonised', 'tissue', 'disease', 'lineage_1', 'is_immune', 'sample_count', 'age_days_original', 'ethnicity_simplified', 'assay_simplified', 'mean_nFeature', 'median_nFeature')
#   ) |>
# 	rename(counts = originalexp) |>
# 	tidybulk::as_SummarizedExperiment(.sample, .feature, counts) |>
# 	
#   saveRDS(output_file)
# 
# # x |> as.SingleCellExperiment(assay = "originalexp") |>
# # 	muscat::aggregateData(
# # 							assay = "originalexp", fun = "sum",
# # 							by = c("sample_")
# # 							)
# # # one sheet per subpopulation
# # assayNames(pb)