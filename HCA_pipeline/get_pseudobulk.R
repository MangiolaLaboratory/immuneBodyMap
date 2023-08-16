library(dplyr)
library(cellxgenedp)
library(tidyverse)
#library(tidySingleCellExperiment)
library(stringr)
library(scMerge)
library(glue)
library(DelayedArray)
library(HDF5Array)
library(openssl)
library(stringr)
library(CuratedAtlasQueryR)
library(purrr)

library(dbplyr)
library(DBI)
library(duckdb)
library(tidySingleCellExperiment)
library(tidySummarizedExperiment)
library(tidyseurat)
library(tidybulk)
library(magrittr)

root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap"
metadata_DB = glue("{root_directory}/sccomp_on_HCA_0.2.3.4/input_relative.rds")
metadata = readRDS(metadata_DB)
  

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


# Read arguments
args = commandArgs(trailingOnly = TRUE)
tissue_harmonised = args[[1]] |> str_replace_all("____", " ")
cell_type_harmonised = args[[2]] |> str_replace_all("____", " ")
is_immune = args[[3]] 
output_file = args[[4]] 

seu  = 
  metadata |>
  filter(
      tissue_harmonised == !!tissue_harmonised & 
      cell_type_harmonised == !!cell_type_harmonised &
      is_immune == !!is_immune
  ) |>
  get_single_cell_experiment(
    cache_directory = "/vast/projects/cellxgene_curated"
  ) |>
	as.Seurat(data = NULL) 

print("calculating stats")

stats = 
	seu |>
	select(sample_, nFeature_originalexp) |>
	with_groups(sample_, ~ .x |>
	summarise(
		mean_nFeature = mean(nFeature_originalexp), 
		median_nFeature = median(nFeature_originalexp)
	)) 

print("aggregating")


seu |> 
	
	# Calculate summary stats
	left_join(stats, by = "sample_") |>
	mutate(sample_ = glue("{sample_}___{cell_type_harmonised}")) |>
  tidyseurat::aggregate_cells( .sample = sample_,
  	.metadata_columns_to_keep = c('cell_type_harmonised', 'sample_', 'orig.ident', 'file_id', 'assay', 'age_days', 'development_stage', 'sex', 'ethnicity', 'tissue_harmonised', 'tissue', 'disease', 'lineage_1', 'is_immune', 'sample_count', 'age_days_original', 'ethnicity_simplified', 'assay_simplified', 'mean_nFeature', 'median_nFeature')
  ) |>
	rename(counts = originalexp) |>
	tidybulk::as_SummarizedExperiment(.sample, .feature, counts) |>
	
  saveRDS(output_file)

# x |> as.SingleCellExperiment(assay = "originalexp") |>
# 	muscat::aggregateData(
# 							assay = "originalexp", fun = "sum",
# 							by = c("sample_")
# 							)
# # one sheet per subpopulation
# assayNames(pb)