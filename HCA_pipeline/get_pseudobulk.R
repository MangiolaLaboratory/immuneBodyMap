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

root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/immuneHealthyBodyMap"
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
#     tissue_harmonised = tissue_harmonised |> str_replace_all(" ", "____"),
#     cell_type_harmonised = cell_type_harmonised |> str_replace_all(" ", "____")
#   ) |>
#   mutate(
#     file_name = glue("{tissue_harmonised}__{cell_type_harmonised}__{is_immune}.rds") |>
#       as.character()|>
#       str_replace_all("/", "____"),
#     output_file_path =
#       glue("{output_directory}/{file_name}") |>
#       as.character()
#   ) |>
# 	mutate(commands = pmap(
# 	  list(output_file_path, tissue_harmonised, cell_type_harmonised, is_immune),
# 		~
# 			c(
# 				glue("CATEGORY=pseudobulk\nMEMORY=50000\nCORES=1\nWALL_TIME=30000"),
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

metadata  = 
  metadata |>
  filter(
      tissue_harmonised == !!tissue_harmonised & 
      cell_type_harmonised == !!cell_type_harmonised &
      is_immune == !!is_immune
  ) |>
  get_single_cell_experiment(
    cache_directory = "/vast/projects/cellxgene_curated"
  ) |>
	as.Seurat(data = NULL) |> 
	
	# Calculate summary stats
	nest(data = -sample_) |> 
	mutate(
		mean_nFeature = map_dbl(data, ~ .x |> pull(nFeature_originalexp) |> mean()), 
		median_nFeature = map_dbl(data, ~ .x |> pull(nFeature_originalexp) |> median())
	) |> 
	unnest(data) |>
  tidyseurat::aggregate_cells(c(sample_, cell_type_harmonised)) |>
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