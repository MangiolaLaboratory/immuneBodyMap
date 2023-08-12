

library(tidyverse)
library(forcats)
library(CuratedAtlasQueryR)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
library(forcats)
library(stringr)

source(
  "https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel"
)

# CURATION OF THE DATASETS
# # Dataset curated for enrichment
# input_absolute |> 
# 	distinct(file_id, tissue_harmonised, tissue) |> 
# 	nest(data = -tissue_harmonised) |> 
# 	mutate(n = map_int(data, nrow)) |>
# 	arrange(n) |> 
# 	#filter(!(n<3 | tissue_harmonised |> str_detect("intestine"))) |> 
# 	unnest(data) |> 
# 	left_join(get_metadata() |> distinct(file_id, collection_id), copy=T) |> 
# 	mutate(collection_url = glue("https://cellxgene.cziscience.com/collections/{collection_id}")) |> 
# 	write_csv("dev/datasets_to_check_in_the_literature_rest.csv")
# 


# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_1 = args[[1]]
lineage_table = args[[2]]
output_common = args[[3]]
output_absolute = args[[4]]
output_relative = args[[5]]

# DATA INPUT
common_data =

  # get_metadata() |>
  readRDS( input_1) |>
	
	# Fix bug
	mutate(tissue_harmonised = tissue_harmonised |> str_replace("vasculaturew", "vasculature")) |>
	mutate(tissue_harmonised = tissue_harmonised |> str_replace("plcenta", "placenta")) |>
	
	
	# I HAVE TO INTEGRATE INTO PIPELINE
	# Correct fishy stem cell labelling
	# If stem for the study's annotation and blueprint is non-immune it is probably wrong, 
	# even because the heart has too many progenitor/stem
	mutate(confidence_class = case_when(
		cell_type_harmonised == "stem" & cell_annotation_blueprint_singler %in% c(
			"skeletal muscle", "adipocytes", "epithelial", "smooth muscle", "chondrocytes", "endothelial"
		) ~ 5,
		TRUE ~ confidence_class
	)) |> 
	
	# Filter at all
	filter(!(cell_type_harmonised == "stem" & cell_annotation_blueprint_singler %in% c(
		"skeletal muscle", "adipocytes", "epithelial", "smooth muscle", "chondrocytes", "endothelial"
	))) |> 

	# # Filter only primary data
	# filter(is_primary_data_x=="TRUE") |>
	
	# Filter so we can add disease to the model
	nest(data = -tissue_harmonised) |>
	mutate(has_normal = map_int(data, ~ .x |> filter(disease == "normal") |> nrow()) > 0) |>
	mutate(disease_count = map_int(data, ~ .x |> distinct(disease) |> nrow())) |>
	filter(!(disease_count==1 & !has_normal)) |>
	unnest(data) |>
	# mutate(tissue_harmonised = tissue_harmonised |> fct_relevel(tissue_harmonised!="normal")) |>
	
  dplyr::select(
    cell_,
    cell_type,
    cell_type_harmonised,
    file_id,
    file_id_db,
    assay,
    age_days,
    development_stage,
    sex,
    ethnicity,
    confidence_class,
    tissue_harmonised,
    tissue,
    sample_,
    disease
  ) |>
  as_tibble() |>

	# Attach lineage
	left_join(read_csv(lineage_table) |> replace_na(list(lineage_1 = "other_non_immune"))) |>
	mutate(is_immune = lineage_1 == "immune") |>
	
	# Fix typo
	mutate(tissue_harmonised = tissue_harmonised |> str_replace("plcenta", "placenta")) |>

  # Fix hematopoietic misclassification
  mutate(
    cell_type_harmonised = if_else(
      !is_immune &
        cell_type |> str_detect("hematopoietic"),
      "stem",
      cell_type_harmonised
    )
  ) |>

  # Filter intestine as it does not fit the small vs large paradigm
  filter(tissue != "intestine") |>
  filter(tissue != "appendix") |>

  # Filter out
  filter(!cell_type |> str_detect("erythrocyte")) |>
	filter(!cell_type |> str_detect("erythroblast")) |>
	filter(!cell_type |> str_detect("erythroid lineage cell")) |>
	filter(!cell_type |> str_detect("erythroid progenitor cell")) |>
	filter(!cell_type |> str_detect("Langerhans cell")) |>
	filter(!cell_type |> str_detect("platelet")) |>
	
  # Filter out samples with less than 30 cells
  add_count(sample_, name = "sample_count") |> 
  filter(sample_count>30) |> 
  
  # Format covatriates
  mutate(assay = assay |> str_replace_all(" ", "_") |> str_replace_all("-", "_")  |> str_remove_all("'")) |>
  mutate(
    ethnicity = case_when(
      ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
      ethnicity |> str_detect("African") ~ "African",
      TRUE ~ ethnicity
    )
  ) |>

  # Fix samples with multiple assays
  unite("sample_", c(sample_ , assay, disease), remove = FALSE) |>

  # Scale age
	mutate(age_days_original = age_days) |>
  mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>

	filter(!is.na(tissue_harmonised)) |>

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
	mutate(disease = disease |> str_replace("normal", "aaa_normal"))
	
	
	

# Save
common_data |> saveRDS(output_common)

# - Immune proportion per tissue
common_data |>

	# Fix typo
	mutate(tissue_harmonised = tissue_harmonised |> str_replace("plcenta", "placenta")) |>

	# Eliminate the samples with immune enrichment/depletion
	anti_join(
		read_csv("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/dev/datasets_to_check_in_the_literature_checked.csv") |> 
			filter(overall_immune_cells_enriched>0) |> 
			distinct(file_id, tissue_harmonised)
	) |> 
	
  # Filter unrepresented organs
  filter(!tissue_harmonised %in% c(
  	"rectum",
  	"thyroid gland",
  	"salival_gland",
  	"testis",
  	"skeletal_muscle",
  	"saliva",
  	"lacrimal gland",
  	"appendix",

  	# Lymphoid
  	"blood",
  	"lymph node",
  	"bone",
  	"spleen",
  	"thymus"
  )) |>

  # Filter Immune enriched dataset
  filter(file_id != "e756c34a-abe7-4822-9a35-55ef12270247") |>
  filter(file_id != "ca4a7d56-739b-4e3c-8ecd-28704914cc14") |>
  filter(
    file_id != "59dfc135-19c1-4380-a9e8-958908273756" |
      !tissue_harmonised %in% c("intestine small", "intestine large")
  ) |>
  
  # Filter extremes
  nest(data = -c(sample_, tissue_harmonised)) |>
  mutate(n_immune = map_int(data, ~ .x |> filter(is_immune) |> nrow())) |> 
  mutate(n__NON_immune = map_int(data, ~ .x |> filter(!is_immune) |> nrow())) |> 
  
  # Filter samples which not include non immune cells
  filter(n_immune > 0) |> 
  
  # Filter samples which include > 90% of immune cells
  filter((n_immune / (n_immune + n__NON_immune )) < 0.75) |> 
  
  unnest(data) |>
  
	# # Filter samples that include too many immune cells, as they are solid tissues
	# nest(data = -sample_) |>
	# mutate(
	# 	immune = map_int(data, ~.x |> filter(is_immune=="TRUE") |> nrow()),
	# 	non_immune = map_int(data, ~.x |> filter(is_immune=="FALSE") |> nrow()),
	# 	immune_proportion = immune / (immune + non_immune)
	# ) |>
	# filter(immune_proportion |> between(0, 0.9)) |>
	# unnest(data) |>

	mutate(is_immune = as.character(is_immune)) |>

  saveRDS(output_absolute)


# Relative
relative_data =
  common_data |>

  # Filter low confidence
  filter(confidence_class %in% c(1, 2, 3)) |>

  # Filter unrepresented organs
  filter(!tissue_harmonised %in% c("rectum", "thyroid gland", "salival_gland")) |>

  # Filter only immune
  filter(is_immune) |>
	mutate(is_immune = as.character(is_immune)) |>
  saveRDS(output_relative)



