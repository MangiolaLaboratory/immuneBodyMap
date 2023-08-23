library(tidyverse)
library(targets)
library(glue)
library(CuratedAtlasQueryR)

# Get input



result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"
root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap"

# metadata = glue("{root_directory}/sccomp_on_HCA_0.2.3.4/input_relative.rds") |> readRDS() |> mutate(age_days = age_days_original) 
# 
# get_metadata() |>
# 	filter(
# 		sample_ %in% (
# 			!!metadata |>
# 				distinct(sample_) |> 
# 				tidyr::extract(sample_, "sample_", "([a-zA-Z0-9]+)_.+") |>
# 				pull(sample_)
# 		)) |>
# 	
# 	# Attach lineage
# 	left_join(
# 		read_csv("~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv") |> 
# 			replace_na(list(lineage_1 = "other_non_immune")) |>
# 			mutate(is_immune = as.character(lineage_1 == "immune")), 
# 		by = join_by(cell_type),
# 		copy = TRUE
# 	) |>
# 	filter(is_immune == "FALSE" & !is.na(lineage_1)) |>
# 	filter(cell_type_harmonised != "platelet") |>
# 	as_tibble() |>
# 	
# 	# Filter rare cell types
# 	add_count(cell_type_harmonised, name = "count_cell_type_harmonised") |>
# 	filter(count_cell_type_harmonised>100) |>
# 	select(-count_cell_type_harmonised) |>
# 	
# 	# Format covatriates
# 	mutate(assay = assay |> str_replace_all(" ", "_") |> str_replace_all("-", "_")  |> str_remove_all("'")) |>
# 	mutate(
# 		ethnicity = case_when(
# 			ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
# 			ethnicity |> str_detect("African") ~ "African",
# 			TRUE ~ ethnicity
# 		)
# 	) |>
# 	
# 	# Mutate days
# 	filter(development_stage!="unknown") |> 
# 	
# 	# Establish the baseline for simplified ethnicity. European as it is the most represented
# 	# This is so I have a tight intercept term for data simulation
# 	mutate(ethnicity_simplified = case_when(
# 		ethnicity %in% c("European", "Chinese", "African", "Hispanic or Latin American") ~ ethnicity,
# 		TRUE ~ "Other"
# 	)) |> 
# 	mutate(
# 		ethnicity_simplified = 
# 			ethnicity_simplified |> 
# 			fct_relevel(c("European", "Chinese", "African", "Hispanic or Latin American", "Other")
# 			)) |> 
# 	
# 	# Establish the baseline for simplified assay
# 	# Summarise assays to get more stable data simulations 
# 	# 10x as baseline
# 	mutate(assay_simplified = if_else(assay |> str_detect("10x"), "10x", assay)) |> 
# 	mutate(assay_simplified = factor(assay_simplified)) |>
# 	
# 	# Establish the baseline for disease
# 	mutate(disease = if_else(disease == "normal", "aaa_normal", disease)) |>
# 	
# 	# Select few columns to make things lighter
# 	select(
# 		cell_, cell_type_harmonised, sample_,  file_id, file_id_db, 
# 		age_days, development_stage, sex, tissue_harmonised, disease, 
# 		ethnicity_simplified, assay_simplified
# 	) |>
# 	# group
# 	nest(data = -c(cell_type_harmonised, tissue_harmonised)) |>
# 	unite( "name", c(cell_type_harmonised, tissue_harmonised), remove = FALSE, sep = "___") |>
# 	mutate(name = name |> str_replace_all("/", "__") |> str_replace_all(" ", "_")) |>
# 	mutate(file_path = glue("{result_directory}/{name}_metadata.qs")) |>
# 	mutate(saved = map2(
# 		data, file_path, ~ .x |> qs::qsave(.y), 
# 		.progress = TRUE
# 	)) |>
# 	select(-data) |>
# 	qs::qsave(glue("{result_directory}/pipeline_input.qs"))

tar_script({
	
	result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"
	
	
	#-----------------------#
	# Input
	#-----------------------#
	library(tidyverse)
	library(targets)
	library(tarchetypes)
	library(tidyseurat)
	library(glue)
	library(qs)
	
	
	#-----------------------#
	# Packages
	#-----------------------#
	tar_option_set( 
		packages = c(
			"CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
			"Seurat", "tidyseurat", "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR"
		), 
		storage = "worker", 
		retrieval = "worker", 
		error = "continue", 		
		format = "qs",
		#debug = "pseudobulk_df_db574b63", # Set the target you want to debug.
		cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
	)
	
	#-----------------------#
	# Future SLURM
	#-----------------------#
	
	library(future)
	library("future.batchtools")
	# slurm <- 
	# 	`batchtools_slurm` |>
	# 	future::tweak( template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
	# 								 resources=list(
	# 								 	ncpus = 20,
	# 								 	memory = 6000,
	# 								 	walltime = 172800
	# 								 )
	# 	)
	# plan(slurm)
	
	small_slurm = 
		tar_resources(
			future = tar_resources_future(
				plan = tweak(
					batchtools_slurm,
					template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
					resources = list(
						ncpus = 2,
						memory = 40000,
						walltime = 172800
					)
				)
			)
		)
	
	big_slurm = 
		tar_resources(
			future = tar_resources_future(
				plan = tweak(
					batchtools_slurm,
					template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
					resources = list(
						ncpus = 19,
						memory = 6000,
						walltime = 172800
					)
				)
			)
		)
	
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
			filter(cell_type_harmonised != "platelet") |>
			as_tibble() |>
			
			# Filter rare cell types
			add_count(cell_type_harmonised, name = "count_cell_type_harmonised") |>
			filter(count_cell_type_harmonised>100) |>
			select(-count_cell_type_harmonised) |>
			
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
			
			# Select few columns to make things lighter
			select(
				cell_, cell_type_harmonised, sample_,  file_id, file_id_db, 
				age_days, development_stage, sex, tissue_harmonised, disease, 
				ethnicity_simplified, assay_simplified
			) |>
			
			# Cell type for non immune are not summarised ernought I'm loosing a lot of samples
			filter(cell_type_harmonised != "animal_cell") |>
			mutate(cell_type_harmonised = case_when(
				cell_type_harmonised |> str_detect("fibro") ~ "stromal_cell",
				cell_type_harmonised |> str_detect("chondro") ~ "connective_tissue_cell",
				cell_type_harmonised |> str_detect("adipoc") ~ "fat_cell",
				cell_type_harmonised |> str_detect("hematopoietic") ~ "hematopoietic_cell",
				cell_type_harmonised |> str_detect("epithe") ~ "epithelial_cell",
				cell_type_harmonised |> str_detect("hepato") ~ "hepatic_cell",
				cell_type_harmonised |> str_detect("keratino") ~ "keratinocyte",
				cell_type_harmonised == "neuron" ~ "neural_cell",
				cell_type_harmonised |> str_detect("progenitor") ~ "hematopoietic_cell",
				cell_type_harmonised |> str_detect("stem") ~ "hematopoietic_cell",
				cell_type_harmonised |> str_detect("tendon") ~ "connective_tissue_cell",
				TRUE ~ cell_type_harmonised
			)) |>
			
			# group
			nest(data = -c(cell_type_harmonised, tissue_harmonised)) 
		
	}
	
	get_pseudobulk = 	function(tissue_cell_type_metadata) {
		
		tissue_cell_type_metadata |>
			mutate(data = pmap(
				list(data, cell_type_harmonised, tissue_harmonised),
				~ ..1 |>
					get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated") |>
					mutate(sample_se = 
								 	glue("{sample_}___{..2}___{..3}") |> 
								 	str_replace_all(" ", "_") |>
								 	str_replace_all("/", "__")
					) |>
					tidySingleCellExperiment::aggregate_cells( .sample = sample_se	), 
				.progress = TRUE
			)) 
		
	}
	
	samples_NOT_complete_confounders_for_ethnicity_assay = function(se){
		
		
		
		se = 
			se |> 
			# distinct(assay_simplified, ethnicity_simplified, .sample) |>
			# 
			nest(se_data = -c(assay_simplified, ethnicity_simplified)) |>
			
			# How many ethnicity per assay
			nest(data = -assay_simplified) |> 
			mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
			unnest(data) |> 
			
			# How many assay per ethnicity
			nest(data = - ethnicity_simplified) |> 
			mutate(n2 = map_int(data, ~ .x |> distinct(assay_simplified) |> nrow())) |> 
			unnest(data) 
		
		# Replace ethnicity
		dummy_assay = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(assay_simplified)
		
		se |>
			mutate(assay_simplified = if_else(n1 + n2 < 3, dummy_assay, assay_simplified)) 	|>
			
			# # Filter
			# filter(!(n1==1 & n2==1)) |>
			select(-n1, -n2) |>
			
			unnest_summarized_experiment(se_data) 
		# |>
		# 	pull(.sample) |>
		# 	unique()
	}
	
	samples_NOT_complete_confounders_for_ethnicity_disease = function(se){
		
		
		
		se = 
			se |> 
			#distinct(disease, ethnicity_simplified, .sample) |>
			
			nest(se_data = -c(disease, ethnicity_simplified)) |>
			
			# How many ethnicity per assay
			nest(data = -disease) |> 
			mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
			unnest(data) |> 
			
			# How many assay per ethnicity
			nest(data = - ethnicity_simplified) |> 
			mutate(n2 = map_int(data, ~ .x |> distinct(disease) |> nrow())) |> 
			unnest(data) 
		
		
		# Replace ethnicity
		dummy_ethnicity = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(ethnicity_simplified)
		
		se |>
			mutate(ethnicity_simplified = if_else(n1 + n2 < 3, dummy_ethnicity, ethnicity_simplified)) 	|>
			
			# # Filter
			# filter(!(n1==1 & n2==1)) |>
			select(-n1, -n2) |>
			
			unnest_summarized_experiment(se_data)
		# |>
		# 	pull(.sample) |>
		# 	unique()
	}
	
	aggregate = function(se_df){
		
		# Read
		se = 
			se_df |>
			
			mutate(data = pmap(
				list(data, cell_type_harmonised, tissue_harmonised),
				~ ..1 |>
					mutate(cell_type_harmonised = ..2, tissue_harmonised = ..3) |>
					select(-any_of(c("file_id_db", ".cell", "original_cell_id")))
			)) |>
			pull(data)
		
		# Merge
		se = do.call(cbind,  se)	
		
		tibble(tissue_harmonised, se = list(!!se))
	}
	
	se_add_dispersion = function(se_df){
		
		se_df |> 
			mutate(se = map(
				se,
				~ {
					counts = .x |> assay("counts_scaled")
					
					.x |> 
						left_join(
							
							# Dispersion data frame
							estimateDisp(counts)$tagwise.dispersion |> 
								setNames(rownames(counts)) |>
								enframe(name = ".feature", value = "dispersion")
						)
				}
			))
		
	}
	
	split_by_gene = function(se_df){
		se_df |> 
			mutate(se = map(
				se,
				~ .x |> 
					mutate(chunk___ = sample(1:10, n(), replace = TRUE)) |>
					nest(se_chunk = chunk___)
			)) |>
			unnest(se) |> 
			select(-chunk___)
	}
	
	analyse = function(se, max_rows_for_matrix_multiplication = NULL, cores = 1){
		
		# Filter
		se = 
			se |> 
			filter(sex != "unknown") |>
			filter(.aggregated_cells > 10) 
		
		# Samble to subset
		# Filter at least 5000 genes recorded
		
		# Samples with many genes
		sample_with_many_genes = 
			se |>
			assay("counts") |> 
			apply(2, function(x) (x>0) |> which() |> length()) |>
			enframe() |>
			mutate(value = as.character(value), name = as.character(name)) |>
			filter(value > 5000) |>
			pull(name)
		
		if(length(sample_with_many_genes)==0) return(se)
		
		se = 
			se |>
			filter(.sample %in% sample_with_many_genes) 
		
		
		# Eliminate complete confounders
		se = 
			se |>
			samples_NOT_complete_confounders_for_ethnicity_assay() |>
			samples_NOT_complete_confounders_for_ethnicity_disease()
		
		# Filter disease
		se = 
			se |>
			filter(disease %in% (
				se |>
					distinct(disease, sex) |>
					count(disease) |>
					filter(n>1) |>
					pull(disease)
			))
		
		# # Filter tissue that has two sexes
		# if(ncol(se)>0)
		# 	se = 
		# 	se |>
		# 	filter(tissue_harmonised %in% (
		# 		se |>
		# 			distinct(tissue_harmonised, sex) |>
		# 			count(tissue_harmonised) |>
		# 			filter(n>1) |>
		# 			pull(tissue_harmonised)
		# 	))
		# 
		# Return prematurely
		if(ncol(se) == 0) return(se)
		if(se |> distinct(sex, ethnicity_simplified) |> count(ethnicity_simplified) |> pull(n) |> max() == 1) return(se)
		
		# Build the formula
		factors = 
			c("sex", "ethnicity_simplified", "assay_simplified", "age_days", ".aggregated_cells", "disease") |>
			enframe(value = "factor") |>
			mutate(n = map_int(
				factor, ~ se |> select(.x) |> distinct() |> nrow()
			)) |>
			filter(n>1) |>
			pull(factor) |>
			str_c(collapse = " + ")
		
		random_effects = 
			c("sex", "ethnicity_simplified") |>
			enframe(value = "factor") |>
			mutate(n = map_int(
				factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
			))   |>
			filter(n>1) |>
			pull(factor) |>
			str_c(collapse = " + ")
		
		# The default
		my_formula = glue("~ {factors}")
		method = "edgeR_quasi_likelihood"
		
		if( 
			se |> distinct(cell_type_harmonised) |> nrow() > 1 &
			length(random_effects) > 0
		) {
			my_formula = glue("{my_formula} + (1 + {random_effects} | cell_type_harmonised)")
			method = "glmmseq_lme4"
		}
		
		if( 	se |> distinct(file_id) |> nrow() > 1	){
			my_formula = glue("{my_formula} + (1 | sample_)")
			method = "glmmseq_lme4"
		}
		
		# Normalise
		se = se |> quantile_normalise_abundance() 
		
		# Select abundant genes within tissues and unite
		abundant_genes = 
			se |>
			nest(data = -cell_type_harmonised) |>
			mutate(abundant_genes = map(
				data,
				~ .x |>
					keep_abundant(.abundance = counts_scaled, factor_of_interest = c(sex)) |>
					rownames(), 
				.progress = TRUE
			)) |>
			pull(abundant_genes) |>
			unlist() |>
			unique()
		
		# Vell types with enough samples
		cell_type_to_keep = 
			se |> 
			distinct(sample_, cell_type_harmonised) |> 
			count(  cell_type_harmonised) |>
			filter(n > 3) |>
			pull(cell_type_harmonised)
		
		if(length(cell_type_to_keep) == 0) return(se)
		
		se |>
			
			# Scale continuous variables
			mutate(age_days = scale(age_days) |> as.numeric()) |>
			
			# Filter abundant genes
			filter(.feature %in% abundant_genes) |>
			
			# Filter cell types to keep
			filter(cell_type_harmonised %in% cell_type_to_keep) |>
			
			# otherwise I get error for some reason
			mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised", "cell_type_harmonised")), as.character)) |>
			mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) |>
			
			# Drop random effect grouping with no enough data
			nest(data = -c(sex, ethnicity_simplified, cell_type_harmonised)) |> 
			add_count(cell_type_harmonised) |> 
			filter(n>1) |> 
			unnest(data) |>
			
			# Test
			test_differential_abundance(
				as.formula(my_formula),
				.abundance = counts_scaled,
				method = method,
				cores = cores, 
				max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
				dispersion = dispersion
			)
	}
	
	#-----------------------#
	# Pipeline
	#-----------------------#
	list(
		
		# Do metadata
		tarchetypes::tar_group_by(
			tissue_cell_type_metadata, 
			split_metadata(glue("{root_directory}/sccomp_on_HCA_0.2.3.4/input_relative.rds")), 
			cell_type_harmonised, tissue_harmonised,
			deployment = "main"
		),
		
		# Get pseudobulk
		tar_target(	
			pseudobulk_df, 
			get_pseudobulk(tissue_cell_type_metadata), 
			pattern = map(tissue_cell_type_metadata),
			resources = small_slurm
		),
		
		# tissue analysis
		tarchetypes::tar_group_by(pseudobulk_df_tissue, pseudobulk_df, tissue_harmonised),
		
		# Group samples
		tar_target(
			pseudobulk_df_tissue_merged, 
			pseudobulk_df_tissue |> aggregate() |> se_add_dispersion(), 
			pattern = map(pseudobulk_df_tissue), 
			resources = small_slurm
		),
		
		# Split in gene chunks
		tar_target(
			pseudobulk_df_tissue_split_by_gene, 
			pseudobulk_df_tissue_merged |> split_by_gene(), 
			pattern = map(pseudobulk_df_tissue_merged), 
			resources = small_slurm
		),
		
		# Analyse
		tar_target(
			estimates, 
			pseudobulk_df_tissue_split_by_gene |> mutate(se = map(
				se,
				~ .x |> analyse(max_rows_for_matrix_multiplication = 10000, cores = 18)
			)), 
			pattern = map(pseudobulk_df_tissue_split_by_gene), 
			resources = big_slurm
		)
	)
	
	
}, ask = FALSE, script = glue("{result_directory}/_targets__pseudobulk_non_immune_split.R"))

job::job({
	result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"
	
	tar_make_future(
		script = glue("{result_directory}/_targets__pseudobulk_non_immune_split.R"),
		store = glue("{result_directory}/_targets__pseudobulk_non_immune_split"), 
		workers = 200, 
		garbage_collection = TRUE
	)
})

tar_read(pseudobulk_df_tissue, store = glue("{result_directory}/_targets__pseudobulk_non_immune_split"), branches = 1)

tar_visnetwork(
	store = glue("{result_directory}/_targets__pseudobulk_non_immune_split")
)

