library(tidyverse)
library(targets)
library(glue)
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
	library(glue)
	library(qs)
	library(crew)
	library(crew.cluster)
	
	#-----------------------#
	# Packages
	#-----------------------#
	tar_option_set(
		packages = c(
			"CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
			"Seurat", "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR", "crew", "magrittr", "digest", "glmmSeq", "readr"
		),
		
		garbage_collection = TRUE,
		#trust_object_timestamps = TRUE,
		memory = "transient",
		storage = "worker",
		retrieval = "worker",
		error = "continue",
		format = "qs",
		
		#-----------------------#
		# SLURM
		#-----------------------#
		controller = crew_controller_group(
			
			crew_controller_slurm(
				name = "slurm_2_20",
				slurm_memory_gigabytes_per_cpu = 20,
				slurm_cpus_per_task = 2,
				workers = 100,
				verbose = T, 
				seconds_timeout = 30
				#,
				#script_lines = "module load R/4.2.1",
				#host = "spartan.hpc.unimelb.edu.au"
			),
			
			crew_controller_slurm(
				name = "slurm_1_40",
				slurm_memory_gigabytes_per_cpu = 40,
				slurm_cpus_per_task = 1,
				workers = 200,
				verbose = T
			),
			
			crew_controller_slurm(
				name = "slurm_1_80",
				slurm_memory_gigabytes_per_cpu = 80,
				slurm_cpus_per_task = 1,
				workers = 20,
				verbose = T
			),
			
			crew_controller_slurm(
				name = "slurm_1_120",
				slurm_memory_gigabytes_per_cpu = 120,
				slurm_cpus_per_task = 1,
				workers = 10,
				verbose = T
			),
			
			crew_controller_slurm(
				name = "slurm_2_400",
				slurm_memory_gigabytes_per_cpu = 400,
				slurm_cpus_per_task = 2,
				workers = 5,
				verbose = T
			),
			
			crew_controller_slurm(
				name = "slurm_3",
				slurm_memory_gigabytes_per_cpu = 10,
				slurm_cpus_per_task = 3,
				workers = 200,
				verbose = T
			),
			
			crew_controller_slurm(
				name = "big_slurm",
				slurm_memory_gigabytes_per_cpu = 4,
				slurm_cpus_per_task = 5,
				workers = 300,
				verbose = T
			),
			
			crew_controller_local(name = "local_30", workers = 30)
		),
		resources = tar_resources(crew = tar_resources_crew("slurm_2_20")) ,
		debug = "pseudobulk_df", # Set the target you want to debug.
		cue = tar_cue(mode = "never")		
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
			#filter(is_immune == "FALSE" & !is.na(lineage_1)) |>
			filter(!is.na(lineage_1)) |>
			filter(!cell_type_harmonised %in% c("platelet", "immune_unclassified")) |>
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
				ethnicity_simplified, assay_simplified, is_immune
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
			
			# Split very big datasets with a lot of samples, with maximum 100 samples
			nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_)) |> 
			with_groups(
				c(cell_type_harmonised, tissue_harmonised, file_id, is_immune), 
				~ .x |> 
					arrange(file_id) |> 
					mutate(sample_chunk = row_number() |>
								 	
								 	# MAX SAMPLES
								 	divide_by(100) |> 
								 	ceiling()
					)
			) |> 
			unnest(data) |> 
			
			# group
			nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_chunk))
		
	}
	
	get_sce = 	function(tissue_cell_type_metadata) {
		
		tissue_cell_type_metadata |>
			mutate(data = pmap(
				list(data, cell_type_harmonised, tissue_harmonised),
				~ ..1 |>
					get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated") |>
					mutate(sample_se =
								 	
								 	# I need to fix Curated CellAtlas with disease sample, duplication for 
								 	# file_id=="cc3ff54f-7587-49ea-b197-1515b6d98c4c", cell_type_harmonised=="stromal_cell"
								 	# for lung
								 	glue("{sample_}___{disease}___{..2}___{..3}") |>
								 	str_replace_all(" ", "_") |>
								 	str_replace_all("/", "__")
					) 
			))
		
	}
	
	get_pseudobulk = 	function(sce_df) {
		
		sce_df |>
			mutate(data = map(
				data, 
				~ {
					.x = tidySingleCellExperiment::aggregate_cells(.x, .sample = sample_se	)
					
					assay(.x, "counts") = assay(.x, "counts") |> as("sparseMatrix")
					
					.x
				}
			))
		
	}
	
	nest_detect_complete_confounder = function(.data, .col1, .col2){
		
		.col1 = enquo(.col1)
		.col2 = enquo(.col2)
		
		.data |>
			
			nest(se_data = -c(!!.col1, !!.col2)) |>
			
			# How many ethnicity per assay
			nest(data = -!!.col1) |>
			mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
			unnest(data) |>
			
			# How many assay per ethnicity
			nest(data = - !!.col2) |>
			mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
			unnest(data)
	}
	
	left_join_detect_complete_confounder = function(.data, .col1, .col2){
		
		.col1 = enquo(.col1)
		.col2 = enquo(.col2)
		
		counts = 
			.data |>
			
			
			# nest(se_data = -c(!!.col1, !!.col2)) |>
			distinct(!!.col1, !!.col2) |> # Quicker
			
			# How many ethnicity per assay
			nest(data = -!!.col1) |>
			mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
			unnest(data) |>
			
			# How many assay per ethnicity
			nest(data = - !!.col2) |>
			mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
			unnest(data)
		
		.data |> 
			left_join(counts) #, by = join_by(!!.col1, !!.col2))
	}
	
	drop_samples_complete_confounder = function(.data, .col1, .col2){
		
		.col1 = enquo(.col1)
		.col2 = enquo(.col2)
		
		.data |> 
			left_join_detect_complete_confounder(!!.col1, !!.col2) |> 
			filter(n1 + n2 > 2) |> 
			select(-n1, -n2) 
		
	}
	

	aggregate = function(se_df){
		
		print("Start aggregate")
		gc()
		
		se_df |>
			
			# Add columns and filter
			mutate(data = pmap(
				list(data, cell_type_harmonised, tissue_harmonised),
				~ {
					# Add columns
					se = 
						..1 |>
						mutate(cell_type_harmonised = ..2, tissue_harmonised = ..3) |>
						select(-any_of(c("file_id_db", ".cell", "original_cell_id")))
					
					
					# Identify samples with many genes
					sample_with_many_genes =
						se |>
						assay("counts") |>
						apply(2, function(x) (x>0) |> which() |> length()) |>
						enframe() |>
						mutate(value = as.character(value), name = as.character(name)) |>
						filter(value > 5000) |>
						pull(name)
					se = se[,sample_with_many_genes, drop=FALSE]
					
					# Filter samples with too few cells
					se = se |> filter(.aggregated_cells > 10)
					
				}, 
				.progress=TRUE
			)) |>
			
			# nest for output
			nest(data = -tissue_harmonised) |> 
			
			# bind
			mutate(data = map(
				data,
				~ {
					
					se = do.call(cbind,  .x |> pull(data))
					
					# Filter very rare gene-transcripts
					all_zeros = assay(se, "counts") |> rowSums() |> equals(0) 
					se = se[!all_zeros,]
					lower_than_total_counts = assay(se, "counts") |> rowSums() < 15 
					se = se[!lower_than_total_counts,]
					
					# Make it sparse
					se@assays@data$counts = as(se@assays@data$counts, "sparseMatrix")
					
					se
				}
			))
		
		
	}
	
	aggregate_cell_type = function(se_df){
		
		print("Start aggregate")
		gc()
		
		se_df |>
			
			# Add columns and filter
			mutate(data = pmap(
				list(data, cell_type_harmonised, tissue_harmonised, file_id),
				~ {
					# Add columns
					se = 
						..1 |>
						mutate(cell_type_harmonised = ..2, tissue_harmonised = ..3, file_id = ..4) |>
						select(-any_of(c("file_id_db", ".cell", "original_cell_id")))
					
					
					# Identify samples with many genes
					sample_with_many_genes =
						se |>
						assay("counts") |>
						apply(2, function(x) (x>0) |> which() |> length()) |>
						enframe() |>
						mutate(value = as.character(value), name = as.character(name)) |>
						filter(value > 5000) |>
						pull(name)
					se = se[,sample_with_many_genes, drop=FALSE]
					
					# Filter samples with too few cells
					se = se |> filter(.aggregated_cells > 10)
					
				}, 
				.progress=TRUE
			)) |>
			
			# nest for output
			nest(data = -cell_type_harmonised) |> 
			
			# bind
			mutate(data = map(
				data,
				~ {
					
					se = do.call(cbind,  .x |> pull(data))
					
					# Filter very rare gene-transcripts
					all_zeros = assay(se, "counts") |> rowSums() |> equals(0) 
					se = se[!all_zeros,]
					lower_than_total_counts = assay(se, "counts") |> rowSums() < 15 
					se = se[!lower_than_total_counts,]
					
					# Make it sparse
					se@assays@data$counts = as(se@assays@data$counts, "sparseMatrix")
					
					se
				}
			))
		
		
	}
	
	map_quantile_scale_abundance = function(se_df){
		
		print("Start scale abundance")
		gc()
		
		se_df |> 
			
			# Quantile tranformation
			mutate(data = map(data,	quantile_normalise_abundance, method = "preprocesscore_normalize_quantiles_use_target")) |> 
			
			# convert to sparse again
			mutate(data = map(data, ~ {
				
				.x@assays@data$counts_scaled = as(.x@assays@data$counts_scaled, "sparseMatrix")
				
				attr(.x, "internals")$tt_columns$.abundance_scaled |> attr(".Environment") = NULL # new_environment()
				
				.x
			}))  
		
	}
	
	map_keep_abundant = function(se_df){
		
		print("Start keep abundant")
		gc()
		
		
		se_df |>
			mutate(data = map(
				data,
				~ {
					
					# For blood
					if(ncol(.x) > 1000)
						.x |>
						keep_abundant(
							.abundance = counts_scaled, 
							factor_of_interest = c(sex, ethnicity_simplified), 
							minimum_counts = 500,
							minimum_proportion = 0.9
						)
					
					else {
						# Select abundant genes within tissues and unite
						abundant_genes =
							.x |>
							nest(data = -cell_type_harmonised) |>
							
							# Filter if only one sex
							mutate(abundant_genes = map(
								data,
								~ {
									
									se = .x
									
									# # Avoid indeterminability
									# if(ncol(se) > 0) {
									#
									# 	ethnicity_to_keep = se |> pivot_sample() |> dplyr::count(ethnicity_simplified) |> filter(n >1) |> pull(ethnicity_simplified)
									# 	se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)
									#
									# 	sex_to_keep = se |> pivot_sample() |> dplyr::count(sex) |> filter(n >1) |> pull(sex)
									# 	se = se |> filter(sex %in% sex_to_keep)
									#
									# 	se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
									#
									# }
									#
									# factors =
									# 	c( "sex", "ethnicity_simplified") |>
									# 	enframe(value = "factor") |>
									# 	mutate(n = map_int(
									# 		factor, ~ se |> select(.x) |> distinct() |> nrow()
									# 	)) |>
									# 	filter(n>1) |>
									# 	pull(factor) |>
									# 	map(sym) |>
									# 	unlist()
									
									se |>
										keep_abundant(.abundance = counts_scaled, factor_of_interest = c(sex, ethnicity_simplified), minimum_counts = 50) |>
										rownames()
									
								},
								.progress = TRUE
							)) |>
							pull(abundant_genes) |>
							unlist() |>
							unique()
						
						.x |> filter(.feature %in% abundant_genes)
					}
					
				} 
			))
		# se_df |> 
		# 	mutate(data = map(
		# 		data,
		# 		~ {
		# 			
		# 			# Select abundant genes within tissues and unite
		# 			abundant_genes =
		# 				.x |>
		# 				nest(data = -cell_type_harmonised) |>
		# 				
		# 				# Filter if only one sex
		# 				mutate(abundant_genes = map(
		# 					data,
		# 					~ {
		# 						
		# 						se = .x
		# 						
		# 						# # Avoid indeterminability
		# 						# if(ncol(se) > 0) {
		# 						# 	
		# 						# 	ethnicity_to_keep = se |> pivot_sample() |> dplyr::count(ethnicity_simplified) |> filter(n >1) |> pull(ethnicity_simplified)
		# 						# 	se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)
		# 						# 	
		# 						# 	sex_to_keep = se |> pivot_sample() |> dplyr::count(sex) |> filter(n >1) |> pull(sex)
		# 						# 	se = se |> filter(sex %in% sex_to_keep)
		# 						# 	
		# 						# 	se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
		# 						# 	
		# 						# }
		# 						# 
		# 						# factors =
		# 						# 	c( "sex", "ethnicity_simplified") |>
		# 						# 	enframe(value = "factor") |>
		# 						# 	mutate(n = map_int(
		# 						# 		factor, ~ se |> select(.x) |> distinct() |> nrow()
		# 						# 	)) |>
		# 						# 	filter(n>1) |>
		# 						# 	pull(factor) |>
		# 						# 	map(sym) |> 
		# 						# 	unlist() 
		# 						
		# 						se |>
		# 						keep_abundant(.abundance = counts_scaled, factor_of_interest = c(sex, ethnicity_simplified), minimum_counts = 50) |>
		# 						rownames()
		# 						
		# 						},
		# 					.progress = TRUE
		# 				)) |>
		# 				pull(abundant_genes) |>
		# 				unlist() |>
		# 				unique()
		# 			
		# 			.x |> filter(.feature %in% abundant_genes)
		# 			
		# 		} 
		# 	))
	}
	
	se_add_dispersion = function(se_df){
		
		print("Start add dispersion")
		gc()
		
		se_df |>
			mutate(data = map(
				data,
				~ {
					
					# Because I have nested map
					se = .x
					
					# Avoid indeterminability
					if(ncol(se) > 0) {
						
						ethnicity_to_keep = se |> pivot_sample() |> dplyr::count(ethnicity_simplified) |> filter(n >1) |> pull(ethnicity_simplified)
						se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)
						
						sex_to_keep = se |> pivot_sample() |> dplyr::count(sex) |> filter(n >1) |> pull(sex)
						se = se |> filter(sex %in% sex_to_keep)
					}
					if(ncol(se) > 0) {
						se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
					}
					
					# If SE empty add dummy dispersion
					if(ncol(se) == 0) rowData(se)$dispersion = rep(NA, nrow(se))
					else if(ncol(se)<1000) {
						
						factors =
							c( "sex", "ethnicity_simplified") |>
							enframe(value = "factor") |>
							mutate(n = map_int(
								factor, ~ se |> select(.x) |> distinct() |> nrow()
							)) |>
							filter(n>1) |>
							pull(factor) |>
							str_c(collapse = " + ")
						
						my_design = glue("~ {factors}") |> as.formula() |> model.matrix( data = colData(se) |> droplevels()) 
						rowData(se)$dispersion =  se |> assay("counts_scaled") |> estimateDisp(design = my_design) %$% tagwise.dispersion  
					} 
					else {
						
						sampled_samples = sample(seq_len(ncol(se)), size = min(ncol(se), 2000))
						
						factors =
							c( "sex", "ethnicity_simplified") |>
							enframe(value = "factor") |>
							mutate(n = map_int(
								factor, ~ se[,sampled_samples, drop=FALSE]  |> select(.x) |> distinct() |> nrow()
							)) |>
							filter(n>1) |>
							pull(factor) |>
							str_c(collapse = " + ")
						
						my_design = model.matrix(~ sex + ethnicity_simplified, data = se[,sampled_samples, drop=FALSE] |> colData() |> droplevels()) 
						
						rowData(se)$dispersion =  
							assay(se[,sampled_samples, drop=FALSE] , "counts_scaled") |> 
							estimateTrendedDisp(design = my_design, subset=1000, rowsum.filter=10)
					}
					
					se
				}
			))
		
	}
	
	add_number_of_gene_chunks =	function(se_df){
		
		print("Start add number of chunks")
		gc()
		
		se_df |> 
			mutate(number_of_chunks = map_int(
				data, 
				~ .x |> 
					distinct(sample_se) |>
					nrow() |>
					
					# Avoid 0 because it will cause error
					# And because ifI have 0 samples I still have one chunk of genes
					max(1) |> 
					divide_by(50) |> 
					multiply_by(10) |> 
					ceiling() |> 
					
					# Otherwise it takes more than 20 minutes
					min(10000)
			))
	}
	
	# SEX TISSUE
	map_avoid_confounders_sex_tissue = function(se_df){
		
		se_df |> 
			mutate(data = map(
				data,
				~ {
					
					if(ncol(.x) <= 1) return(NULL)
					
					# Filter
					se = 
						.x |> 
						filter(sex != "unknown") |>
						
						
						# Eliminate complete confounders
						resolve_complete_confounders_of_non_interest(assay_simplified, disease, ethnicity_simplified) 
						
					rm(.x)
					gc()
					
					# Filter disease
					se =
						se |>
						filter(disease %in% (
							se |>
								distinct(disease, sex) |>
								dplyr::count(disease) |>
								filter(n>1) |>
								pull(disease)
						))
					
					# Return prematurely
					if(ncol(se) == 0) return(NULL)
					if(se |> distinct(sex, ethnicity_simplified) |> dplyr::count(ethnicity_simplified) |> pull(n) |> max() == 1) return(NULL)
					
					
					# Vell types with enough samples
					cell_type_to_keep =
						se |>
						distinct(sample_, cell_type_harmonised) |>
						dplyr::count(  cell_type_harmonised) |>
						filter(n > 3) |>
						pull(cell_type_harmonised)
					if(length(cell_type_to_keep)== 0) return(NULL)
					
					
					se = 
						se |>
						
						# Scale continuous variables
						mutate(age_days = scale(age_days) |> as.numeric()) |>
						
						# Filter cell types to keep
						filter(cell_type_harmonised %in% cell_type_to_keep) |>
						
						# otherwise I get error for some reason
						mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised", "cell_type_harmonised")), as.character)) |>
						mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) 
					
					# Drop random effect grouping with no enough data
					combinations_to_keep = 
						se |>
						distinct(sex, ethnicity_simplified, cell_type_harmonised) |>
						add_count(cell_type_harmonised) |>
						filter(n>1)
					
					se |>
						right_join(combinations_to_keep) 
					
				}))
		
	}
	
	map_create_formula_sex_tissue = function(se_df){
		
		se_df |> 
			mutate(formula = map(
				data,
				~ {
					
					if(is.null(.x)) return(~1)
					if(ncol(.x) <= 1) return(~1)
					if(.x |> distinct(sex, ethnicity_simplified) |> dplyr::count(ethnicity_simplified) |> pull(n) |> max() == 1) return(~1)
					
					se = .x
					
					# Build the formula
					factors =
						c("age_days", "sex", "ethnicity_simplified", "assay_simplified",  ".aggregated_cells", "disease") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(.x) |> distinct() |> nrow()
						)) |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					random_effects =
						c("age_days", "sex", "ethnicity_simplified") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
						))   |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					# The default
					my_formula = glue("~ {factors}")
					
					if(
						se |> distinct(cell_type_harmonised) |> nrow() > 1 &
						length(random_effects) > 0
					) 
						my_formula = glue("{my_formula} + (1 + {random_effects} | cell_type_harmonised)")
					
					if( 
						se |> distinct(sample_) |> nrow() > 1	&& 
						se |> distinct(cell_type_harmonised) |> nrow() > 1
					)
						my_formula = glue("{my_formula} + (1 | sample_)")
					
					
					# Add the interaction
					if(se |> 
						 nest_detect_complete_confounder(age_days, sex) |> 
						 filter(n1 + n2 <= 2) |> 
						 nrow() == 0
					) 
						my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")
					
					
					return(my_formula)
					
				}))
		
	}
	
	# SEX CELL TYPE
	
	map_avoid_confounders_sex_cell_type = function(se_df){
		
		se_df |> 
			mutate(data = map(
				data,
				~ {
					
					if(ncol(.x) <= 1) return(NULL)
					
					# Filter
					se = 
						.x |> 
						filter(sex != "unknown") |>
						
						
						# Eliminate complete confounders
						resolve_complete_confounders_of_non_interest(assay_simplified, disease, ethnicity_simplified) 
					
					
					rm(.x)
					gc()
					
					# Filter disease
					se =
						se |>
						filter(disease %in% (
							se |>
								distinct(disease, sex) |>
								dplyr::count(disease) |>
								filter(n>1) |>
								pull(disease)
						))
					
					# Return prematurely
					if(ncol(se) == 0) return(NULL)
					if(se |> distinct(sex, ethnicity_simplified) |> dplyr::count(ethnicity_simplified) |> pull(n) |> max() == 1) return(NULL)
					
					
					# Vell types with enough samples
					tissues_to_keep =
						se |>
						distinct(sample_, tissue_harmonised) |>
						dplyr::count(  tissue_harmonised) |>
						filter(n > 3) |>
						pull(tissue_harmonised)
					if(length(tissues_to_keep)== 0) return(NULL)
					
					se =
						se |>
						
						# Scale contninuous variables
						mutate(age_days = scale(age_days) |> as.numeric()) |>
						
						# Filter cell types to keep
						filter(tissue_harmonised %in% tissues_to_keep) |>
						
						# otherwise I get error for some reason
						mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised")), as.character)) |>
						mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) 
					
					# Drop random effect grouping with no enough data
					combinations_to_keep = 
						se |>
						distinct(sex, ethnicity_simplified, tissue_harmonised) |>
						add_count(tissue_harmonised) |>
						filter(n>1)
					
					
					se |>
						right_join(combinations_to_keep) 
					
				}))
	}
	
	map_create_formula_sex_cell_type = function(se_df){
		
		se_df |> 
			mutate(formula = map(
				data,
				~ {
					
					if(is.null(.x)) return(~1)
					if(ncol(.x) <= 1) return(~1)
					if(.x |> distinct(sex, ethnicity_simplified) |> dplyr::count(ethnicity_simplified) |> pull(n) |> max() == 1) return(~1)
					
					se = .x
					
					# Build the formula
					factors = 
						c("age_days", "sex", "ethnicity_simplified", "assay_simplified",  ".aggregated_cells", "disease") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(.x) |> distinct() |> nrow()
						)) |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					random_effects =
						c("age_days", "sex", "ethnicity_simplified") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
						))   |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					
					# The default
					my_formula = glue("~ {factors}")
					
					if( 
						se |> distinct(tissue_harmonised) |> nrow() > 1 &
						length(random_effects) > 0
					) 
						my_formula = glue("{my_formula} + (1 + {random_effects} | tissue_harmonised)")
					
					if( 	se |> distinct(file_id) |> nrow() > 1	)
						my_formula = glue("{my_formula} + (1 | file_id)")
					
					# Add the interaction
					if(se |> 
						 nest_detect_complete_confounder(age_days, sex) |> 
						 filter(n1 + n2 <= 2) |> 
						 nrow() == 0
					) 
						my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")
					
					my_formula
					
				}))
	}
	
	# 
	
	samples_NOT_complete_confounders_for_age_ethnicity = function(se){
		
			
		clean = 
			se |> 
			#distinct(ethnicity_simplified, assay_simplified, .sample) |>
			
			nest(se_data = -c(ethnicity_simplified, age_days)) |>
			
			# How many ethnicity per assay
			nest(data = -ethnicity_simplified) |> 
			mutate(n1 = map_int(data, ~ .x |> distinct(age_days) |> nrow())) |> 
			unnest(data) |> 
			
			# How many age days per ethnicity
			nest(data = - age_days) |> 
			mutate(n2 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
			unnest(data) |>
			
			filter(n1+n2>2) |>
			select(-n1, -n2) 
		
		if(
			nrow(clean) == 0 || clean |>
			distinct(ethnicity_simplified) |>
			nrow() == 1
		){
			
			# If I just have one observation
			if(	se |> 
					#distinct(ethnicity_simplified, assay_simplified, .sample) |>
					mutate(age_days = age_days > 1) |>
					distinct(ethnicity_simplified, age_days) |> 
					nrow() < 2
			) return(se)
			
			clean = 
				se |> 
				#distinct(ethnicity_simplified, assay_simplified, .sample) |>
				mutate(age_days = age_days > 1) |>
				nest(se_data = -c(ethnicity_simplified, age_days)) |>
				
				# How many ethnicity per assay
				nest(data = -ethnicity_simplified) |> 
				mutate(n1 = map_int(data, ~ .x |> distinct(age_days) |> nrow())) |> 
				unnest(data) |> 
				
				# How many assay per ethnicity
				nest(data = - age_days) |> 
				mutate(n2 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
				unnest(data) |>
				
				filter(n1+n2>2) |>
				select(-n1, -n2) |> 
				mutate(age_days = age_days |> as.integer())
			
			
		}
		
		clean |> unnest_summarized_experiment(se_data) 
			
	}
	
	map_avoid_confounders_ethnicity_tissue = function(se_df){
		
		se_df |> 
			mutate(data = map(
				data,
				~ {
					
					if(ncol(.x) <= 1) return(NULL)
					
					# Filter
					se = 
						.x |> 
						
						# Eliminate complete confounders
						resolve_complete_confounders_of_non_interest(assay_simplified, disease, sex) |> 
						samples_NOT_complete_confounders_for_age_ethnicity()
					
					
					rm(.x)
					gc()
					
					# ethnicity to keep
					ethnicity_to_keep = se |> pivot_sample() |> dplyr::count(ethnicity_simplified) |> filter(n>1) |> pull(ethnicity_simplified)
					se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)
					
					# Filter disease
					se =
						se |>
						filter(disease %in% (
							se |>
								distinct(disease, ethnicity_simplified) |>
								dplyr::count(disease) |>
								filter(n>1) |>
								pull(disease)
						))
					
					# Return prematurely
					if(ncol(se) == 0) return(NULL)
					if(se |> distinct(sex, ethnicity_simplified) |> dplyr::count(sex) |> pull(n) |> max() == 1) return(NULL)
					
					
					# Vell types with enough samples
					cell_type_to_keep =
						se |>
						distinct(sample_, cell_type_harmonised) |>
						dplyr::count(  cell_type_harmonised) |>
						filter(n > 3) |>
						pull(cell_type_harmonised)
					if(length(cell_type_to_keep)== 0) return(NULL)
					
					se = 
						se |>
						
						# Scale continuous variables
						mutate(age_days = scale(age_days) |> as.numeric()) |>
						
						# Filter cell types to keep
						filter(cell_type_harmonised %in% cell_type_to_keep) |>
						
						# otherwise I get error for some reason
						mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised", "cell_type_harmonised")), as.character)) |>
						mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) 
					
					# Drop random effect grouping with no enough data
					combinations_to_keep = 
						se |>
						distinct(sex, ethnicity_simplified, cell_type_harmonised) |>
						add_count(cell_type_harmonised) |>
						filter(n>1)
					
					
					se |>
						right_join(combinations_to_keep) 
					
					
				}))
		
	}
	
	map_create_formula_ethnicity_tissue = function(se_df){
		
		se_df |> 
			mutate(formula = map(
				data,
				~ {
					
					if(is.null(.x)) return(~1)
					if(ncol(.x) <= 1) return(~1)
					if(.x |> distinct(sex, ethnicity_simplified) |> dplyr::count(sex) |> pull(n) |> max() == 1) return(~1)
					
					se = .x
					
					# Build the formula
					factors =
						c("age_days", "sex", "ethnicity_simplified", "assay_simplified",  ".aggregated_cells", "disease") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
						)) |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					random_effects =
						c("age_days", "sex", "ethnicity_simplified") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
						))   |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					# The default
					my_formula = glue("~ {factors}")
					
					if(
						se |> distinct(cell_type_harmonised) |> nrow() > 1 &
						length(random_effects) > 0
					) 
						my_formula = glue("{my_formula} + (1 + {random_effects} | cell_type_harmonised)")
					
					
					if( 	
						se |> distinct(sample_) |> nrow() > 1	&& 
						se |> distinct(cell_type_harmonised) |> nrow() > 1
					)
						my_formula = glue("{my_formula} + (1 | sample_)")
					
					
					# Add the interaction
					if(se |> 
						 nest_detect_complete_confounder(age_days, sex) |> 
						 filter(n1 + n2 <= 2) |> 
						 nrow() == 0
					) 
						my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")
					
					my_formula
					
				}))
		
	}
	
	map_avoid_confounders_ethnicity_cell_type = function(se_df){
		
		se_df |> 
			mutate(data = map(
				data,
				~ {
					
					if(ncol(.x) <= 1) return(NULL)
					
					# Filter
					se = 
						.x |> 
						
						# Eliminate complete confounders
						resolve_complete_confounders_of_non_interest(assay_simplified, disease, sex) 
					
					
					rm(.x)
					gc()
					
					# Filter disease
					se =
						se |>
						filter(disease %in% (
							se |>
								distinct(disease, ethnicity_simplified) |>
								dplyr::count(disease) |>
								filter(n>1) |>
								pull(disease)
						))
					
					
					
					# Return prematurely
					if(ncol(se) == 0) return(NULL)
					if(se |> distinct(sex, ethnicity_simplified) |> dplyr::count(sex) |> pull(n) |> max() == 1) return(NULL)
					
					# Vell types with enough samples
					tissues_to_keep =
						se |>
						distinct(sample_, tissue_harmonised) |>
						dplyr::count(  tissue_harmonised) |>
						filter(n > 3) |>
						pull(tissue_harmonised)
					
					if(length(tissues_to_keep)== 0) return(NULL)
					
					se =
						se |>
						
						# Scale contninuous variables
						mutate(age_days = scale(age_days) |> as.numeric()) |>
						
						# Filter cell types to keep
						filter(tissue_harmonised %in% tissues_to_keep) |>
						
						# otherwise I get error for some reason
						mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised")), as.character)) |>
						mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) 
					
					# Drop random effect grouping with no enough data
					combinations_to_keep = 
						se |>
						distinct(sex, ethnicity_simplified, tissue_harmonised) |>
						add_count(tissue_harmonised) |>
						filter(n>1)
					
					
					se |>
						right_join(combinations_to_keep) 
					
				}))
	}
	
	map_create_formula_ethnicity_cell_type = function(se_df){
		
		se_df |> 
			mutate(formula = map(
				data,
				~ {
					
					if(is.null(.x)) return(~1)
					if(ncol(.x) <= 1) return(~1)
					if(.x |> distinct(sex, ethnicity_simplified) |> dplyr::count(sex) |> pull(n) |> max() == 1) return(~1)
					
					se = .x
					
					# Build the formula
					factors = 
						c("age_days", "sex", "ethnicity_simplified", "assay_simplified", ".aggregated_cells", "disease") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(.x) |> distinct() |> nrow()
						)) |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					random_effects =
						c("age_days", "sex", "ethnicity_simplified") |>
						enframe(value = "factor") |>
						mutate(n = map_int(
							factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
						))   |>
						filter(n>1) |>
						pull(factor) |>
						str_c(collapse = " + ")
					
					# The default
					my_formula = glue("~ {factors}")
					
					if( 
						se |> distinct(tissue_harmonised) |> nrow() > 1 &
						length(random_effects) > 0
					) 
						my_formula = glue("{my_formula} + (1 + {random_effects} | tissue_harmonised)")
					
					
					if( 	se |> distinct(file_id) |> nrow() > 1	)
						my_formula = glue("{my_formula} + (1 | file_id)")
					
					
					# Add the interaction
					if(se |> 
						 nest_detect_complete_confounder(age_days, sex) |> 
						 filter(n1 + n2 <= 2) |> 
						 nrow() == 0
					) 
						my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")
					
					my_formula		
					
				}))
	}
	
	#-----------------------#
	# Pipeline
	#-----------------------#
	list(
		
		# Do metadata
		tarchetypes::tar_group_by(
			tissue_cell_type_metadata,
			glue("{root_directory}/sccomp_on_HCA_0.2.3.4/input_relative.rds") |> 
				split_metadata(),
			cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_chunk, 
			deployment = "main",
			resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
		),
		
		# Get SCE
		tar_target(
			pseudobulk_df,
			tissue_cell_type_metadata |> get_sce() |> get_pseudobulk(),
			pattern = map(tissue_cell_type_metadata),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("local_30"))
		),
		
		
		#-----------------------#
		# TISSUE ANALYSES
		#-----------------------#
		
		# Group samples
		tarchetypes::tar_group_by(
			pseudobulk_df_tissue, 
			pseudobulk_df, 
			tissue_harmonised, is_immune,
			resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
		),
		
		# Add dispersion
		tar_target(
			pseudobulk_df_scaled_abundant_tissue,
			pseudobulk_df_tissue |> 
				aggregate() |> 
				map_quantile_scale_abundance() |> 
				map_keep_abundant() ,
			pattern = map(pseudobulk_df_tissue),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
		),
		
		# Make data SEX TISSUE
		tar_target(
			pseudobulk_df_scaled_abundant_curated_formula_tissue_sex,
			pseudobulk_df_scaled_abundant_tissue |> 
				map_avoid_confounders_sex_tissue() |> 
				map_create_formula_sex_tissue(),
			pattern = map(pseudobulk_df_scaled_abundant_tissue),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("slurm_2_20"))
		),
		
		# Make data SEX TISSUE
		tar_target(
			pseudobulk_df_scaled_abundant_curated_formula_tissue_ethnicity,
			pseudobulk_df_scaled_abundant_tissue |> 
				map_avoid_confounders_ethnicity_tissue() |> 
				map_create_formula_ethnicity_tissue(),
			pattern = map(pseudobulk_df_scaled_abundant_tissue),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("slurm_2_20"))
		),
		
		#-----------------------#
		# IMMUNE CELL TYPE ANALYSES
		#-----------------------#
		
		# Group samples
		tarchetypes::tar_group_by(
			pseudobulk_df_cell_type, 
			pseudobulk_df |> 
				filter(is_immune == "TRUE"), 
			cell_type_harmonised, is_immune,
			resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
		),
		
		# Add dispersion
		tar_target(
			pseudobulk_df_scaled_abundant_cell_type,
			pseudobulk_df_cell_type |> 
				aggregate_cell_type() |> 
				map_quantile_scale_abundance() |> 
				map_keep_abundant(),
			pattern = map(pseudobulk_df_cell_type),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("slurm_1_40"))
		),
		
		# Make data SEX TISSUE
		tar_target(
			pseudobulk_df_scaled_abundant_curated_formula_cell_type_sex,
			pseudobulk_df_scaled_abundant_cell_type |> 
				map_avoid_confounders_sex_cell_type() |> 
				map_create_formula_sex_cell_type(),
			pattern = map(pseudobulk_df_scaled_abundant_cell_type),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("slurm_2_20"))
		),
		
		# Make data SEX TISSUE
		tar_target(
			pseudobulk_df_scaled_abundant_curated_formula_cell_type_ethnicity,
			pseudobulk_df_scaled_abundant_cell_type |> 
				map_avoid_confounders_ethnicity_cell_type() |> 
				map_create_formula_ethnicity_cell_type(),
			pattern = map(pseudobulk_df_scaled_abundant_cell_type),
			iteration = "group",
			resources = tar_resources(crew = tar_resources_crew("slurm_1_120"))
		)
		
	)
	
	
}, ask = FALSE, script = glue("{result_directory}/_targets__pseudobulk_non_immune_split_Dec2023.R"))



#job::job({
# result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"

tar_make(
	 callr_function = NULL,
	script = glue("{result_directory}/_targets__pseudobulk_non_immune_split_Dec2023.R"),
	store = glue("{result_directory}/_targets__pseudobulk_non_immune_split_Dec2023")
)



#  "/stornext/General/scratch/GP_Transfer/michael_targets/"
