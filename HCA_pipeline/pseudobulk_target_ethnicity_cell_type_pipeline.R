
library(targets)
library(glue)
result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4"



tar_script(	{
	
	library(targets, tarchetypes, crew, "crew.cluster")
	library(tarchetypes)
	library(crew)
	library(crew.cluster)
	library(tidyverse)
	library(stringr)
	library(glue)
	library(tidySummarizedExperiment)
	
	result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4"
	
	#-----------------------#
	# Future SLURM
	#-----------------------#
	
	library(future)
	library("future.batchtools")
	slurm <- 
		`batchtools_slurm` |>
		future::tweak( template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
									 resources=list(
									 	ncpus = 20,
									 	memory = 6000,
									 	walltime = 172800
									 )
		)
	plan(slurm)
	
	#-----------------------#
	# Future SLURM
	#-----------------------#
	
	# small_slurm =
	# 	crew_controller_slurm(
	# 		name = "small_slurm",
	# 		slurm_memory_gigabytes_per_cpu = 8,
	# 		slurm_cpus_per_task = 1,
	# 		workers = 200,
	# 		verbose = T,
	# 		#script_lines = "module load R/4.2.1", 
	# 		host = "spartan.hpc.unimelb.edu.au"
	# 	)
	# 
	# big_slurm =
	# 	crew_controller_slurm(
	# 		name = "big_slurm",
	# 		slurm_memory_gigabytes_per_cpu = 4,
	# 		slurm_cpus_per_task = 11,
	# 		workers = 200,
	# 		verbose = T,
	# 		#script_lines = "module load R/4.2.1", 
	# 		host = "spartan.hpc.unimelb.edu.au"
	# 	)
	
	#-----------------------#
	# Packages
	#-----------------------#
	tar_option_set( 
		packages = c("tidybulk", "stringr", "tibble", "tidySummarizedExperiment", "purrr", "glue", "dplyr"),
		storage = "worker", 
		retrieval = "worker", 
		error = "continue"		
		# ,
		# controller = crew_controller_group(small_slurm, big_slurm),
		# resources = tar_resources(crew = tar_resources_crew("small_slurm"))
	)
	
	#-----------------------#
	# Functions
	#-----------------------#
	
	build_pseudobulk = function(metadata, tissue_harmonised, cell_type_harmonised, is_immune){
		
		# Get seurat from CuratedAtlasQuery
		seu = 				
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
		
		# Calculate median number of genes genes per samples
		stats = 
			seu |>
			select(sample_, nFeature_originalexp) |>
			with_groups(sample_, ~ .x |>
										summarise(
											mean_nFeature = mean(nFeature_originalexp), 
											median_nFeature = median(nFeature_originalexp)
										)) 
		
		
		seu |> 
			
			# Calculate summary stats
			left_join(stats, by = "sample_") |>
			mutate(sample_ = glue("{sample_}___{cell_type_harmonised}")) |>
			tidyseurat::aggregate_cells( .sample = sample_,
																	 .metadata_columns_to_keep = c('cell_type_harmonised', 'sample_', 'file_id', 'assay', 'age_days', 'development_stage', 'sex', 'ethnicity', 'tissue_harmonised', 'tissue', 'disease', 'lineage_1', 'is_immune', 'sample_count', 'age_days_original', 'ethnicity_simplified', 'assay_simplified', 'mean_nFeature', 'median_nFeature')
			) |>
			rename(counts = originalexp) |>
			tidybulk::as_SummarizedExperiment(.sample, .feature, counts) 
		
	}
	
	samples_NOT_complete_confounders_for_sex_assay = function(se){
		
		
		
		se = 
			se |> 
			# distinct(assay_simplified, sex, .sample) |>
			# 
			nest(se_data = -c(assay_simplified, sex)) |>
			
			# How many ethnicity per assay
			nest(data = -assay_simplified) |> 
			mutate(n1 = map_int(data, ~ .x |> distinct(sex) |> nrow())) |> 
			unnest(data) |> 
			
			# How many assay per ethnicity
			nest(data = - sex) |> 
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
	
	samples_NOT_complete_confounders_for_sex_disease = function(se){
		
		
		
		se = 
			se |> 
			#distinct(disease, sex, .sample) |>
			
			nest(se_data = -c(disease, sex)) |>
			
			# How many ethnicity per assay
			nest(data = -disease) |> 
			mutate(n1 = map_int(data, ~ .x |> distinct(sex) |> nrow())) |> 
			unnest(data) |> 
			
			# How many assay per ethnicity
			nest(data = - sex) |> 
			mutate(n2 = map_int(data, ~ .x |> distinct(disease) |> nrow())) |> 
			unnest(data) 
		
		
		# Replace ethnicity
		dummy_ethnicity = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(sex)
		
		se |>
			mutate(sex = if_else(n1 + n2 < 3, dummy_ethnicity, sex)) 	|>
			
			# # Filter
			# filter(!(n1==1 & n2==1)) |>
			select(-n1, -n2) |>
			
			unnest_summarized_experiment(se_data)
		# |>
		# 	pull(.sample) |>
		# 	unique()
	}
	
	merge_and_filter = function(se_files){
		se = 
		do.call(
			cbind, 
			glue("{result_directory}/{se_files}") |> 
				map(
					~ .x |>
						readRDS() |>
						select(-any_of("orig.ident")) |>
						select(
							.aggregated_cells, .feature,median_nFeature, .sample, age_days, assay_simplified,cell_type_harmonised,
							counts, disease, ethnicity_simplified , file_id, sex, tissue_harmonised
						)
				) 
		)	
		
		# Filter
		se = 
			se |> 
			filter(sex != "unknown") |>
			filter(.aggregated_cells > 10) 
		
		# Samble to subset
		# Filter at least 5000 genes recorded
		se = 
			se |>
			filter(.sample %in% (
				se |>
					assay("counts") |> 
					apply(2, function(x) (x>0) |> which() |> length()) |>
					enframe() |>
					mutate(value = as.character(value), name = as.character(name)) |>
					filter(value > 5000) |>
					pull(name)
			)) 
		
		# Eliminate complete confounders
		se = 
			se |>
			samples_NOT_complete_confounders_for_sex_assay() |>
			samples_NOT_complete_confounders_for_sex_disease()
		
		# Filter ethnicity_simplified
		se = 
			se |>
			filter(disease %in% (
				se |>
					distinct(disease, ethnicity_simplified) |>
					count(disease) |>
					filter(n>1) |>
					pull(disease)
			))
		
		# Filter tissue that has two sexes
		se = 
			se |>
			filter(tissue_harmonised %in% (
				se |>
					distinct(tissue_harmonised, ethnicity_simplified) |>
					count(tissue_harmonised) |>
					filter(n>1) |>
					pull(tissue_harmonised)
			))
		
		se
		
	}
	
	analyse = function(se, max_rows_for_matrix_multiplication = NULL, cores = 1){
		
		if(ncol(se) == 0) return(se)
		
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
			se |> distinct(tissue_harmonised) |> nrow() > 1 &
			length(random_effects) > 0
		) {
			my_formula = glue("{my_formula} + (1 + {random_effects} | tissue_harmonised)")
			method = "glmmseq_lme4"
		}
		
		if( 	se |> distinct(file_id) |> nrow() > 1	){
			my_formula = glue("{my_formula} + (1 | file_id)")
			method = "glmmseq_lme4"
		}
		
		# Normalise
		se = se |> quantile_normalise_abundance() 
		
		# Select abundant genes within tissues and unite
		abundant_genes = 
			se |>
			nest(data = -tissue_harmonised) |>
			mutate(abundant_genes = map(
				data,
				~ .x |>
					keep_abundant(.abundance = counts_scaled, factor_of_interest = c(ethnicity_simplified)) |>
					rownames(), 
				.progress = TRUE
			)) |>
			pull(abundant_genes) |>
			unlist() |>
			unique()
		
		se =
			se |>
			
			# Scale contninuous variables
			mutate(age_days = scale(age_days) |> as.numeric()) |>
			
			# Filter abundant genes
			filter(.feature %in% abundant_genes) |>
			
			# otherwise I get error for some reason
			mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised")), as.character)) |>
			mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) |>
			# Test
			test_differential_abundance(
				as.formula(my_formula),
				.abundance = counts_scaled,
				method = method,
				cores = cores, 
				max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication
			)
	}
	

	
	
	
	#-----------------------#
	# Workflow
	#-----------------------#
	tar_map(
		
		# Input
		values = 
			result_directory |>
			dir( pattern = ".rds") |>
			enframe(value = "file") |>
			tidyr::extract(col = file, into = c("tissue", "cell_type"), regex = "/?([a-zA-Z]+)____(.+)____TRUE\\.rds", remove = FALSE) |> 
			filter(!is.na(cell_type)) |>
			select(cell_type, file) |>
			nest(cell_type_files = -cell_type) |>
			mutate(cell_type_files = map(cell_type_files, ~.x$file)),
		
		# Names jobs
		names = "cell_type",
		
		# Filter
		tar_target(se_filtered, merge_and_filter(cell_type_files)),
		
		# Dispersion
		tar_target(dispersion, se_to_dispersion(se_filtered)),
		
		# Estimate
		tar_target(
			data, analyse(se_filtered, max_rows_for_matrix_multiplication = 10000, cores = 18)		
			# ,
			# resources = tar_resources(crew = tar_resources_crew("big_slurm"))
			
		)
	)
}, ask = FALSE, script = glue("{result_directory}/_targets__ethnicity.R") )


# job::job({

tar_make_future(
	script = glue("{result_directory}/_targets__ethnicity.R"),
	store = glue("{result_directory}/_targets__ethnicity"), 
	workers = 200, 
	garbage_collection = TRUE
)

# tar_make(
# 	script = glue("{result_directory}/_targets__ethnicity.R"), 
# 	store = glue("{result_directory}/_targets__ethnicity")
# )
# 	
# })

