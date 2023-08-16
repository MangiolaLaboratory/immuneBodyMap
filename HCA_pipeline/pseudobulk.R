
library(targets)
result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4"


tar_script({

	library(targets)
	library(tarchetypes)
	library(targets)
	library(stringr)
	library(crew)
	library(glue)
	library(tidyverse)
	library(crew.cluster)

	result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4"


	#-----------------------#
	# Future SLURM
	#-----------------------#

	# library(future)
	# library("future.batchtools")
	# slurm <- future::tweak(batchtools_slurm,
	# 											 template = glue("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/third_party_sofware/slurm_batchtools.tmpl"),
	# 											 resources=list(
	# 											 	ncpus = 20,
	# 											 	memory = 5000,
	# 											 	walltime = 6000
	# 											 )
	# )
	# plan(slurm)


	small_slurm =
		crew_controller_slurm(
			name = "small_slurm",
			slurm_memory_gigabytes_per_cpu = 30,
			slurm_cpus_per_task = 1,
			workers = 200,
			verbose = TRUE,
			seconds_wall = 6000,
			#script_directory = glue("{result_directory}/_targets.R"),
			script_lines = "module load R/4.2.1",
		)

	big_slurm =
		crew_controller_slurm(
			name = "big_slurm",
			slurm_memory_gigabytes_per_cpu = 2,
			slurm_cpus_per_task = 20,
			workers = 200,
			verbose = TRUE,
			seconds_wall = 6000,
			#script_directory = glue("{result_directory}/_targets.R"),
			script_lines = "module load R/4.2.1",
		)

	#-----------------------#
	# Packages
	#-----------------------#
	tar_option_set(
		packages = c("tidybulk", "stringr", "tibble", "tidySummarizedExperiment", "purrr", "glue"),
		storage = "worker",
		retrieval = "worker",
		error = "continue",
		controller = crew_controller_group(small_slurm, big_slurm),
		resources = tar_resources(crew = tar_resources_crew("small_slurm"))
	)

	#-----------------------#
	# Functions
	#-----------------------#

	filter_samples_complete_confounders_for_ethnicity_assay = function(se){
		se |>
			nest(se_data = -c(assay_simplified, ethnicity_simplified)) |>

			# How many ethnicity per assay
			nest(data = -assay_simplified) |>
			mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |>
			unnest(data) |>

			# How many assay per ethnicity
			nest(data = - ethnicity_simplified) |>
			mutate(n2 = map_int(data, ~ .x |> distinct(assay_simplified) |> nrow())) |>
			unnest(data) |>

			# Filter
			filter(!(n1==1 & n2==1)) |>
			select(-n1, -n2) |>

			unnest_summarized_experiment(se_data)
	}

	merge_and_filter = function(se_files){
		se =
			glue("{result_directory}/{se_files}") |>
			map(
				~ .x |>
					readRDS() |>
					select(-any_of("orig.ident")) |>
					select(
						.aggregated_cells, .feature,median_nFeature, .sample, age_days, assay_simplified,cell_type_harmonised,
						counts, disease, ethnicity_simplified , file_id, sex, tissue_harmonised
					)
			) |>
			reduce(cbind)


		# Samble to subset
		# Filter at least 5000 genes recorded
		sample_to_subset =
			se |>
			assay("counts") |>
			apply(2, function(x) (x>0) |> which() |> length()) |>
			enframe() |>
			mutate(value = as.character(value), name = as.character(name)) |>
			filter(value > 5000) |>
			pull(name)

		# Filter
		se |>
			filter(sex != "unknown") |>
			filter(.sample %in% sample_to_subset) |>
			filter(.aggregated_cells > 10) |>
			filter_samples_complete_confounders_for_ethnicity_assay()

	}

	normalise = function(se){

		if(ncol(se) == 0) return(se)

		se |>
			mutate(age_days = scale(age_days) |> as.numeric()) |>
			mutate(sex_etnicity = glue("{sex}__{ethnicity_simplified}")) |>
			identify_abundant(factor_of_interest = "sex_etnicity") |>
			quantile_normalise_abundance()

	}

	analyse = function(se){

		if(ncol(se) == 0) return(se)

		factors =
			c("sex", "ethnicity_simplified", "assay_simplified", "age_days", ".aggregated_cells") |>
			enframe(value = "factor") |>
			mutate(n = map_int(
				factor, ~ se |> select(.x) |> distinct() |> nrow()
			)) |>
			filter(n>1) |>
			pull(factor) |>
			str_c(collapse = " + ")

		random_effects =
			c("sex") |> #, "ethnicity_simplified") |>
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

		se[sample(1:2000, size = 2000),] |>

			# otherwise I get error for some reason
			mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised")), as.character)) |>

			test_differential_abundance(
					as.formula(my_formula),
					.abundance = counts_scaled,
					method = method,
					cores = 20
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

		# Normalise
		tar_target(	se_normalised, normalise(se_filtered)),

		# Estimate
		tar_target(
			data, analyse(se_normalised),
			resources = tar_resources(crew = tar_resources_crew("big_slurm"))

		)
	)
}, ask = FALSE, script = glue("{result_directory}/_targets.R") )

# tar_make_future(
# 	script = glue("{result_directory}/_targets.R"),
# 	store = glue("{result_directory}/_targets"),
# 	workers = 200
# )

tar_make(
	script = glue("{result_directory}/_targets.R"),
	store = glue("{result_directory}/_targets")
)

se = tar_read(se_normalised_plasma, store = glue("{result_directory}/_targets")) |> analyse()


library(tidyverse)
library(tidybulk)
library(tidySummarizedExperiment)

# glue("{result_directory}/_targets__sex/objects") |>
# dir( pattern = "data_") |>
# 	map(~ .x |> tar_read(store = glue("{result_directory}/_targets__sex")))
#
#
# se = tar_read(se_filtered_granulocyte, store = glue("{result_directory}/_targets__sex")) |> analyse(max_rows_for_matrix_multiplication = 10000, cores = 18)
#
#
# se = tar_read(data_granulocyte, store = glue("{result_directory}/_targets__sex"))

de_sex =
	tar_meta(store = glue("{result_directory}/_targets__sex"), starts_with("data_")) |>
	filter(!is.na(data)) |>
	mutate(se = map(
		name,
		~ tar_read_raw(.x, store = glue("{result_directory}/_targets__sex")) |>
			pivot_transcript(),
		.progress=T
	))



# se =
# 	tar_read(se_filtered_granulocyte, store = glue("{result_directory}/_targets__sex")) |>
# 	samples_NOT_complete_confounders_for_ethnicity_disease() |>
# 	analyse(max_rows_for_matrix_multiplication = 1000, cores = 18)

se |>
	pivot_transcript() |>
	filter(P_sex_adjusted < 0.05) |>
	select(.feature, contains("sexmale")) |>

	pivot_longer(
		-c(.feature , sexmale),
		names_pattern = "([a-zA-Z\\.]+)_tissue_harmonised.+(lower|mode|upper)",
		names_to = c("tissue", "stat")
	) |>
	pivot_wider(names_from = stat, values_from = value) |>
	ggplot(aes(mode, tissue)) +
	geom_point() +
	facet_wrap(~.feature, scales = "free_x")


se |>
	pivot_transcript() |>
	select(.feature, contains("sex"), P_sex_adjusted) |>
	mutate(feature_to_print = if_else(P_sex_adjusted < 0.05, .feature, "")) |>
	arrange(P_sex_adjusted) |>
	ggplot(aes(sexmale, P_sex, label = feature_to_print)) +
	geom_point() +
	ggrepel::geom_text_repel() +
	scale_y_continuous(trans = log10_reverse_trans())

se |>
	filter(.feature== "NCOA7") |>
	ggplot(aes(sex, counts_scaled + 1)) +
	geom_boxplot(fill ="white", outlier.shape = NA) +
	geom_jitter(aes(color = tissue_harmonised), size = 0.2, height =0) +
	facet_wrap(~ tissue_harmonised, scales = "free_y") +
	scale_y_log10() +
	theme_bw()
