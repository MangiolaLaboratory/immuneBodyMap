
result_directory = "~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.1"

# AGE
differential_composition_age = readRDS(glue("{result_directory}/age_absolute_FALSE.rds"))

# Significance global statistics
count_significance_age_immune_load =
	differential_composition_age |>
	test_contrasts(test_composition_above_logit_fold_change = 0.017) |>
	filter(factor=="age_days") |>
	filter(is_immune=="TRUE") |>
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)


count_significance_age_immune_load_tissue =
	differential_composition_age |>
	
	# Find stats of random effect with groups
	test_contrasts(
		contrasts =
			differential_composition_age |>
			filter(parameter |> str_detect("___age_days")) |>
			distinct(parameter) |>
			mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
			deframe( ),
		test_composition_above_logit_fold_change = 0.1
	) |>
	filter(is_immune == "TRUE") |>
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

rm(differential_composition_age)
gc()

differential_composition_age_relative = readRDS(glue("{result_directory}/age_relative_FALSE.rds"))

count_significance_age_cell_type = 
	differential_composition_age_relative |>
	test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
	filter(parameter == "age_days") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

count_significance_age_cell_type_tissue = 
	differential_composition_age_relative |>
	test_contrasts(
		contrasts =
			differential_composition_age_relative |>
			filter(parameter |> str_detect("___age_days")) |>
			distinct(parameter) |>
			mutate(contrast = glue("age_days + {parameter}") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
			deframe( ),
		test_composition_above_logit_fold_change = 0.1
	)  |> 
	filter(parameter == "age_days") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

rm(differential_composition_age_relative)
gc()


# ETHNICITY
differential_composition_ethnicity_absolute = glue("{result_directory}/ethnicity_absolute_FALSE.rds") |> readRDS()


count_significance_ethnicity_immune_load = 
	differential_composition_ethnicity_absolute |> 
	test_contrasts(
		c(
			African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
			Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
			European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
			Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
		),
		test_composition_above_logit_fold_change = 0.1
	) |>
	
	filter(parameter %in% c("African",  "Hispanic", "European", "Chinese")) |>
	filter(is_immune == "TRUE") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)


count_significance_ethnicity_immune_load_tissue = 
	differential_composition_ethnicity_absolute |>
	
	# Find stats of random effect with groups
	test_contrasts(
		contrasts =
			differential_composition_ethnicity_absolute |>
			filter(parameter |> str_detect("___ethnicity")) |>
			distinct(parameter) |>
			tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
			mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
			mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
			select(tissue_harmonised, contrast) |> 
			separate(tissue_harmonised, c("tissue_harmonised", "ethnicity"), sep="_") |> 
			add_count(tissue_harmonised) |> mutate(contrast = glue("({contrast})")) |>  
			pivot_wider(names_from = ethnicity, values_from = contrast) |> 
			
			# Filter where there are at least 3 ethnicities
			mutate(n_minus_1 = n-1) |> 
			filter(n_minus_1>1) |>
			
			# Build contrasts
			rowwise() |> 
			mutate(
				other_Chinese = c(European, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
				other_European = c(Chinese, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
				other_African = c(European, Chinese, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
				`other_Hispanic or Latin American` = c(European, African, Chinese) |> str_subset(".+") |> str_c(collapse=" + "),
			) |> 
			mutate(
				contrast_Chinese = glue("1/{n_minus_1}*({other_Chinese}) - {Chinese}"),
				contrast_European = glue("1/{n_minus_1}*({other_European}) - {European}"),
				contrast_African = glue("1/{n_minus_1}*({other_African}) - {African}"),
				contrast_Hispanic= glue("1/{n_minus_1}*({`other_Hispanic or Latin American`}) - {`Hispanic or Latin American`}"),
			) |> 
			rename(
				this_Chinese = Chinese,
				this_European = European,
				this_African = African,
				this_Hispanic = `Hispanic or Latin American`
			) |> 
			select(tissue_harmonised, starts_with(c("contrast_", "this_"))) |>
			pivot_longer(cols = starts_with(c("contrast_", "this_")), names_to = c("type", "ethnicity"), names_sep = "_") |> 
			pivot_wider(names_from = type, values_from = value) |>
			filter(!is.na(this)) |> 
			select(tissue_harmonised, ethnicity, contrast) |> 
			unite("tissue_harmonised", c(tissue_harmonised, ethnicity)) |> 
			deframe( )
	) |> 
	filter(is_immune == "TRUE") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

differential_composition_ethnicity_relative = glue("{result_directory}/ethnicity_relative_FALSE.rds") |> readRDS()

count_significance_ethnicity_cell_type = 
	differential_composition_ethnicity_relative |>
	test_contrasts(
		c(
			African = "1/3*(`ethnicityHispanic or Latin American` + ethnicityEuropean + ethnicityChinese) - ethnicityAfrican",
			Hispanic = "1/3*(ethnicityAfrican + ethnicityEuropean + ethnicityChinese) - `ethnicityHispanic or Latin American`",
			European = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityChinese) - ethnicityEuropean",
			Chinese = "1/3*(ethnicityAfrican + `ethnicityHispanic or Latin American` + ethnicityEuropean) - ethnicityChinese"
		),
		test_composition_above_logit_fold_change = 0.1
	) |>
	
	filter(parameter %in% c("African",  "Hispanic", "European", "Chinese"))
		mutate(is_significnt = c_FDR<0.05) |> 
		count(is_significnt)

		

count_significance_ethnicity_cell_type_tissue = 
	differential_composition_ethnicity_relative |>
	
	# Find stats of random effect with groups
	test_contrasts(
		contrasts =
			differential_composition_ethnicity_relative |>
			filter(parameter |> str_detect("___ethnicity")) |>
			distinct(parameter) |>
			tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
			mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
			mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
			select(tissue_harmonised, contrast) |> 
			separate(tissue_harmonised, c("tissue_harmonised", "ethnicity"), sep="_") |> 
			add_count(tissue_harmonised) |> mutate(contrast = glue("({contrast})")) |>  
			pivot_wider(names_from = ethnicity, values_from = contrast) |> 
			
			# Filter where there are at least 3 ethnicities
			mutate(n_minus_1 = n-1) |> 
			filter(n_minus_1>1) |>
			
			# Build contrasts
			rowwise() |> 
			mutate(
				other_Chinese = c(European, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
				other_European = c(Chinese, African, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
				other_African = c(European, Chinese, `Hispanic or Latin American`) |> str_subset(".+") |> str_c(collapse=" + "),
				`other_Hispanic or Latin American` = c(European, African, Chinese) |> str_subset(".+") |> str_c(collapse=" + "),
			) |> 
			mutate(
				contrast_Chinese = glue("1/{n_minus_1}*({other_Chinese}) - {Chinese}"),
				contrast_European = glue("1/{n_minus_1}*({other_European}) - {European}"),
				contrast_African = glue("1/{n_minus_1}*({other_African}) - {African}"),
				contrast_Hispanic= glue("1/{n_minus_1}*({`other_Hispanic or Latin American`}) - {`Hispanic or Latin American`}"),
			) |> 
			rename(
				this_Chinese = Chinese,
				this_European = European,
				this_African = African,
				this_Hispanic = `Hispanic or Latin American`
			) |> 
			select(tissue_harmonised, starts_with(c("contrast_", "this_"))) |>
			pivot_longer(cols = starts_with(c("contrast_", "this_")), names_to = c("type", "ethnicity"), names_sep = "_") |> 
			pivot_wider(names_from = type, values_from = value) |>
			filter(!is.na(this)) |> 
			select(tissue_harmonised, ethnicity, contrast) |> 
			unite("tissue_harmonised", c(tissue_harmonised, ethnicity)) |> 
			deframe( )
	) |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)


# SEX
differential_composition_sex_absolute = glue("{result_directory}/sex_absolute_FALSE.rds") |> readRDS()

count_significance_sex_immune_load =
	differential_composition_sex_absolute |> 
	test_contrasts(		test_composition_above_logit_fold_change = 0.1	) |>
	filter(is_immune == "TRUE") |> 
	filter(parameter == "sexmale") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)


count_significance_sex_immune_load_tissue =
	differential_composition_sex_absolute |> 
	test_contrasts(
		contrasts =
			differential_composition_sex_absolute |>
			filter(parameter |> str_detect("___sex")) |>
			distinct(parameter) |>
			mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
			filter(contrast |> str_detect("_female", negate = TRUE)) |> 
			deframe( ),
		test_composition_above_logit_fold_change = 0.1
	) |>
	filter(is_immune == "TRUE") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

rm(differential_composition_sex_absolute)
gc()

differential_composition_sex_relative = glue("{result_directory}/sex_relative_FALSE.rds") |> readRDS()
 

count_significance_sex_cell_type = 
	differential_composition_sex_relative |>
	test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
	filter(parameter == "sexmale") |> 
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

count_significance_sex_cell_type_tissue = 
	differential_composition_sex_relative |>
	
	# Find stats of random effect with groups
	test_contrasts(
		contrasts =
			differential_composition_sex_relative |>
			filter(parameter |> str_detect("_male___sex")) |>
			distinct(parameter) |>
			mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
			deframe( ),
		test_composition_above_logit_fold_change = 0.1
	)  |>
	mutate(is_significnt = c_FDR<0.05) |> 
	count(is_significnt)

rm(differential_composition_sex_relative)
gc()

plot_significance_overall =
	count_significance_age_immune_load |>
	mutate(name = c(
		"count_significance_age_immune_load"
	)) |>
	bind_rows(
		count_significance_age_immune_load_tissue |>
			mutate(name = c(
				"count_significance_age_immune_load_tissue"
			))
	) |>
	bind_rows(count_significance_age_cell_type |>
							mutate(name = c(
								
								"count_significance_age_cell_type"
							))) |>
	bind_rows(count_significance_age_cell_type_tissue |>
							mutate(name = c(
								
								"count_significance_age_cell_type_tissue"
							))) |>
	
	bind_rows(count_significance_sex_immune_load |>
							mutate(name = c(
								
								"count_significance_sex_immune_load"
							))) |>
	bind_rows(count_significance_sex_immune_load_tissue |>
							mutate(name = c(
								
								"count_significance_sex_immune_load_tissue"
							))) |>
	bind_rows(count_significance_sex_cell_type |>
							mutate(name = c(
								
								"count_significance_sex_cell_type"
							))) |>
	bind_rows(count_significance_sex_cell_type_tissue |>
							mutate(name = c(
								
								"count_significance_sex_cell_type_tissue"
							))) |>
	
	bind_rows(count_significance_ethnicity_immune_load |>
							mutate(name = c(
								
								"count_significance_ethnicity_immune_load"
							))) |>
	bind_rows(count_significance_ethnicity_immune_load_tissue |>
							mutate(name = c(
								
								"count_significance_ethnicity_immune_load_tissue"
							))) |>
	bind_rows(count_significance_ethnicity_cell_type |>
							mutate(name = c(
								
								"count_significance_ethnicity_cell_type"
							))) |>
	bind_rows(count_significance_ethnicity_cell_type_tissue |>
							mutate(name = c(
								
								"count_significance_ethnicity_cell_type_tissue"
							))) |>
	tidyr::extract(name, c("factor", "variable", "resolution"), "count_significance_([a-zA-Z]+)_([a-zA-Z]+_[a-zA-Z]+)_?(.*)", remove = FALSE) |>
	mutate(resolution = if_else(resolution == "", "overall", resolution)) |>
	unite("xlab", c(variable, resolution), remove = FALSE) |>
	mutate(xlab = xlab |> fct_relevel(c("immune_load_overall", "immune_load_tissue", "cell_type_overall", "cell_type_tissue"))) |>
	with_groups(name, ~ .x |> mutate(sum_n = sum(n))) |>
	mutate(proportion = n/sum_n) |>
	mutate(factor = factor |> str_to_sentence()) |>
	ggplot(aes(xlab, proportion, fill=is_significnt)) +
	geom_bar(stat = "identity")+
	geom_text(aes(y = 0.5, label = sum_n), size = 2.5, angle=90) +
	facet_wrap( ~ factor,  nrow=1) +
	scale_fill_manual(values = c("FALSE"="grey", "TRUE"="#D5C711")) +
	ylab("Proportion of significant tests") +
	xlab("Hypotheses") +
	theme_multipanel +
	theme(axis.text.x = element_text(angle=20, hjust = 1, vjust = 1))
