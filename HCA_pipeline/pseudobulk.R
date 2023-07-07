library(tidySummarizedExperiment)
library(tidyverse)
library(tidybulk)

# watch out for this 5d53ce534b4d5d37418b5906cd706d84

# counts = 
# 	dir("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.3", full.names = TRUE) |>
# 	enframe(value="file") |>
# 	arrange(file) |>
# 	filter(file |> str_detect("rds$")) |>
# 	mutate(se = imap(	file, readRDS	)) |> 
# 	
# 	mutate(se = imap(
# 		se,
# 		~ .x |>
# 			select(
# 				.feature, .sample, counts,  sample_, 
# 				cell_type_harmonised, tissue_harmonised, 
# 				.aggregated_cells, sex, ethnicity, 
# 				age_days, assay, file_id, disease, mean_nFeature, median_nFeature, contains("simplified")
# 			) , 
# 		.progress = TRUE
# 	)) |>
# 	
# 	# # Establish the baseline for simplified ethnicity. European as it is the most represented
# 	# # This is so I have a tight intercept term for data simulation
# 	# mutate(ethnicity_simplified = case_when(
# 	# 	ethnicity %in% c("European", "Chinese", "African", "Hispanic or Latin American") ~ ethnicity,
# 	# 	TRUE ~ "Other"
# 	# )) |> 
# 	# mutate(
# 	# 	ethnicity_simplified = 
# 	# 		ethnicity_simplified |> 
# 	# 		fct_relevel(c("European", "Chinese", "African", "Hispanic or Latin American", "Other")
# 	# 		)) |> 
# 	# 
# 	# # Establish the baseline for simplified assay
# 	# # Summarise assays to get more stable data simulations 
# 	# # 10x as baseline
# 	# mutate(assay_simplified = if_else(assay |> str_detect("10x"), "10x", assay)) |> 
# 	# mutate(assay_simplified = factor(assay_simplified)) |>
# 	
# 	# Filter a minimum number of samples
# 	filter(map_int(se, ncol)>0) |>
# 	mutate(cell_type_DE  = map_chr(
# 		se,
# 		~ .x |> pull(cell_type_harmonised) |> unique(), 
# 		.progress = TRUE
# 	)) |> 
# 	mutate(tissue_harmonised_DE  = map_chr(
# 		se,
# 		~ .x |> pull(tissue_harmonised) |> unique(), 
# 		.progress = TRUE
# 	)) |>
# 	select(-name, -file) |>
# 	nest(data = -c(cell_type_DE, tissue_harmonised_DE)) |>
# 	mutate(data = map(data, ~ .x |> pull(se), .progress = TRUE))
# 
# counts |> saveRDS("immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.3/pseudobulk_for_fernando.rds")
# 
# counts = counts |> 
# 	mutate(data = map(data, ~ .x |> unnest_summarized_experiment(se), .progress = TRUE)) |>
# 	filter(map_int(data, ncol)>6) 
	

x = 
	dir("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4", full.names = TRUE) |>
	enframe(value = "file") |>
	extract(col = file, into = c("tissue", "cell_type"), regex = ".+pseudobulk_0.2.3.4/([a-zA-Z]+)__(.+)__TRUE\\.rds", remove = FALSE) |> 
	nest(cell_type_data = -cell_type)
	
	# Select only one cell type
	filter(cell_type=="b____memory") |>
	mutate(se = map(file, readRDS))

se = x |> 
	mutate(se = map(se, ~ .x |> select(-any_of("orig.ident"))))  |> 
	mutate(se = map(se, ~ .x |>  select(.aggregated_cells, .feature,median_nFeature, .sample, age_days, assay_simplified,cell_type_harmonised, counts, disease, ethnicity_simplified , file_id, sex, tissue_harmonised) )) |> 
	pull(se) %>% 
	reduce(cbind) |>
	
	# Filter based on cell
	filter(.aggregated_cells > 30) |>

	# Filter at least 5000 genes recorded
	filter(.sample %in% (
		se |>
			assay("counts") |> 
			apply(2, function(x) (x>0) |> which() |> length()) |>
			enframe() |>
			filter(value > 5000) |>
			pull(name)
	))

se = 
	se |>
	mutate(age_days = scale(age_days) |> as.numeric()) |>
	mutate(sex_etnicity = glue("{sex}__{ethnicity_simplified}")) |>
	identify_abundant(factor_of_interest = "sex") |>
	quantile_normalise_abundance() 

se |> saveRDS("~/temp.rds")

# se = 
	se[sample(1:2000, size = 100),] |>
	filter(tissue_harmonised != "blood") |>
	test_differential_abundance(
		~ 1 + sex + ethnicity_simplified + assay_simplified + age_days + file_id + ( 1 + sex + ethnicity_simplified | tissue_harmonised ), 
		.abundance = "counts_scaled",
		method = "glmmseq_lme4" , 
		cores=10
	)
 
