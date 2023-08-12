
library(targets)
library(glue)
result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.4"




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

