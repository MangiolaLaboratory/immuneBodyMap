library(tidyverse)
library(forcats)
library(CuratedAtlasQueryR)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
library(tidyHeatmap)
library(ComplexHeatmap)
library(scales)
library(ggplot2)
library(targets)
library(tidybulk)
library(tidySummarizedExperiment)

home = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s"
result_directory = glue("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age")


# Calculate softmax from an array of reals
softmax <- function (x) {
	logsumexp <- function (x) {
		y = max(x)
		y + log(sum(exp(x - y)))
	}

	exp(x - logsumexp(x))
}

# This function is used to format ggplots to save space
dropLeadingZero <-
	function(l) {
		stringr::str_replace(l, '0(?=.)', '')
	}

# Scale axis for ggplot
S_sqrt <- function(x) {
	sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
	x ^ 2 * sign(x)
}
S_sqrt_trans <-
	function()
		scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

# This function shorten names of cell types for visualisation purposes
clean_names = function(x){
	x |>  mutate(
		tissue_harmonised =
			tissue_harmonised |>
			str_remove("tissue_harmonised") |>
			str_replace_all("_", " ") |>
			str_replace("gland", "gld") |>
			str_replace("node", "nd") |>
			str_replace("skeletal", "sk")
	)
}


# Read inout files
data_for_immune_proportion_absolute_file = glue("{result_directory}/input_absolute.rds")
data_for_immune_proportion = readRDS(data_for_immune_proportion_absolute_file)
data_for_immune_proportion_relative_file = glue("{result_directory}/input_relative.rds")
data_for_immune_proportion_relative = readRDS(data_for_immune_proportion_relative_file)

# Color coding for tissue
tissue_color =
	data_for_immune_proportion_relative |>
	distinct(tissue_harmonised ) |>
	arrange(tissue_harmonised) |>
	mutate(color = dittoSeq::dittoColors()[1:n()]) |>
	deframe()

# Load cell type colors
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/f0b6bf9f59847c8b9f0a638262a6b8dd697affb7/color_cell_types.R")
cell_type_color =
	data_for_immune_proportion |>
	pull(cell_type_harmonised) |>
	unique() |>
	get_cell_type_color() %>%
	c(macro = .["macrophage"] |> as.character())

# Threshold that equates to a linear increase of 1% from 20% to 21%
# This convoluted threshold gives a lay meaning to a increase in the softmax space
# which is not linear and hard to grasp the meanin of a lay audience
# We use this threshold to test for significant effect bigger that it
FDR_threshold_1_percent_change_at_20_percent_baseline = 0.017

#------------------------------#
# Ethnicity analyses for immune cellularity
#------------------------------#

# Read results
differential_composition_ethnicity_absolute_file = glue("{result_directory}/ethnicity_absolute_FALSE.rds")
proportions_ethnicity_absolute_file = glue("{result_directory}/ethnicity_absolute_FALSE_proportion_adjusted.rds")
differential_composition_ethnicity_absolute = readRDS(differential_composition_ethnicity_absolute_file)
proportions_ethnicity_absolute_adjusted = readRDS(proportions_ethnicity_absolute_file)

# # save csv for SUPPLEMENTARY
# differential_composition_ethnicity_absolute |>
# 	sccomp_test(test_composition_above_logit_fold_change = 0.1) |>
# 	select(-count_data) |>
# 	write_csv("sccomp_on_HCA_0.2.3.4/SUPPLEMENTARY_ethnicity_cellularity_estimates.csv")


# Plot effect of composition and variability
# For immune cellularity (proportion of immune cells)
# Wtih uncertainty
plot_ethnicity_absolute_1D =
	differential_composition_ethnicity_absolute   |>
	  sccomp_remove_unwanted_variation(~ ethnicity_simplified) |>
	  left_join(
	    data_for_immune_proportion_relative |> distinct(sample_, ethnicity_simplified)
	  ) |>
	  filter(!is.na(ethnicity_simplified)) |>
	  filter(is_immune == "TRUE") |>
	  with_groups(ethnicity_simplified, ~ .x |> mutate(median_proportion = median(adjusted_proportion))) |>
	  mutate(ethnicity_simplified = ethnicity_simplified |> fct_reorder(median_proportion, .desc = TRUE)) |>
	  with_groups(ethnicity_simplified, ~ .x |> mutate(sample_dummy = 1:n())) |>
	  select(ethnicity_simplified, adjusted_proportion, sample_dummy) |>

	  # I dont have many samples
	  filter(ethnicity_simplified != "Chinese") |>
	  pivot_wider(names_from = ethnicity_simplified, values_from = adjusted_proportion) |>

	  drop_na() |>

	  select(sample_dummy,    Other,`Hispanic or Latin American`, African, European,) |>
	  tidybulk::as_matrix(rownames = sample_dummy) |>
	  bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
	  coord_flip() +
	  xlab("Immune cellularity adjusted") +
	  theme_multipanel +
	  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))




# Create dataset to create the mannequin heatmap of the
# Tissues with differential immune cellularity
ethnicity_absolute_organ_tissue =
	differential_composition_ethnicity_absolute |>
	sccomp_test(
		contrasts =
			differential_composition_ethnicity_absolute |>
			filter(parameter |> str_detect("ethnicity") ) |>
			distinct(parameter) |>
			tidyr::extract( parameter, "ethnicity", "ethnicity_simplified(.+)___", remove = FALSE) |>
			filter(!is.na(ethnicity)) |>
			mutate(contrast = glue("`ethnicity_simplified{ethnicity}`  + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue", ".+___(.+)", remove = FALSE) |>
			unite("name", ethnicity, tissue) |>
			select(-parameter) |>
			deframe( ),
		test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
	) |>
	filter(is_immune == "TRUE")

# # save csv for SUPPLEMENTARY
# differential_composition_ethnicity_absolute |>
# 	sccomp_test(
# 		contrasts =
# 			differential_composition_ethnicity_absolute |>
# 			filter(parameter |> str_detect("sexmale___")) |>
# 			distinct(parameter) |>
# 			mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
# 			tidyr::extract(parameter, "tissue_harmonised", ".+___(.+)") |>
# 			filter(contrast |> str_detect("_female", negate = TRUE)) |>
# 			deframe( ),
# 		test_composition_above_logit_fold_change = 0.1
# 	) |>
# 	select(-count_data) |>
# 	select(1, 2, 4, 5, 6, 7, 8) |>
# 	write_csv("sccomp_on_HCA_0.2.3.5/SUPPLEMENTARY_ethnicity_cellularity_tissue_estimates_contrasts.csv")


# Draw the color palette for the mannequin heatmap of the
# Tissues with differential immune cellularity
colors_palette_for_organ_abundance =
	ethnicity_absolute_organ_tissue |>
	#filter(parameter |> str_detect("Other", negate = TRUE)) |>
	filter(c_FDR<0.05) |>
	select(parameter, c_effect, c_FDR) |>
	mutate(color = circlize::colorRamp2(
		seq(1.45,-1.45, length.out = 11),
		RColorBrewer::brewer.pal(11, "RdBu")
	)(c_effect)) |>
	mutate(rgb = map_chr(
		color,
		~ .x |>
			col2rgb() |>
			paste(collapse = " ")
	)) |>
	arrange(parameter) |>
	pull(color) |>
	scales::show_col(	cex_label = 0.5	)


#
# # Draw the boxplot for significant changes of the
# # Tissues with differential immune cellularity
# plot_ethnicity_absolute_organ_boxoplot_adjusted =
#   proportions_ethnicity_absolute_adjusted |>
#   left_join(
#     data_for_immune_proportion |>
#       distinct(sample_, tissue_harmonised, ethnicity, sex,tissue, file_id)
#   ) |>
#   inner_join(
#     ethnicity_absolute_organ_tissue |>
#       filter(c_FDR<0.07) |>
#       separate(parameter, c("tissue_harmonised", "sex"), sep="_") |>
#       distinct(tissue_harmonised)
#   ) |>
#
#   mutate(tissue_harmonised = tissue_harmonised |> str_to_sentence()) |>
#   filter(is_immune =="TRUE") |>
#   ggplot(aes(sex, adjusted_proportion )) +
#   geom_boxplot(outlier.shape = NA, lwd = 0.2, fatten = 0.2) +
#   geom_jitter(aes(color = file_id), width = 0.1, size=0.1) +
#   facet_wrap(~ tissue_harmonised, scale="free_x", nrow = 1) +
#   guides(color = "none") +
#   ylab("Adjusted proportion") +
#   theme_multipanel +
#   theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))



#------------------------------#
# Ethnicity analyses for immune composition
#------------------------------#

differential_composition_ethnicity_relative_file = glue("{result_directory}/ethnicity_relative_FALSE.rds")
proportions_ethnicity_relative_file = glue("{result_directory}/ethnicity_relative_FALSE_proportion_adjusted.rds")
differential_composition_ethnicity_relative = readRDS(differential_composition_ethnicity_relative_file)
proportions_ethnicity_relative_adjusted = readRDS(proportions_ethnicity_relative_file)
gene_chr = read_csv("symbol_chr.csv")

# # save csv for SUPPLEMENTARY
# differential_composition_ethnicity_relative |>
# 	sccomp_test(test_composition_above_logit_fold_change = 0.1) |>
# 	select(-count_data) |>
# 	write_csv("sccomp_on_HCA_0.2.3.4/SUPPLEMENTARY_ethnicity_composition_estimates.csv")

# # Test difference in immune composition at the body
# # Comparing each ethnicity with the average of the others
# data_for_ethinicity_relative_plot =
# 	differential_composition_ethnicity_relative |>
# 	sccomp_test(
# 		c(
# 			european = "`(Intercept)`",
# 			asian = "`(Intercept)` + ethnicity_simplifiedChinese",
# 			african = "`(Intercept)` + ethnicity_simplifiedAfrican",
# 			hispanic = "`(Intercept)` + `ethnicity_simplifiedHispanic or Latin American`",
# 			other = "`(Intercept)` + ethnicity_simplifiedOther"
# 		),
# 		test_composition_above_logit_fold_change = 0.1
# 	) |>
#
# 	filter(parameter %in% c("european", "asian", "african", "hispanic", "other"))
#
#
#
# # Plot volcano of the differences at the immune composition level
# # (for each cell type)
# volcano_relative_ethnicity =
# 	data_for_ethinicity_relative_plot |>
# 	filter(cell_type_harmonised != "non_immune") |>
# 	mutate(naive_experienced = case_when(
# 		cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
# 		cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
# 	)) |>
# 	mutate(significant = c_FDR<0.05) |>
# 	mutate(cell_type_harmonised = case_when(c_FDR<1e-2~cell_type_harmonised)) |>
# 	ggplot(aes(c_effect, c_FDR)) +
# 	geom_vline(xintercept = -0.1, color="grey", linetype = "dashed", size=0.2) +
# 	geom_vline(xintercept = 0.1, color="grey", linetype = "dashed", size=0.2) +
# 	geom_hline(yintercept = 0.05, color="grey", linetype = "dashed", size=0.2) +
#
# 	geom_point(aes(color=naive_experienced, size=significant)) +
# 	ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
# 	facet_wrap(~parameter, nrow=1) +
# 	scale_y_continuous(trans = tidybulk::log10_reverse_trans()) +
# 	scale_x_continuous(trans="S_sqrt") +
# 	scale_color_brewer(palette="Set1", na.value = "grey50") +
# 	scale_size_discrete(range = c(0, 0.5)) +
# 	ylab("False-discovery rate") +
# 	xlab("Effect from baseline (average ethnicites)") +
# 	theme_multipanel




# Scatter plot of the significant cell types which compositon change through ageing
# The proportions are adjusted to exclude other effects including
# Sex, ethnicity, random effects (e.g. datasets) and technology

# outliers_df =
#   differential_composition_ethnicity_relative |>
#   select(cell_type_harmonised, count_data) |>
#   unnest(count_data) |>
#   distinct() |>
#   filter(outlier) |>
#   select(cell_type_harmonised, sample_)




# # Plot for presentation B memory across tissues
# res_relative_proportions_ethnicity_tissue =
#   differential_composition_ethnicity_relative |>
#   sccomp_remove_unwanted_variation(~ 1 + ethnicity_simplified + ( 1 + ethnicity_simplified | tissue_harmonised ), ~ ethnicity_simplified)  |>
#   inner_join(data_for_immune_proportion_relative |> tidybulk::pivot_sample(sample_) )
#
# res_relative_proportions_ethnicity_tissue |> saveRDS(glue("{result_directory}/res_relative_proportions_ethnicity_tissue.rds"))

res_relative_proportions_ethnicity_tissue = readRDS(glue("{result_directory}/res_relative_proportions_ethnicity_tissue.rds"))


# # Get trend line to be used in the scatter plot of global compositional changes, plot_ethnicity_relative
# line_ethnicity_relative_mean_b_memory_per_tisue =
#
#   differential_composition_ethnicity_relative |>
#   filter(cell_type_harmonised=="b memory") |>
#   nest(data = -cell_type_harmonised) |>
#
#   # Add tissue
#   mutate(
#     tissue_harmonised =
#       list(
#         differential_composition_ethnicity_relative |>
#           filter(factor=="tissue_harmonised") |>
#           distinct(parameter) |>
#           pull(1) |>
#           str_remove("tissue_harmonised")
#       )
#   ) |>
#   unnest(tissue_harmonised) |>
#
#   mutate(x = list(seq(-3, 3, by = 0.1))) |>
#   mutate(y = pmap(list(data, x, tissue_harmonised), ~ {
#     ( ..2 * ..1 |> filter(parameter == "sexmale") |> pull(c_effect))  + # SLOPE
#       (..1 |> filter(parameter == "(Intercept)") |> pull(c_effect)) + # INTERCEPT
#       (..1 |> filter(parameter == glue("tissue_harmonised{..3}")) |> pull(c_effect)) + # INTERCEPT TISSUE
#       (..2 * ..1 |> filter(parameter == glue("{..3}___ethnicity_days")) |> pull(c_effect))   # GROUP-LEVEL SLOPE
#
#   })) |>
#   dplyr::select(-data) |>
#   unnest(c(x, y)) |>
#   with_groups(x, ~ .x |> mutate(proportion = softmax(y))) |>
#   mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
#   filter(x_corrected |> between(30.0 , 30295.0))
#


# res_relative_proportions_ethnicity_tissue |>
#
#   filter(cell_type_harmonised=="cd14 mono") |>
#   ggplot(aes(sexmale, adjusted_proportion,colour = tissue_harmonised)) +
#   geom_point(
#     aes(fill = tissue_harmonised),
#     shape = 21,
#     stroke = 0,
#     size = 0.6
#   ) +
#   # geom_line(
#   # 	aes(x_corrected, proportion, color = tissue_harmonised),
#   # 	data =
#   # 		line_ethnicity_relative_mean_b_memory_per_tisue
#   # ) +
#   geom_smooth(
#     method = "glm",
#     method.args = list(family = "binomial"),
#     se = FALSE, size=1) +
#   scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
#   scale_x_continuous(
#     labels = function(x)
#       round(x / 356)
#   ) +
#   scale_fill_manual(values = tissue_color) +
#   scale_color_manual(values = tissue_color) +
#   xlab("Years") +
#   ylab("Adjusted proportions") +
#   guides(fill = "none", color = "none") +
#   theme_multipanel

# res_relative_proportions_ethnicity_tissue |>




# # Fold change for cell type by organ
# # These statistics are used in the paper result section
# differential_composition_ethnicity_relative |>
#
#   print_estimate_plus_minus_relative(
#     contrasts_baseline =
#       differential_composition_ethnicity_relative |>
#       filter(parameter |> str_detect("sexmale___")) |>
#       tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
#       mutate(parameter = glue("{parameter}")) |>
#       distinct(parameter, tissue_harmonised) |>
#       mutate(contrast = glue("`(Intercept)` + `tissue_harmonised{tissue_harmonised}`") |> as.character()) |>
#       select(parameter, contrast) |>
#       filter(parameter |> str_detect("adipose", negate = TRUE)) |>
#       deframe(),
#
#     contrasts =
#       differential_composition_ethnicity_relative |>
#       filter(parameter |> str_detect("___ethnicity_days")) |>
#       tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
#       distinct(parameter, tissue_harmonised) |>
#       mutate(contrast = glue("sexmale  + `{parameter}`") |> as.character()) |>
#       mutate(contrast = glue("`(Intercept)` + `tissue_harmonised{tissue_harmonised}` + (1.73 * ({contrast}))")) |>
#       mutate(parameter = glue("__{parameter}")) |>
#       select(parameter, contrast) |>
#       filter(parameter |> str_detect("adipose", negate = TRUE)) |>
#
#
#       deframe(),
#
#     contrasts_uncertainty =
#       differential_composition_ethnicity_relative |>
#       filter(parameter |> str_detect("___ethnicity_days")) |>
#       distinct(parameter) |>
#       mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
#       tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE)  |>
#       filter(parameter |> str_detect("adipose", negate = TRUE)) |>
#       mutate(parameter = glue("__{parameter}")) |>
#
#       deframe()
#   )

# Prepare the dataset for drawing the heatmap of changes by tissues and celltypes, Figure 3
ethnicity_tissue_sample_size =
	readRDS(glue("{result_directory}/input_relative.rds")) |>
	distinct(ethnicity, ethnicity_simplified, tissue_harmonised, sample_) |>
	dplyr::count(ethnicity_simplified, tissue_harmonised)

FDR_non_european =
	differential_composition_ethnicity_relative |>
	sccomp_test(
		contrasts =
			differential_composition_ethnicity_relative |>
			filter(parameter |> str_detect("ethnicity") ) |>
			distinct(parameter) |>
			tidyr::extract( parameter, "ethnicity", "ethnicity_simplified(.+)___", remove = FALSE) |>
			filter(!is.na(ethnicity)) |>
			mutate(contrast = glue("`ethnicity_simplified{ethnicity}`  + `{parameter}`") |> as.character()) |>
			tidyr::extract(parameter, "tissue", ".+___(.+)", remove = FALSE) |>
			filter(ethnicity != "Other") |>
			filter(tissue != "nose") |>
			unite("name", ethnicity, tissue, remove = FALSE, sep="___") |>
			select(name, contrast) |>
			deframe()
		#,
		#test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
	) |>
	separate(parameter, c("ethnicity", "tissue"), sep="___", remove = FALSE) |>
	select(ethnicity, tissue, c_FDR, c_effect, cell_type = cell_type_harmonised)



contrasts_df_non_european =
	differential_composition_ethnicity_relative |>
	filter(parameter |> str_detect("ethnicity") ) |>
	distinct(parameter) |>
	tidyr::extract( parameter, "ethnicity", "ethnicity_simplified(.+)___", remove = FALSE) |>
	filter(!is.na(ethnicity)) |>
	tidyr::extract(parameter, "tissue", ".+___(.+)", remove = FALSE) |>
	mutate(contrast = glue("`(Intercept)` + `ethnicity_simplified{ethnicity}`  + `(Intercept)___{tissue}` + `{parameter}`") |> as.character()) |>
	tidyr::extract(parameter, "tissue", ".+___(.+)", remove = FALSE) |>
	filter(ethnicity != "Other") |>
	filter(tissue != "nose") |>
	unite("name", ethnicity, tissue, remove = FALSE, sep="___")

contrasts_european =
	differential_composition_ethnicity_relative |>
	filter(parameter |> str_detect("\\(Intercept\\)___") ) |>
	distinct(parameter) |>
	mutate(ethnicity = "European") |>
	filter(!is.na(ethnicity)) |>
	mutate(contrast = glue("`(Intercept)`  + `{parameter}`") |> as.character()) |>
	tidyr::extract(parameter, "tissue", ".+___(.+)", remove = FALSE) |>
	filter(tissue %in% (contrasts_df_non_european |> pull(tissue))) |>
	unite("name", ethnicity, tissue, sep="___") |>
	select(name, contrast) |>
	deframe()

df_heatmap_ethnicity_relative_organ_cell_type =

	differential_composition_ethnicity_relative |>

	# Find stats of random effect with groups
	sccomp_test(
		contrasts =
			contrasts_df_non_european |>
			select(name, contrast) |>
			deframe( ) |>
			c(contrasts_european),
		test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
	)  |>

	add_count(cell_type_harmonised) |>
	arrange(parameter, desc(n)) |>

	dplyr::rename(tissue = parameter) |>
	dplyr::rename(cell_type = cell_type_harmonised) |>

	filter(!cell_type %in% c("non_immune")) |>

	# To be fixed in the model
	mutate(is_treg = cell_type =="treg") |>
	nest(data = -is_treg) |>
	mutate(data = map2(
		data, is_treg,
		~ {
			if(.y) .x |> mutate(c_effect = c_effect/10 )
			else(.x)
		}
	)) |>
	unnest(data) |>

	# Cell type abundance
	with_groups(cell_type, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(cell_type_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

	# Tissue diversity
	with_groups(tissue, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(tissue_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

	# First rank
	with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

	dplyr::rename(`Mean diff` = cell_type_mean_change) |>
	mutate(`Mean diff tissue` = -tissue_mean_change) |>

	# Color
	left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
	left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  )  |>

	# Counts
	left_join(
		data_for_immune_proportion_relative |>
			dplyr::count(tissue_harmonised, ethnicity_simplified, name = "count_tissue") |>
			unite("tissue", ethnicity_simplified, tissue_harmonised, sep="___") |>
			mutate(count_tissue = log(count_tissue))
	) |>


	# Join_intercept
	left_join(

		differential_composition_ethnicity_relative |>
			select(cell_type_harmonised, count_data) |>
			unnest(count_data) |>
			distinct() |>
			with_groups(sample_, ~ .x |> mutate(proportion = count/sum(count))) |>
			with_groups(c(tissue_harmonised, cell_type_harmonised), ~ .x |> summarise(mean_proportion = mean(proportion))) |>
			mutate(mean_proportion = pmax(mean_proportion, 1e-7)) |>
			mutate(mean_proportion_logit = boot::logit(mean_proportion)) |>
			dplyr::rename(tissue = tissue_harmonised, cell_type = cell_type_harmonised) |>
			mutate(mean_proportion_logit = (mean_proportion_logit - min(mean_proportion_logit))/4)

	) |>

	# Format
	dplyr::rename(tissue_ethnicity = tissue) |>
	mutate(tissue_ethnicity = tissue_ethnicity |> str_remove(" or Latin American")) |>
	separate(tissue_ethnicity, c("ethnicity", "tissue"), sep="___", remove = FALSE) |>

	# Add stats
	select(-c_FDR) |>
	left_join(FDR_non_european |> dplyr::rename(difference = c_effect)) |>
	replace_na(list(c_FDR = 1)) |>

	# Shorten names
	mutate(cell_type = cell_type |>
				 	str_replace("macrophage", "macro") |>
				 	str_replace("megakaryocytes", "mega") |>
				 	str_remove("phage") |>
				 	str_replace("th1/th17", "th1/17") |>
				 	str_replace("mono", "mn") |>
				 	str_replace("tcm", "cm") |>
				 	str_replace("cd4 th", "Th") |>
				 	str_replace("memory", "mem") |>
				 	str_replace("naive", "nv") |>
				 	str_replace("terminal effector cd4 t", "cd4 eff") |>
				 	str_remove("cyte")
	) |>
	mutate(tissue = tissue |>
				 	str_replace("intestine", "int") |>
				 	str_replace("large", "lrg") |>
				 	str_replace("small", "sml") |>
				 	str_replace("node", "nd") |>
				 	str_replace("prostate", "prost")
	) |>

	# Order
	mutate(tissue = fct_reorder(tissue, `Mean diff tissue`)) |>
	mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))

# Plot heatmap with tidyHeatmap, Figure 3
plot_heatmap_ethnicity_relative_organ_cell_type =

	df_heatmap_ethnicity_relative_organ_cell_type |>


	with_groups(tissue_ethnicity, ~ .x |> mutate(proportion = softmax(c_effect))) |>
	mutate(proportion_label = proportion |> round(3) |> dropLeadingZero())  |>

	# Filter
	filter(!cell_type %in% c( "dnt", "thymo", "mega" )) |>
	filter(cell_type != "immune_unclassified") |>

	# Filte relevant tissues that have enough samples for non europeans
	filter(tissue %in% (
		ethnicity_tissue_sample_size |>
			filter(!ethnicity_simplified %in% c("European", "Other")) |>
			filter(n>10) |>
			pull(tissue_harmonised) |> unique()
	)) |>

	left_join(
		ethnicity_tissue_sample_size |>
			select(tissue = tissue_harmonised, ethnicity = ethnicity_simplified)
	) |>


	group_by(tissue) |>
	mutate(ethnicity = ethnicity |> fct_relevel(c("European", "African", "Chinese", "Hispanic or Latin American"))) |>
	arrange(ethnicity) |>
	mutate(tissue_ethnicity  = tissue_ethnicity |> fct_reorder(as.integer(ethnicity))) |>

	# Heatmap
	heatmap(
		tissue_ethnicity, cell_type, proportion,
		palette_value = circlize::colorRamp2(
			c(seq(0, 0.03, length=5), 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
			c(
				viridis::magma(5),
				viridis::viridis(7) |> rev()
			)
		),
		cluster_columns = FALSE,
		cluster_rows = FALSE,
		row_names_gp = gpar(fontsize = 6),
		column_names_gp = gpar(fontsize = 6),
		column_title_gp = gpar(fontsize = 0),
		row_title_gp = gpar(fontsize = 0),
		show_heatmap_legend = FALSE,
		#transform = car::logit,
		clustering_distance_rows = "manhattan",
		clustering_distance_columns = "manhattan",
		clustering_method_rows = "ward.D",
		clustering_method_columns = "ward.D",
		column_dend_height = unit(0.3, "cm"),
		row_dend_width = unit(0.5, "cm")
	) |>

	annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
	#annotation_bar(`Mean diff tissue`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
	# annotation_tile(
	#   tissue, show_legend = FALSE,
	#   palette =
	#     df_heatmap_ethnicity_relative_organ_cell_type |>
	#     distinct(tissue, tissue_color) |>
	#     arrange(tissue) |>
	#     deframe(),
	#   size = unit(0.2, "cm")
	# ) |>
	annotation_tile(
		cell_type, show_legend = FALSE,
		palette =
			df_heatmap_ethnicity_relative_organ_cell_type |>
			distinct(cell_type, cell_type_color)  |>
			arrange(cell_type) |>
			deframe(),
		size = unit(0.2, "cm")
	) |>

	annotation_tile(
		count_tissue, show_legend = FALSE,
		size = unit(0.2, "cm"),
		palette = c( "white", "black")
	) |>
	layer_arrow_down(c_FDR < 0.05 & n > 5 & difference < 0, .size = 1.5) |>
	layer_arrow_up(c_FDR < 0.05 & n > 5 & difference > 0, .size = 1.5) #, .size = sqrt(mean_proportion  ) * 5 )

# Save heatmap separately
plot_heatmap_ethnicity_relative_organ_cell_type |>
	save_pdf(
		filename = glue("{result_directory}/plot_heatmap_ethnicity_relative_organ_cell_type.pdf"),
		width = 80*1.5, height = 60*1.5, units = "mm"
	)


# DE
library(targets)
result_directory_pseudobulk = "pseudobulk_0.2.3.5_non_immune"
store = glue("{result_directory_pseudobulk}/_targets__pseudobulk_non_immune_split3")
# Plot of importance of composition vs transcription


# de_ethnicity_cell_type =
# 	tar_meta( store = store	) |>
# 	dplyr::filter(name |> str_detect("estimates_ethnicity_cell_type_")) |>
# 	filter(!is.na(data)) |>
# 	mutate(se = map(
# 		name,
# 		~ .x |>
# 			tar_read_raw(store=store ) |>
# 			mutate(is_immune = map_lgl(data, ~ .x |> tidySummarizedExperiment::pull(cell_type_harmonised) %in% c(
# 				"b memory",     "cd8 tcm"  ,    "cd8 tem",      "plasma" ,
# 				"b naive" ,     "cd14 mono" ,   "cd4 fh"   ,    "cd4 naive"  ,
# 				"cd4 th1/th17", "cd4 th17",     "cd4 th2" ,     "cdc"  ,
# 				"ilc" ,         "macrophage",   "mait" ,        "nk" ,
# 				"tgd"  ,        "treg"
# 			) |> any())) |>
# 			mutate(data = map(data, tidybulk::pivot_transcript)),
# 		.progress=T
# 		# ,
# 		# .env_globals = environment()
# 	)) |>
# 	select(-data) |>
# 	unnest(se)
# de_ethnicity_cell_type |> saveRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_ethnicity_cell_type.rds")

de_ethnicity_cell_type = readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_ethnicity_cell_type.rds")

# de_ethnicity =
# 	tar_meta(store = glue("{result_directory_de}/_targets__ethnicity"), starts_with("data_")) |>
# 	filter(!is.na(data)) |>
# 	mutate(se = map(
# 		name,
# 		~ tar_read_raw(.x, store = glue("{result_directory_de}/_targets__ethnicity")) |>
# 			pivot_transcript(),
# 		.progress=T
# 	))
# de_ethnicity |> saveRDS(glue("{result_directory}/de_ethnicity.rds"))


rank_de_cell_type =
	readRDS(glue("{result_directory}/de_ethnicity.rds")) |>
	unnest(se) |>
	dplyr::count(name, P_ethnicity_simplified_adjusted < 0.05) |>
	drop_na() |>
	spread(`P_ethnicity_simplified_adjusted < 0.05`, n) |>
	mutate(proportion_significant = `TRUE` / (`FALSE` + `TRUE`)) |>
	arrange(desc(proportion_significant)) |>
	dplyr::rename(cell_type = name) |>

	# Parse cell type
	mutate(cell_type =
				 	cell_type |>
				 	str_remove("data_") |>
				 	str_replace_all("__", " ") |>
				 	str_replace("th1 th17", "th1/th17")
	) |>

	filter(!cell_type %in% c("immune_unclassified", "dnt")) |>
	rowid_to_column("rank_de")



plot_ranks_cell_type_barplot =
	rank_de_cell_type |>
	ggplot(aes(`TRUE`, cell_type |> fct_reorder(proportion_significant))) +
	geom_bar(aes(fill = cell_type), stat = "identity") +
	scale_fill_manual(values = cell_type_color) +
	guides(fill="none") +
	xlab("Number of significant genes") +
	theme_multipanel +
	theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_ranks_cell_type =

		differential_composition_ethnicity_relative |>
			sccomp_test(contrasts = c("non_european" = "1/3*(ethnicity_simplifiedChinese + ethnicity_simplifiedAfrican + `ethnicity_simplifiedHispanic or Latin American`)")) |>
			filter(parameter=="non_european") |>
			arrange(c_FDR) |>
			filter(!cell_type_harmonised %in% c("immune_unclassified", "dnt")) |>
			rowid_to_column("rank_composition") |>
			select(cell_type = cell_type_harmonised, rank_composition) |>

			full_join(

				# DE
				rank_de_cell_type ,
				by = join_by(cell_type)

			) |>
			select(cell_type, contains("rank")) |>
			pivot_longer(contains("rank")) |>
			mutate(name = name |> str_remove("rank")) |>

			# Parse cell types
			mutate(cell_type_formatted = cell_type |>
						 	str_replace("megakaryocytes", "mega") |>
						 	str_remove("phage") |>
						 	str_replace("th1 th17", "th1/17") |>
						 	str_replace("mono", "mn") |>
						 	str_replace("tcm", "cm") |>
						 	str_replace("cd4 th", "Th") |>
						 	str_replace("memory", "mem") |>
						 	str_replace("naive", "nv") |>
						 	str_replace("terminal effector cd4 t", "cd4 eff") |>
						 	str_remove("cyte")
			) |>


			# Plot
			ggplot(aes(name, value, group=cell_type,  label=cell_type_formatted)) +
			geom_line(aes(color = cell_type)) +
			geom_point(aes(color = cell_type)) +
			geom_text(size = 2.5) +
			scale_color_manual(values = cell_type_color) +
			scale_y_reverse() +
			guides(color = "none") +
			ylab("Rank of tissue with highest significance") +
			theme_multipanel +
			theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())



de_ethnicity_tissue =
	tar_meta( store = store	) |>
	dplyr::filter(name |> str_detect("estimates_ethnicity_tissue_")) |>
	filter(!is.na(data)) |>
	mutate(se = future_map(
		name,
		~ .x |>
			tar_read_raw(store=store ) |>
			mutate(is_immune = map_lgl(data, ~ .x |> tidySummarizedExperiment::pull(cell_type_harmonised) %in% c(
				"b memory",     "cd8 tcm"  ,    "cd8 tem",      "plasma" ,
				"b naive" ,     "cd14 mono" ,   "cd4 fh"   ,    "cd4 naive"  ,
				"cd4 th1/th17", "cd4 th17",     "cd4 th2" ,     "cdc"  ,
				"ilc" ,         "macrophage",   "mait" ,        "nk" ,
				"tgd"  ,        "treg"
			) |> any())) |>
			mutate(data = map(data, tidybulk::pivot_transcript)),
		.progress=T
		# ,
		# .env_globals = environment()
	)) |>
	select(-data) |>
	unnest(se)
de_ethnicity_tissue |> saveRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_ethnicity_tissue.rds")

# de_ethnicity_tissue_non_immune =
# 	tar_meta(store = glue("pseudobulk_0.2.3.5_non_immune/_targets__pseudobulk_non_immune"), starts_with("estimates_")) |>
# 	filter(!is.na(data)) |>
# 	mutate(se = map(
# 		name,
# 		~ {
# 			se = tar_read_raw(.x, store = glue("pseudobulk_0.2.3.5_non_immune/_targets__pseudobulk_non_immune"))
#
# 			se |>
# 				pivot_transcript() |>
# 				mutate(tissue_harmonised = se |> pull(tissue_harmonised) |> unique())
# 			},
# 		.progress=T
# 	))
# de_ethnicity_tissue_non_immune |> saveRDS("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue_non_immune.rds")

de_ethnicity_tissue_non_immune =
	readRDS("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue_non_immune.rds") |>
	mutate(tissue = map_chr(se, ~ .x |> pull(tissue_harmonised) |> unique())) |>
	filter(!is.na(tissue)) |>
	filter(map_lgl(se, ~ "P_ethnicity_simplified_adjusted" %in% colnames(.x), .progress = TRUE)) |>
	mutate(se = map(se, ~ .x |> select(.feature, P_ethnicity_simplified_adjusted)))

de_ethnicity_tissue =
	readRDS("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue.rds") |>
	dplyr::rename(tissue = name) |>
	mutate(tissue = tissue |> str_remove("data_")) |>
	# Parse cell type
	mutate(tissue = case_when(
		tissue == "node" ~ "lymph node",
		tissue == "large" ~ "intestine large",
		tissue == "small" ~ "intestine small",
		tissue == "gland" ~ "adrenal gland",
		TRUE ~ tissue
	)) |>
	filter(map_lgl(se, ~ "P_ethnicity_simplified_adjusted" %in% colnames(.x))) |>
	mutate(se = map(se, ~ .x |> select(.feature, P_ethnicity_simplified_adjusted)))

rank_de_tissue =
	de_ethnicity_tissue |>

	left_join(
		de_ethnicity_tissue_non_immune |>
			select(tissue, se_non_immune = se)
	) |>

	mutate(n = map2(
		se,se_non_immune,
		~ {
			.x = .x |>
			filter(!.feature %in% gene_chr$ID)

			if(!is.null(.y)) 	.x = .x |>	filter(!.feature %in% (.y |> pull(.feature)))
			.x |>
				dplyr::count( P_ethnicity_simplified_adjusted < 0.05)

		}
	)) |>
	unnest(n) |>
	filter(!is.na(`P_ethnicity_simplified_adjusted < 0.05`)) |>
	select(tissue,`P_ethnicity_simplified_adjusted < 0.05`, n) |>
	spread(`P_ethnicity_simplified_adjusted < 0.05`, n) |>
	mutate(proportion_significant = `TRUE` / (`FALSE` + `TRUE`)) |>
	arrange(desc(proportion_significant)) |>
	rowid_to_column("rank_de")



plot_ranks_tissue_barplot =
	rank_de_tissue |>
	ggplot(aes(`TRUE`, tissue |> fct_reorder(proportion_significant))) +
	geom_bar(aes(fill = tissue), stat = "identity") +
	scale_fill_manual(values = tissue_color) +
	guides(fill="none") +
	xlab("Number of significant genes") +
	theme_multipanel +
	theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())



plot_ranks_tissue =

		differential_composition_ethnicity_relative |>
			sccomp_test() |>
			filter(parameter |> str_detect("sexmale___")) |>
			with_groups(parameter, ~ .x |> summarise(median_FDR= median(c_FDR))) |>
			arrange(median_FDR) |>
			rowid_to_column("rank_composition") |>
			mutate(tissue = parameter |> str_remove("sexmale___")) |>
			select(tissue , rank_composition) |>


			# DE
			full_join( rank_de_tissue,	by = join_by(tissue)	) |>
			select(tissue, contains("rank")) |>
			pivot_longer(contains("rank")) |>
			mutate(name = name |> str_remove("rank")) |>

			# Parse cell types
			mutate(
				tissue =
					tissue |>
					str_replace_all("_", " ") |>
					str_replace("gland", "gld") |>
					str_replace("skeletal", "sk")
			) |>

			# Plot
			ggplot(aes(name, value, group=tissue,  label=tissue)) +
			geom_line(aes(color = tissue)) +
			geom_point(aes(color = tissue)) +
			geom_text(size = 2.5) +
			scale_color_manual(values = tissue_color) +
			scale_y_reverse() +
			guides(color = "none") +
			theme_multipanel


# UMAP
adjusted_ethnicity_fixed =
  differential_composition_ethnicity_relative |>
  sccomp_remove_unwanted_variation(~ ethnicity_simplified + group)

replicate_ethnicity_fixed =
  differential_composition_ethnicity_relative |>
  select(cell_type_harmonised, count_data) |>
  unnest(count_data) |>
  distinct() |>
  with_groups(sample_, ~ .x |> mutate(adjusted_proportion = count/sum(count)))

y =
  replicate_ethnicity_fixed |>
  tidybulk::as_SummarizedExperiment(
    sample_,
    cell_type_harmonised,
    c(adjusted_proportion)
  ) |>
  reduce_dimensions(method = "UMAP", scale=FALSE)

# UMAP by dataset and ethnicity
# This did not produce any good result
plot_PCA_ethnicity =
  adjusted_ethnicity_fixed |>
  tidybulk::as_SummarizedExperiment(
    sample_,
    cell_type_harmonised,
    c(adjusted_proportion)
  ) |>
  left_join(
    data_for_immune_proportion_relative |>
      distinct(sample_, ethnicity_simplified, tissue_harmonised, file_id)
  ) |>

  filter(
    tissue_harmonised=="blood"
  ) |>

  nest(data = -file_id) |>

  filter(map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow()) > 1) |>
  arrange(map_int(data, ncol) |> desc()) |>
  slice_head(n=50) |>

  # Filter
  mutate(data = map(
    data,
    ~ .x |>
      # Parse cell types
      mutate(cell_type_harmonised = cell_type_harmonised |>
               str_replace("megakaryocytes", "mega") |>
               str_remove("phage") |>
               str_replace("th1 th17", "th1/17") |>
               str_replace("mono", "mn") |>
               str_replace("tcm", "cm") |>
               str_replace("cd4 th", "Th") |>
               str_replace("memory", "mem") |>
               str_replace("naive", "nv") |>
               str_replace("terminal effector cd4 t", "cd4 eff") |>
               str_remove("cyte")
      )   |>
      inner_join(
        df_heatmap_ethnicity_relative_organ_cell_type |>
          filter(c_FDR < 0.05 & n > 5 ) |>
          select(tissue, cell_type) |>
          filter(!cell_type %in% c("thymo", "immune_unclassified")) |>
          distinct(tissue, cell_type),
        by =  c("tissue_harmonised" = "tissue", "cell_type_harmonised" = "cell_type")
      )
  )) |>
  slice_head(n = 7) |>
  dplyr::slice(c(1, 2, 5, 7)) |>
  mutate(data = map(data, reduce_dimensions, method = "PCA", scale = TRUE))  |>
  mutate(plot = map(
    data,
    ~ {

      # Do PCA
      PCA <- .x |> assay("adjusted_proportion") |> t() |>   prcomp(scale. = TRUE)

      # Extract PC axes for plotting
      PCAvalues <- data.frame(cell_type = colnames(.x), PCA$x, Ethnicity = .x |> pivot_sample() |> pull(ethnicity_simplified))

      # Extract loadings of the variables
      PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)

      # Plot
      gg =
        ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = Ethnicity)) +
        geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*5),
                                             yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
                     color = "black") +
        geom_point(size = 0.5) +
        annotate("text", x = (PCAloadings$PC1*5), y = (PCAloadings$PC2*5),
                 label = PCAloadings$Variables) +
        scale_color_manual(values = c("European" = "#E41A1C", "Chinese" = "#377EB8", "African" = "#4DAF4A", "Other" = "#984EA3")) +
        guides(color="none") +
        theme_multipanel

      ggExtra::ggMarginal(gg, groupColour = TRUE, groupFill = TRUE)
      },
    .progress = TRUE
  )) |>
  pull(plot) |>
  wrap_plots()


# # Which genes are shared across tissue programs
#
# readRDS("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue.rds") |>
# 	filter(map_lgl(se, ~ "P_ethnicity_simplified_adjusted" %in% colnames(.x))) |>
# 	mutate(se = map(se, ~ .x |> select(P_ethnicity_simplified_adjusted, .feature))) |>
# 	select(se, name) |>
# 	unnest(se) |>
# 	nest(tissues = -.feature) |>
# 	mutate(tissues_significant = map(
# 		tissues,
# 		~ .x |> filter(P_ethnicity_simplified_adjusted<0.05)
# 	)) |>
#
# 	mutate(
# 		n_tissue = map_int(tissues, nrow),
# 		n = map_int(tissues_significant, nrow)
# 	) |>
# 	filter(n_tissue >= 5) |>
# 	mutate(proportion = n/n_tissue) |>
# 	filter(proportion > 0.30) |>
# 	arrange(desc(proportion)) |>
# 	filter(!.feature %in% gene_chr$ID)

de_ethnicity_cell_type = readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_ethnicity_cell_type.rds")

plot_genes_most_altered =

  de_ethnicity_cell_type |>
  filter(map_lgl(data, ~ "P_ethnicity_simplified_adjusted" %in% colnames(.x))) |>
  mutate(data = map(data, ~ .x |> select(P_ethnicity_simplified_adjusted, .feature))) |>
  select(data, name) |>
  unnest(data) |>
  filter(name |> str_detect("immune_unclassified", negate = TRUE)) |>
  nest(cell_types = -.feature) |>

  mutate(cell_types_significant = map(
    cell_types,
    ~ .x |> filter(P_ethnicity_simplified_adjusted<0.05)
  )) |>

  mutate(
    n_cell_type = map_int(cell_types, nrow),
    n = map_int(cell_types_significant, nrow)
  ) |>
  filter(n_cell_type >= 5) |>
  mutate(proportion = n/n_cell_type) |>
  #filter(proportion > 0.30) |>
  arrange(desc(n)) |>
  mutate(cell_type_names = map_chr(cell_types_significant, ~ .x |> pull(name) |> str_c(collapse = ", ") |> str_remove_all("data_"))) |>

  # Plot
  mutate(class = case_when(
    row_number() <= 10 ~ "shared",
    row_number() == 11 ~ "CD4",
    row_number() %in% c(21, 23, 34) ~ "naive",
    row_number() == 41 ~ "Monocyte"
  )) |>

  mutate(label = if_else(!is.na(class), .feature, NA)) |>
  #mutate(n = as.character(n)) |>
  ggplot(aes(proportion, n, label = label)) +
  geom_point(aes(size = n_cell_type, color = class, alpha = !is.na(class))) +
  ggrepel::geom_text_repel(min.segment.length = 0, size=2.5, max.overlaps = 100) +
  scale_color_brewer(palette = "Set1", na.value="grey") +
  scale_size_continuous(range = c(0.2, 2)) +
  #scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1)) +
  ylab("Number of cell type with significant changes") +
  xlab("Proportion of significant cell types on all tested cell types") +
  xlim(c(NA, 0.4)) +
  theme_multipanel

# # Tissue DE
# job::job({
#
	# library(future)
	# library(furrr)
	# library("future.batchtools")
	#
	# slurm <-
	# 	`batchtools_slurm` |>
	# 	future::tweak( template = glue("{home}/third_party_sofware/slurm_batchtools.tmpl"),
	# 								 resources=list(
	# 								 	ncpus = 2,
	# 								 	memory = 4000,
	# 								 	walltime = 72800
	# 								 )
	# 	)
	# plan(slurm)
#
#
#
# 	se_adjust =
# 		tar_meta(store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue"), starts_with("data_")) |>
# 		filter(!is.na(data)) |>
#
# 		mutate(se = future_map(
# 			name,
# 			~ {
# 				se = tar_read_raw(.x, store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue"))
#
# 				if(SummarizedExperiment::ncol(se)==0) return(NULL)
# 				if(!"P_ethnicity_simplified_adjusted" %in% colnames(SummarizedExperiment::rowData(se))) return(NULL)
#
# 				# Get formula
# 				factor_unwanted =
# 					se |>
# 					attr("internals") %$%
# 					glmmseq_lme4 %>%
# 					`@` (formula) |>
# 					as.character() %>%
# 					.[[3]] |>
# 					str_split("\\+")  %>%
# 					.[[1]] |>
# 					str_subset("\\(|\\||sex", negate = TRUE) |>
# 					str_trim()
#
# 				# Make object lighter
# 				attr(se, "internals")$glmmseq_lme4 = NULL
#
# 				se |>
# 					adjust_abundance(
# 						.factor_of_interest = c(sex),
# 						.factor_unwanted = any_of(factor_unwanted)  ,
# 						.abundance = counts_scaled,
# 						method = "limma_remove_batch_effect"
# 					)  |>
# 					tidySummarizedExperiment::filter(P_ethnicity_simplified_adjusted<0.05) |>
# 					tidySummarizedExperiment::select(
# 						.feature, .sample, sample_,
# 						contains("counts"), contains("sexmale"),
# 						tissue_harmonised, cell_type_harmonised,
# 						sex, P_ethnicity_simplified_adjusted, P_ethnicity
# 					)
#
#
# 			},
# 			.progress=T
# 		)) |>
#
# 		# Get the effects
# 		filter(!map_lgl(se, is.null)) |>
#
# 		# Parse tissue
# 		mutate(tissue = name |> str_remove("data_")) |>
# 		mutate(tissue = case_when(
# 			tissue == "node" ~ "lymph node",
# 			tissue == "large" ~ "intestine large",
# 			tissue == "small" ~ "intestine small",
# 			tissue == "gland" ~ "adrenal gland",
# 			TRUE ~ tissue
# 		)) |>
#
# 		# join the non immune DE
# 		left_join(
# 			readRDS("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue_non_immune.rds") |>
# 				mutate(tissue = map_chr(se, ~ .x |> pull(tissue_harmonised) |> unique())) |>
# 				filter(!is.na(tissue)) |>
# 				filter(map_lgl(se, ~ "P_ethnicity_simplified_adjusted" %in% colnames(.x), .progress = TRUE)) |>
# 				mutate(de_genes = map(se, ~ .x |> filter( P_ethnicity_simplified_adjusted<0.05) |> pull(.feature))) |>
# 				select(tissue, de_genes),
# 			by = "tissue"
# 		) |>
#
# 		mutate(effects = map(
# 			se,
# 			~ .x |>
# 				# Significant
# 				filter(P_ethnicity_simplified_adjusted<0.05)  |>
#
# 				# Non sex genes
# 				filter(!.feature %in% gene_chr$ID) |>
#
# 				# Filter cell types for which we have both sexes
# 				as_tibble() |>
# 				nest(data = -c(.feature, cell_type_harmonised)) |>
# 				filter(map_int(data, ~ .x |> distinct(sex) |> nrow()) >1) |>
# 				filter(map_int(data, ~ .x |> filter(counts_scaled >1) |> nrow()) >0) |>
# 				unnest(data) |>
#
# 				nest(data = -c(.feature, contains("mode"), sexmale, P_ethnicity_simplified_adjusted)) |>
# 				mutate(cell_type_to_select = map(data, ~ .x |> pull(cell_type_harmonised) |> unique())) |>
# 				select(-data) |>
#
#
# 				pivot_longer(cols = contains("mode"), names_sep = "_cell_type_harmonised__", names_to = c("cell_type", "stat")) |>
# 				filter(map2_lgl(cell_type, cell_type_to_select, ~ .x %in% .y)) |>
# 				with_groups(c(.feature,P_ethnicity_simplified_adjusted, sexmale), ~ .x |> summarise(mode_sd = sd(value), n = n())) |>
# 				arrange(mode_sd) |>
#
# 				# Only filter confident
# 				#filter(n>5) |>
# 				filter(P_ethnicity_simplified_adjusted<0.001) |>
# 				filter(abs(sexmale)>1),
#
# 			.progress = TRUE
# 		)) |>
#
# 		# Separate by up or down regulation
# 		mutate(se = map(
# 			se,
# 			~ .x |>
# 				mutate(enriched_m_f = if_else(sexmale>0, "male", "female")) |>
# 				nest(se = -enriched_m_f),
# 			.progress = TRUE
# 		)) |>
# 		mutate(effects = map(
# 			effects,
# 			~ .x |>
# 				mutate(enriched_m_f = if_else(sexmale>0, "male", "female")) |>
# 				mutate(generic_vc_cell_specific = if_else(n>3, "generic", "specific")) |>
# 				nest(effects = -c(enriched_m_f, generic_vc_cell_specific))
# 		)) |>
# 		filter(map_int(effects, nrow) > 0) |>
# 		mutate(se = map2(
# 			se, effects,
# 			~ .x |>
# 				inner_join(.y)
# 		)) |>
# 		select(-effects) |>
# 		unnest(se) |>
#
# 		# Pathway analyses
# 		mutate(enriched_pathways = map2(
# 			se, effects,
# 			~ .x |>
# 				mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
# 																							keys = .feature,
# 																							keytype = "SYMBOL",
# 																							column = "ENTREZID",
# 																							multiVals = "first"
# 				)) |>
# 				test_gene_overrepresentation(
# 					.do_test = .feature %in% .y$.feature,
# 					species = "Homo sapiens", .entrez = entrez) ,
# 			.progress=TRUE
# 		))
#
# })
#
# se_adjust |> qs::qsave("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue_adjusted.qs")

# # Cell type DE
# job::job({
#
# 	library(future)
# 	library(furrr)
# 	library("future.batchtools")
#
# 	slurm <-
# 		`batchtools_slurm` |>
# 		future::tweak( template = glue("{home}/third_party_sofware/slurm_batchtools.tmpl"),
# 									 resources=list(
# 									 	ncpus = 2,
# 									 	memory = 40000,
# 									 	walltime = 72800
# 									 )
# 		)
# 	plan(slurm)
#
#
# 	se_adjust =
# 		tar_meta(store = glue("pseudobulk_0.2.3.4/_targets__ethnicity"), starts_with("data_")) |>
# 		filter(!is.na(data)) |>
#
# 		mutate(se = future_map(
# 			name,
# 			~ {
# 				se = tar_read_raw(.x, store = glue("pseudobulk_0.2.3.4/_targets__ethnicity"))
#
# 				if(SummarizedExperiment::ncol(se)==0) return(NULL)
# 				if(!"P_ethnicity_simplified_adjusted" %in% colnames(SummarizedExperiment::rowData(se))) return(NULL)
#
# 				# Get formula
# 				factor_unwanted =
# 					se |>
# 					attr("internals") %$%
# 					glmmseq_lme4 %>%
# 					`@` (formula) |>
# 					as.character() %>%
# 					.[[3]] |>
# 					str_split("\\+")  %>%
# 					.[[1]] |>
# 					str_subset("\\(|\\||sex", negate = TRUE) |>
# 					str_trim()
#
# 				# Make object lighter
# 				attr(se, "internals")$glmmseq_lme4 = NULL
#
# 				se |>
# 					adjust_abundance(
# 						.factor_of_interest = c(sex),
# 						.factor_unwanted = any_of(factor_unwanted)  ,
# 						.abundance = counts_scaled,
# 						method = "limma_remove_batch_effect"
# 					)  |>
# 					tidySummarizedExperiment::filter(P_ethnicity_simplified_adjusted<0.05) |>
# 					tidySummarizedExperiment::select(
# 						.feature, .sample,
# 						contains("counts"), contains("sexmale"),
# 						tissue_harmonised, cell_type_harmonised,
# 						sex, P_ethnicity_simplified_adjusted, P_ethnicity
# 					)
#
#
# 			},
# 			.progress=T
# 		)) |>
#
# 		# Get the effects
# 		filter(!map_lgl(se, is.null)) |>
#
# 		# If I have random effects
# 		filter(map_lgl(se, ~ .x[, 1] |> as_tibble() |> colnames() |> str_subset("mode") |> length() >0 )) |>
# 		mutate(effects = map(
# 			se,
# 			~ .x |>
# 				# Significant
# 				filter(P_ethnicity_simplified_adjusted<0.05)  |>
#
# 				# Non sex genes
# 				filter(!.feature %in% gene_chr$ID) |>
#
# 				# Filter cell types for which we have both sexes
# 				as_tibble() |>
# 				nest(data = -c(.feature, tissue_harmonised)) |>
# 				filter(map_int(data, ~ .x |> distinct(sex) |> nrow()) >1) |>
# 				filter(map_int(data, ~ .x |> filter(counts_scaled >1) |> nrow()) >0) |>
# 				unnest(data) |>
#
# 				nest(data = -c(.feature, contains("mode"), sexmale, P_ethnicity_simplified_adjusted)) |>
# 				mutate(cell_type_to_select = map(data, ~ .x |> pull(tissue_harmonised) |> unique())) |>
# 				select(-data) |>
#
#
# 				pivot_longer(cols = contains("mode"), names_sep = "_tissue_harmonised__", names_to = c("cell_type", "stat")) |>
# 				filter(map2_lgl(cell_type, cell_type_to_select, ~ .x %in% .y)) |>
# 				with_groups(c(.feature,P_ethnicity_simplified_adjusted, sexmale), ~ .x |> summarise(mode_sd = sd(value), n = n())) |>
# 				arrange(mode_sd) |>
#
# 				# Only filter confident
# 				#filter(n>5) |>
# 				filter(P_ethnicity_simplified_adjusted<0.001) |>
# 				filter(abs(sexmale)>1),
#
# 			.progress = TRUE
# 		)) |>
#
# 		# Separate by up or down regulation
# 		mutate(se = map(
# 			se,
# 			~ .x |>
# 				mutate(enriched_m_f = if_else(sexmale>0, "male", "female")) |>
# 				nest(se = -enriched_m_f),
# 			.progress = TRUE
# 		)) |>
# 		mutate(effects = map(
# 			effects,
# 			~ .x |>
# 				mutate(enriched_m_f = if_else(sexmale>0, "male", "female")) |>
# 				mutate(generic_vc_cell_specific = if_else(n>3, "generic", "specific")) |>
# 				nest(effects = -c(enriched_m_f, generic_vc_cell_specific))
# 		)) |>
# 		filter(map_int(effects, nrow) > 0) |>
# 		mutate(se = map2(
# 			se, effects,
# 			~ .x |>
# 				inner_join(.y)
# 		)) |>
# 		select(-effects) |>
# 		unnest(se) |>
#
# 		# Pathway analyses
# 		mutate(enriched_pathways = map2(
# 			se, effects,
# 			~ .x |>
# 				mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
# 																							keys = .feature,
# 																							keytype = "SYMBOL",
# 																							column = "ENTREZID",
# 																							multiVals = "first"
# 				)) |>
# 				test_gene_overrepresentation(
# 					.do_test = .feature %in% .y$.feature,
# 					species = "Homo sapiens", .entrez = entrez) ,
# 			.progress=TRUE
# 		))
#
# })
#
# se_adjust_cell_type |> qs::qsave("sccomp_on_HCA_0.2.3.5/de_ethnicity_cell_type_adjusted.qs")




# Select tissue pathways
xx =
	qs::qread("sccomp_on_HCA_0.2.3.5/de_ethnicity_tissue_adjusted.qs") |>
	select(name, enriched_m_f, enriched_pathways, generic_vc_cell_specific) |>
	mutate(enriched_pathways = map(
		enriched_pathways,
		~ .x |>
			filter(p.adjust<0.05) |>
			mutate(GeneRatio = map_dbl(GeneRatio, ~ parse(text = .x) |> eval())) |>
			arrange(desc(GeneRatio)) |>
			head(10) |>

			# Convert back to entrez
			mutate(symbol = map(
				entrez,
				~ AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
																keys = .x,
																keytype = "ENTREZID",
																column = "SYMBOL",
																multiVals = "first"
				) |> suppressMessages()
			))
	)) |>
	unnest(enriched_pathways)



# Compose plot with patchwork
plot =

	# Row 1
	((
		((
			plot_ethnicity_absolute_1D |
			  plot_spacer() |
			  plot_PCA_ethnicity
		) +
			plot_layout(width = c(0.5, 2, 1))
		) /

		# Row 2
		((
			wrap_heatmap(
				plot_heatmap_ethnicity_relative_organ_cell_type,
				padding = grid::unit(c(-50, -20, -0, -5), "points" ),
				clip = FALSE
			) |
				plot_ranks_cell_type |
				plot_ranks_cell_type_barplot |
				plot_ranks_tissue | # Add bar plot
				plot_ranks_tissue_barplot
		) + plot_layout(width = c(1.5, 0.9, 0.2, 1.1, 0.2))
	) /

		# Row 3
		((
			(		wrap_plots(venn_0_3) / wrap_plots(venn_0_5) ) |
			plot_genes_most_altered |
			plot_spacer() |
			plot_spacer()
		) + plot_layout(width = c(1, 0.8, 1, 1))
	)

	) + plot_layout( height = c(55, 60, 55))
	) +
	plot_layout(width = c(84, 97)) &
	theme(
		plot.margin = margin(0, 0, 0, 0, "pt"),
		legend.key.size = unit(0.2, 'cm'),
		legend.position = "bottom"
	)


ggsave(
	glue("{result_directory}/figure_ethnicity.pdf"),
	plot = plot,
	units = c("mm"),
	width = 183 ,
	height = 200 ,
	limitsize = FALSE
)




# Pathway analyses of consistent vs diverse genes across cell type
genes_consistent = xx|> head(100) |> pull(.feature)
genes_diverse = xx |> tail(100) |> pull(.feature)

x1 =
	se_adjust |>
	pull(se) %>%
	.[[3]] |>
	filter(sexmale>0) |>
	mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
																				keys = .feature,
																				keytype = "SYMBOL",
																				column = "ENTREZID",
																				multiVals = "first"
	)) |>
	test_gene_overrepresentation(
		.do_test = .feature %in% genes_consistent,
		species = "Homo sapiens", .entrez = entrez)

x2 = se_adjust |>
	pull(se) %>%
	.[[3]] |>
	filter(sexmale>0) |>
	mutate(entrez = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
																				keys = .feature,
																				keytype = "SYMBOL",
																				column = "ENTREZID",
																				multiVals = "first"
	)) |>
	test_gene_overrepresentation(
		.do_test = .feature %in% genes_diverse,
		species = "Homo sapiens", .entrez = entrez)

# Plot the distribution of a gene across the tissue results
se_adjust |>
	select(name, se) |>
	filter(!map_lgl(se, is.null)) |>
	mutate(se = map(se,
									~ .x |>
										filter(.feature== "LONRF1") |>
										as_tibble()
	)) |>
	unnest(se) |>
	ggplot(aes(sex, counts_scaled_adjusted + 1)) +
	geom_boxplot(fill ="white", outlier.shape = NA) +
	geom_jitter(aes(color = tissue_harmonised), size = 0.2, height =0) +
	facet_grid(tissue_harmonised ~ cell_type_harmonised , scales = "free_y") +
	scale_y_log10() +
	theme_bw()

# Pathway analyses for tissues unique
random_effects_fold_changes

# Export single cell
bone_pseudobulk_samples =
	tar_read_raw("data_bone", store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue")) |>
	distinct(sample_) |>
	tidyr::extract(col = sample_, into = "sample_", regex = "([a-zA-Z0-9]+)_*") |>
	pull(sample_)


library(DelayedArray)
library(HDF5Array)


get_metadata() |>
	#	filter(cell_type_harmonised == "cd4 th17") |>
	filter(tissue_harmonised %in% c("bone")) |>
	filter(sample_ %in% bone_pseudobulk_samples) |>
	add_count(sample_, name = "count_sample") |>
	filter(count_sample > 200) |>
	as_tibble()  |>
	#with_groups(sample_, ~ .x |> sample_n(min(1000, n() ))) |>
	get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts") |>
	saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.5/bone")

# kidney
kidney_pseudobulk_samples =
	tar_read_raw("data_kidney", store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue")) |>
	distinct(sample_) |>
	tidyr::extract(col = sample_, into = "sample_", regex = "([a-zA-Z0-9]+)_*") |>
	pull(sample_)

get_metadata() |>
	#	filter(cell_type_harmonised == "cd4 th17") |>
	filter(tissue_harmonised %in% c("kidney")) |>
	filter(sample_ %in% kidney_pseudobulk_samples) |>
	add_count(sample_, name = "count_sample") |>
	filter(count_sample > 200) |>
	as_tibble()  |>
	with_groups(sample_, ~ .x |> sample_n(min(2000, n() ))) |>
	get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts") |>
	saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.5/kidney")


# adipose
adipose_pseudobulk_samples =
	tar_read_raw("data_adipose", store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue")) |>
	distinct(sample_) |>
	tidyr::extract(col = sample_, into = "sample_", regex = "([a-zA-Z0-9]+)_*") |>
	pull(sample_)


get_metadata() |>
	#	filter(cell_type_harmonised == "cd4 th17") |>
	filter(tissue_harmonised %in% c("adipose")) |>
	filter(sample_ %in% adipose_pseudobulk_samples) |>
	add_count(sample_, name = "count_sample") |>
	filter(count_sample > 200) |>
	as_tibble()  |>
	#with_groups(sample_, ~ .x |> sample_n(min(1000, n() ))) |>
	get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts") |>
	saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.5/adipose")





do.call(
	cbind,
	tar_meta(store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue"), starts_with("se_filtered")) |>
		filter(!is.na(data)) |>
		mutate(se = map(
			name,
			~ tar_read_raw(.x, store = glue("pseudobulk_0.2.3.4/_targets__ethnicity_tissue")),
			.progress=T
		)) |>
		pull(se)
)





se |>
	filter(.feature== "NCOA7") |>
	ggplot(aes(sex, counts_scaled + 1)) +
	geom_boxplot(fill ="white", outlier.shape = NA) +
	geom_jitter(aes(color = tissue_harmonised), size = 0.2, height =0) +
	facet_wrap(~ tissue_harmonised, scales = "free_y") +
	scale_y_log10() +
	theme_bw()


# Clean environment
rm(differential_composition_ethnicity_relative , differential_composition_ethnicity_absolute )
gc()


# Differential Expression
library(ggExtra)
library(magrittr)

cell_types = readRDS(glue("{result_directory}/input_common.rds")) |> distinct(tissue_harmonised) |> pull(tissue_harmonised) |> str_replace(" ", "_")
de_table =
	dir(glue("{result_directory}/DGEA"), full.names = TRUE) |>
	enframe(value="file") |>
	mutate(de = map(file, read_csv)) |>
	mutate(tissue_harmonised = map(
		file,
		~ {
			for(ct in cell_types){
				if(.x |> tolower()  |> str_subset(ct |> tolower()) |> length() |> equals(1) )
					return(ct)
			}
		}
	)) |>
	unnest(tissue_harmonised) |>
	mutate(tissue_harmonised = tissue_harmonised |> toupper()) |>
	tidyr::extract(  file, "cell_type", ".+/.+/(.+)\\.csv") |>
	mutate(cell_type = cell_type |> str_remove(tissue_harmonised) |> str_remove("^_") ) |>

	# Add cell typoe lineage
	left_join(
		tribble(
			~cell_type, ~lineage,
			"B_MEMORY",	"b",
			"B_NAIVE",	"b",
			"CD14_MONO"	,"myeloid",
			"CD16_MONO",	"myeloid",
			"CD4_FH",	"cd4",
			"CD4_NAIVE",	"cd4",
			"CD4_TH1"	,"cd4",
			"CD4_TH17",	"cd4",
			"CD4_TH1_TH17",	"cd4",
			"CD4_TH2"	,"cd4",
			"CD8_TCM"	,"cd8",
			"CD8_TEM"	,"cd8",
			"CDC"	,"myeloid",
			"GRANULOCYTE"	,"granulocyte",
			"ILC"	,"ilc",
			"IMMUNE_UNCLASSIFIED"	,"nn",
			"MACROPHAGE",	"myeloid",
			"MAIT","mait",
			"MAST",	"mast",
			"PDC",	"pdc",
			"PLASMA",	"b",
			"TERMINAL_EFFECTOR_CD4_T"	,"cd4",
			"TGD"	,"tgd",
			"THYMOCYTE"	,"thymocyte",
			"TREG"	,"cd4"
		)
	) |>
	filter(cell_type != "IMMUNE_UNCLASSIFIED") |>
	unnest(de)





# de_table |>
#   distinct(cell_type, tissue_harmonised, symbol, lineage) |>
#   filter(!symbol %in% gene_chr$ID) |>
#   mutate(exists = 1) |>
#   complete(nesting(cell_type,tissue_harmonised, lineage), symbol, fill = list(exists = 0)) |>
#   mutate(cell_type_tissue = glue("{cell_type} {tissue_harmonised}") ) |>
#
#   heatmap(
#     cell_type_tissue,
#     symbol,
#     exists
#   ) |>
#   annotation_tile(cell_type) |>
#   annotation_tile(lineage) |>
#   annotation_tile(tissue_harmonised)

de_table |>
	with_groups(tissue_harmonised, ~ .x |> dplyr::count(cell_type) |> summarise(mean_de = median(n)))
dplyr::count(cell_type, tissue_harmonised)

de_table |>
	dplyr::count(cell_type, tissue_harmonised) |>
	arrange(desc(n)) |>
	pull(n) |> hist()

de_table |>
	dplyr::count(symbol) |>
	arrange(desc(n)) |>
	pull(n) |> hist()

de_table_tissue =
	de_table |>
	nest(data = -c(tissue_harmonised)) |>
	mutate(n_celltype_in_tissue = map_int(data, ~ .x |> distinct(cell_type) |> nrow())) |>
	mutate(shared_genes = map2(
		data, n_celltype_in_tissue,
		~ .x |>
			dplyr::count(symbol) |>
			filter(n/.y>=0.75) |>
			pull(symbol)
	))


de_table_celltype =
	de_table |>
	nest(data = -c(cell_type)) |>
	mutate(n_tissue_in_celltype = map_int(data, ~ .x |> distinct(tissue_harmonised) |> nrow())) |>
	mutate(shared_genes = map2(
		data, n_tissue_in_celltype,
		~ .x |>
			dplyr::count(symbol) |>
			filter(n/.y>=0.75) |>
			pull(symbol)
	))

# Select the specific and unique
de_table_celltype |>
	mutate(shared_genes_specific = map(
		shared_genes,
		~ .x |> setdiff(
			de_table_tissue |>
				filter(n_celltype_in_tissue>2) |>
				dplyr::select(shared_genes) |>
				unnest(shared_genes) |>
				pull(1)
		)
	)) |>
	arrange(map_int(shared_genes_specific, length) |> desc()) |>
	filter(n_tissue_in_celltype > 2)


de_table_tissue |>
	mutate(shared_genes_specific = map(
		shared_genes,
		~ .x |> setdiff(
			de_table_celltype |>
				filter(n_tissue_in_celltype>2) |>
				dplyr::select(shared_genes) |>
				unnest(shared_genes) |>
				pull(1)
		)
	)) |>
	arrange(map_int(shared_genes_specific, length) |> desc()) |>
	filter(n_celltype_in_tissue > 2)


ilc =
	get_metadata() |>
	filter(cell_type_harmonised == "ilc") |>
	filter(tissue_harmonised %in% c("lymph node", "heart", "blood", "kidney", "liver", "lung")) |>
	add_count(sample_, name = "count_sample") |>
	filter(count_sample > 200) |>
	as_tibble()  |>
	with_groups(sample_, ~ .x |> sample_n(min(1000, n() ))) |>
	get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts")

library(DelayedArray)
library(HDF5Array)
ilc |> saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.5/ilc")


th17 =
	get_metadata() |>
	filter(cell_type_harmonised == "cd4 th17") |>
	filter(tissue_harmonised %in% c("lymph node", "heart", "blood", "kidney", "liver", "lung")) |>
	add_count(sample_, name = "count_sample") |>
	filter(count_sample > 200) |>
	as_tibble()  |>
	with_groups(sample_, ~ .x |> sample_n(min(1000, n() ))) |>
	get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts")

th17 |> saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.5/th17")


plot_count_de =
	de_table |>
	unnest(de) |>
	filter(adj.P.Val<0.05) |>
	dplyr::count(cell_type, tissue_harmonised) |>

	with_groups(cell_type, ~ .x |> mutate(median_n = median(n))) |>

	ggplot(aes(fct_reorder(cell_type, -median_n), n, color = tissue_harmonised)) +
	geom_point(aes()) +
	scale_y_log10() +
	theme_multipanel +
	theme(axis.text.x = element_text(angle=90))


ggMarginal(plot_count_de, groupColour = TRUE, groupFill = TRUE)
