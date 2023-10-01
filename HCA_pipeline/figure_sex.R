library(tidyverse)
library(forcats)
#library(CuratedAtlasQueryR)
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
result_directory = glue("{home}/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.7_double_interaction_sex_age")
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
# Sex analyses for immune cellularity
#------------------------------#

# Read results
differential_composition_sex_absolute_file = glue("{result_directory}/sex_absolute_FALSE.rds")
proportions_sex_absolute_file = glue("{result_directory}/sex_absolute_FALSE_proportion_adjusted.rds")
differential_composition_sex_absolute = readRDS(differential_composition_sex_absolute_file)
proportions_sex_absolute_adjusted = readRDS(proportions_sex_absolute_file)

# # save csv for SUPPLEMENTARY
# differential_composition_sex_absolute |>
# 	test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
# 	select(-count_data) |>
# 	write_csv("sccomp_on_HCA_0.2.3.4/SUPPLEMENTARY_sex_cellularity_estimates.csv")


# Get parameter draws for relative abundance 
# to manually plot the uncertainty
draws_abundance_sex = 
  differential_composition_sex_absolute   |> 
  sccomp:::get_abundance_contrast_draws(contrasts = c(sexmale = "sexmale")) |> 
  filter(is_immune=="TRUE")

# Get parameter draws for variability 
# to manually plot the uncertainty
draws_variability_sex = 
  differential_composition_sex_absolute   |> 
  sccomp:::get_variability_contrast_draws( contrasts = c(sexmale = "sexmale")) |> 
  filter(is_immune=="TRUE") |> 
  filter(parameter=="sexmale")

# Plot effect of composition and variability
# For immune cellularity (proportion of immune cells)
# Wtih uncertainty
plot_sex_absolute_1D =
  tibble(
    Variability = draws_variability_sex |> pull(.value),
    Abundance = draws_abundance_sex |> pull(.value)
  ) |>
  tidybulk::as_matrix() |>
  bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
  coord_flip() +
  xlab("Effect male immune cellularity") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 0.5))

draws_abundance_sex = 
  differential_composition_sex_absolute   |> 
  sccomp:::get_abundance_contrast_draws(contrasts = c(sexmale_ = "`age_days:sexmale`")) |> 
  filter(is_immune=="TRUE") |> 
  filter(parameter=="sexmale_")

# Get parameter draws for variability 
# to manually plot the uncertainty
draws_variability_sex = 
  differential_composition_sex_absolute   |> 
  sccomp:::get_variability_contrast_draws( contrasts = c(sexmale_ = "`age_days:sexmale`")) |> 
  filter(is_immune=="TRUE") |> 
  filter(parameter=="sexmale_")

# Plot effect of composition and variability
# For immune cellularity (proportion of immune cells)
# Wtih uncertainty
plot_sex_absolute_1D_interaction =
  tibble(
    Variability = draws_variability_sex |> pull(.value),
    Abundance = draws_abundance_sex |> pull(.value)
  ) |>
  tidybulk::as_matrix() |>
  bayesplot::mcmc_intervals(point_size = 1, inner_size = 0.5, outer_size = 0.25) +
  coord_flip() +
  xlab("Effect male immune cellularity") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 0.5))


# Plot for presentation B memory across tissues
# data_adjusted_absolutesex_interation =
# 	differential_composition_sex_absolute |>
# 	remove_unwanted_variation(~ age_days*sex , ~ 1)  |>
# 	inner_join(data_for_immune_proportion |> tidybulk::pivot_sample(sample_) )
# 
# data_adjusted_absolutesex_interation |> saveRDS(glue("{result_directory}/data_adjusted_absolutesex_interation.rds"))

data_adjusted_absolutesex_interation = readRDS(glue("{result_directory}/data_adjusted_absolutesex_interation.rds"))

#res_relative_proportions_sex_tissue |> saveRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/res_relative_proportions_sex_tissue.rds")

# Get trend line to be used in the scatter plot of global compositional changes, plot_sex_relative
prediction_df = 
  expand_grid(sex = c("male", "female"), age_days = seq(-3, 3, by = 0.1)) |> 
  mutate(sample_ = 1:n() |> as.character())

line_sex_relative_mean =
  
  differential_composition_sex_absolute |>
  sccomp_predict(
    ~ age_days*sex, 
    new_data = prediction_df, 
    number_of_draws = 1
  ) |> 
  left_join(prediction_df) |> 
  mutate(x_corrected = (age_days * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0)) |>
  filter(is_immune == "TRUE") 

data_adjusted_absolutesex_interation |>
  filter(is_immune == "TRUE") |> 
  ggplot(aes(age_days_original, adjusted_proportion)) +
  geom_point(
    aes(fill = sex),
    shape = 21,
    stroke = 0,
    size = 0.6
  ) +
  geom_line(
    aes(x_corrected, proportion_mean, color = sex, group=sex),
    data = line_sex_relative_mean
  ) +
  geom_smooth(
    aes(color = sex), 
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE, size=1) +
  # scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  scale_fill_manual(values = c(female = "red", male = "blue")) +
  scale_color_manual(values = c(female = "red", male = "blue")) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none", color = "none") +
  theme_multipanel


# Create dataset to create the mannequin heatmap of the 
# Tissues with differential immune cellularity
sex_absolute_organ_tissue =
  differential_composition_sex_absolute |> 
  test_contrasts(
    contrasts =
      differential_composition_sex_absolute |>
      filter(parameter |> str_detect("^sexmale___")) |>
      distinct(parameter) |>
      mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", ".+___(.+)") |>
      filter(contrast |> str_detect("_female", negate = TRUE)) |> 
      deframe( ),
    test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
  ) |>
  filter(is_immune == "TRUE") |> 
  filter(parameter != "nose") |> 
  arrange(desc(abs(c_effect)))


# Create dataset to create the mannequin heatmap of the 
# Tissues with differential immune cellularity
sex_absolute_organ_tissue_interaction =
  differential_composition_sex_absolute |> 
  test_contrasts(
    contrasts =
      differential_composition_sex_absolute |>
      filter(parameter |> str_detect("age_days:sexmale___")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days:sexmale + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", ".+___(.+)") |>
      filter(contrast |> str_detect("_female", negate = TRUE)) |> 
      deframe( ),
    test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
  ) |>
  filter(is_immune == "TRUE") |> 
  filter(parameter != "nose") |> 
  arrange(desc(abs(c_effect)))


# # Plot for presentation B memory across tissues
# data_adjusted_absolutesex_interation_tissue =
# 	differential_composition_sex_absolute |>
# 	remove_unwanted_variation(~ age_days*sex + (age_days*sex | tissue_harmonised) , ~ 1)  |>
# 	inner_join(data_for_immune_proportion |> tidybulk::pivot_sample(sample_) )
# 
# data_adjusted_absolutesex_interation_tissue |> saveRDS(glue("{result_directory}/data_adjusted_absolutesex_interation_tissue.rds"))

data_adjusted_absolutesex_interation_tissue = readRDS(glue("{result_directory}/data_adjusted_absolutesex_interation_tissue.rds"))

prediction_df = 
  expand_grid(
    sex = c("male", "female"), 
    age_days = seq(-3, 3, by = 0.1), 
    tissue_harmonised = c("eye", "heart", "kidney", "lung")
  ) |> 
  mutate(sample_ = 1:n() |> as.character()) |> 
  expand_grid()

line_sex_relative_mean_per_tisue =
  
  differential_composition_sex_absolute |>
  sccomp_predict(
    ~ age_days*sex + (age_days*sex | tissue_harmonised), 
    new_data = prediction_df, 
    number_of_draws = 1
  ) |> 
  left_join(prediction_df) |> 
  mutate(x_corrected = (age_days * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0)) |>
  filter(is_immune == "TRUE") 

data_adjusted_absolutesex_interation_tissue |> 
  filter(tissue_harmonised %in% c("heart", "kidney", "lung")) |> 
  filter(is_immune == "TRUE") |> 
  ggplot(aes(age_days_original, adjusted_proportion)) +
  geom_point(
    aes(fill = sex),
    shape = 21,
    stroke = 0,
    size = 0.6
  ) +
  geom_smooth(
    aes(color=sex),
    method = "glm",
    method.args = list(family = "binomial"),
    se = FALSE, size=0.2) +
  # geom_line(
  # 	aes(x_corrected, proportion_mean, color = sex, group=sex),
  # 	data = line_sex_relative_mean_per_tisue
  # ) +
  facet_wrap(~tissue_harmonised, scales = "free") +
  
  # scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  scale_fill_manual(values = c(female = "red", male = "blue")) +
  scale_color_manual(values = c(female = "red", male = "blue")) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none", color = "none") +
  theme_multipanel

#res_relative_proportions_sex_tissue |> saveRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/res_relative_proportions_sex_tissue.rds")

# Get trend line to be used in the scatter plot of global compositional changes, plot_sex_relative
prediction_df = 
  expand_grid(sex = c("male", "female"), age_days = seq(-3, 3, by = 0.1)) |> 
  mutate(sample_ = 1:n() |> as.character()) 

line_sex_relative_mean_b_memory_per_tisue =
  
  differential_composition_sex_absolute |>
  sccomp_predict(
    ~ age_days*sex, 
    new_data = prediction_df, 
    number_of_draws = 1
  ) |> 
  left_join(prediction_df) |> 
  mutate(x_corrected = (age_days * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0)) |>
  filter(is_immune == "TRUE") 

data_adjusted_absolutesex_interation |>
  filter(is_immune == "TRUE") |> 
  ggplot(aes(age_days_original, adjusted_proportion)) +
  geom_point(
    aes(fill = sex),
    shape = 21,
    stroke = 0,
    size = 0.6
  ) +
  geom_line(
    aes(x_corrected, proportion_mean, color = sex, group=sex),
    data = line_sex_relative_mean_b_memory_per_tisue
  ) +
  # geom_smooth(
  #   method = "glm",
  #   method.args = list(family = "binomial"),
  #   se = FALSE, size=1) +
  # scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  scale_fill_manual(values = c(female = "red", male = "blue")) +
  scale_color_manual(values = c(female = "red", male = "blue")) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none", color = "none") +
  theme_multipanel



# # save csv for SUPPLEMENTARY
# differential_composition_sex_absolute |> 
# 	test_contrasts(
# 		contrasts =
# 			differential_composition_sex_absolute |>
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
# 	write_csv("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/SUPPLEMENTARY_sex_cellularity_tissue_estimates_contrasts.csv")


# Draw the color palette for the mannequin heatmap of the 
# Tissues with differential immune cellularity
colors_palette_for_organ_abundance =
  sex_absolute_organ_tissue |>
  select(parameter, c_effect) |> 
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
  pull(color) |>
  scales::show_col(	cex_label = 0.5	)


# 
# # Draw the boxplot for significant changes of the 
# # Tissues with differential immune cellularity
# plot_sex_absolute_organ_boxoplot_adjusted =
#   proportions_sex_absolute_adjusted |>
#   left_join(
#     data_for_immune_proportion |>
#       distinct(sample_, tissue_harmonised, ethnicity, sex,tissue, file_id)
#   ) |>
#   inner_join(
#     sex_absolute_organ_tissue |> 
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
# Sex analyses for immune composition
#------------------------------#

differential_composition_sex_relative_file = glue("{result_directory}/sex_relative_FALSE.rds")
proportions_sex_relative_file = glue("{result_directory}/sex_relative_FALSE_proportion_adjusted.rds")
differential_composition_sex_relative = readRDS(differential_composition_sex_relative_file)  
proportions_sex_relative_adjusted = readRDS(proportions_sex_relative_file)
gene_chr = read_csv("symbol_chr.csv")

# # save csv for SUPPLEMENTARY
# differential_composition_sex_relative |>
# 	test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
# 	select(-count_data) |>
# 	write_csv("sccomp_on_HCA_0.2.3.4/SUPPLEMENTARY_sex_composition_estimates.csv")


# # Volcano plot of cell type difference overall between sexes
# volcano_relative_sex = 
#   differential_composition_sex_relative |>
#   
#   test_contrasts(test_composition_above_logit_fold_change = 0.1) |> 
#   filter(cell_type_harmonised != "non_immune") |> 
#   filter(cell_type_harmonised != "immune_unclassified") |> 
#   
#   filter(parameter == "sexmale") |> 
#   mutate(naive_experienced = case_when(
#     cell_type_harmonised |> str_detect("naive") ~ "Antigen naive lymphcites",
#     cell_type_harmonised |> str_detect("cd4|cd8|memory") ~ "Antigen experienced lymphcites"
#   )) |> 
#   mutate(significant = c_FDR<0.05) |> 
#   mutate(cell_type_harmonised = case_when(c_FDR<0.05~cell_type_harmonised)) |> 
#   ggplot(aes(c_effect, c_FDR)) + 
#   
#   geom_vline(xintercept = -0.1, color="grey", linetype = "dashed", size=0.2) +
#   geom_vline(xintercept = 0.1, color="grey", linetype = "dashed", size=0.2) +
#   geom_hline(yintercept = 0.05, color="grey", linetype = "dashed", size=0.2) +
#   
#   geom_point(aes(color=naive_experienced, size=significant)) +
#   ggrepel::geom_text_repel(aes(label = cell_type_harmonised), size = 1.5 ) +
#   scale_y_continuous(trans = tidybulk::log10_reverse_trans()) + 
#   scale_x_continuous(trans="S_sqrt") +
#   scale_color_brewer(palette="Set1", na.value = "grey50") +
#   scale_size_discrete(range = c(0, 0.5)) +
#   guides(size="none", color="none") +
# 	ylab("False-discovery rate") +
# 	xlab("Effect from baseline (female)") + 
#   theme_multipanel

tissue_baseline_plot_sex_relative = 
  differential_composition_sex_relative |> 
  sccomp_replicate(formula_composition = ~ 1 + age_days  + (1 + age_days  | tissue_harmonised)) |>
  left_join(data_for_immune_proportion_relative |> distinct(sample_, tissue_harmonised)) |> 
  with_groups(c(cell_type_harmonised, tissue_harmonised), ~ .x |> summarise(
    lower = generated_proportions |> quantile(0.05),
    mean = generated_proportions |> mean(),
    upper = generated_proportions |> quantile(0.95),
  )) |> 
  nest(data = -cell_type_harmonised) |> 
  filter(map_int(data, ~ .x |> filter(mean > 0.001) |> nrow()) > 3) |> 
  unnest(data)

variability_abundance_plot = 
  differential_composition_sex_relative |> 
  test_contrasts(test_composition_above_logit_fold_change = 0.2) |> 
  filter(parameter=="sexmale") |> 
  filter(! cell_type_harmonised %in% c("dnt", "immune_unclassified")) |> 
  mutate(
    c_lower = case_when(c_FDR<0.05 ~ c_lower),
    c_upper = case_when(c_FDR<0.05 ~ c_upper),
    v_lower = case_when(v_FDR<0.05 ~ v_lower),
    v_upper = case_when(v_FDR<0.05 ~ v_upper),
  ) |> 
  
  # Filter cell types that are more than 3 tissues, 
  # otherwise a body level inference does not make much sense
  inner_join(    tissue_baseline_plot_sex_relative |> distinct(cell_type_harmonised) ) |> 
  
  # Shorten names
  mutate(cell_type_harmonised = cell_type_harmonised |> 
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
  
  ggplot(aes(c_effect, v_effect, label=cell_type_harmonised)) + 
  
  geom_vline(xintercept = c(-0.2, 0.2), color="grey", linetype="dashed") + 
  geom_hline(yintercept = c(-0.2, 0.2), color="grey", linetype="dashed") + 
  
  geom_linerange(aes(ymin = v_lower, ymax = v_upper), color="#D12424", alpha = 0.8) +
  geom_linerange(aes(xmin = c_lower, xmax = c_upper), color = "#4297C6", alpha=0.8) +
  
  geom_point() + 
  ggrepel::geom_text_repel(size = 2.5) +
  coord_cartesian(ylim=c(NA, 1.2)) +
  theme_multipanel


# Scatter plot of the significant cell types which compositon change through ageing
# The proportions are adjusted to exclude other effects including
# Sex, ethnicity, random effects (e.g. datasets) and technology

# outliers_df =
#   differential_composition_sex_relative |>
#   select(cell_type_harmonised, count_data) |>
#   unnest(count_data) |>
#   distinct() |>
#   filter(outlier) |>
#   select(cell_type_harmonised, sample_)


significant_cell_types_plot_sex_relative = 
  differential_composition_sex_relative |>
  test_contrasts(test_composition_above_logit_fold_change = 0.2) |>
  filter(parameter == "sexmale") |>
  filter(c_FDR<0.05 | v_FDR < 0.05) |>
  mutate(is_compositionally_different = c_FDR<0.05) |> 
  distinct(cell_type_harmonised, c_effect, is_compositionally_different) |> 
  
  # Filter cell types that are more than 3 tissues, 
  # otherwise a body level inference does not make much sense
  inner_join(    tissue_baseline_plot_sex_relative |> distinct(cell_type_harmonised) ) 

confidence_interval_plot_sex_relative = 
  differential_composition_sex_relative |> 
  sccomp_replicate(formula_composition = ~ sex) |>
  left_join(data_for_immune_proportion_relative |> distinct(sample_, sex)) |> 
  with_groups(c(sex, cell_type_harmonised), ~ .x |> summarise(
    lower = generated_proportions |> quantile(0.05),
    mean = generated_proportions |> mean(),
    upper = generated_proportions |> quantile(0.95),
  )) |> 
  inner_join(significant_cell_types_plot_sex_relative) 

plot_sex_relative =
  ggplot() +
  geom_jitter(
    aes(sex, adjusted_proportion, fill = tissue_harmonised),
    shape = 21,
    stroke = 0,
    size = 0.4,
    data = 
      proportions_sex_relative_adjusted |>
      
      # Attach sex
      inner_join(
        data_for_immune_proportion_relative |>
          tidybulk::pivot_sample(sample_)
      ) |>
      
      # Filter cell types that are more than 3 tissues, 
      # otherwise a body level inference does not make much sense
      inner_join(    tissue_baseline_plot_sex_relative |> filter(mean>0.005) |>  distinct(tissue_harmonised, cell_type_harmonised)  ) |> 
      
      # Drop outliers
      # anti_join( outliers_df ) |>
      
      # Filter
      filter(cell_type_harmonised != "immune_unclassified")  |>
      
      # Filter granulocytes because different in blood or solid
      #filter(cell_type_harmonised != "dnt") |>
      
      inner_join(significant_cell_types_plot_sex_relative ) |>
      
      filter(development_stage != "unknown") |>
      filter(!cell_type_harmonised %in% c("immune_unclassified", "dnt")) |>
      #filter(tissue_harmonised != "blood") |>
      # Fix samples with multiple assays
      unite("sample_", c(sample_ , assay), remove = FALSE) |>
      
      # Fix groups
      unite("group", c(tissue_harmonised , file_id), remove = FALSE) 
    
  ) +
  geom_point(aes(sex, mean),
             data = confidence_interval_plot_sex_relative,
             color = "red"
  ) +
  geom_line(aes(sex, mean, group = cell_type_harmonised),
            data = confidence_interval_plot_sex_relative,
            color = "red"
  ) +
  geom_linerange(aes(x = sex,ymin = lower, ymax = upper), color = "grey29", size = 0.5, data = confidence_interval_plot_sex_relative) +
  facet_wrap( ~ fct_reorder(cell_type_harmonised, !is_compositionally_different), nrow=2, scale="free_y") +
  scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_fill_manual(values = tissue_color) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Years") +
  ylab("Cell proportion (sqrt scale)") +
  guides(fill = "none", color = "none") +
  theme_multipanel

# # Plot for presentation B memory across tissues
# res_relative_proportions_sex_tissue = 
#   differential_composition_sex_relative |>
#   remove_unwanted_variation(~ 1 + sex + tissue_harmonised + ( 1 + sex | tissue_harmonised ), ~ sex)  |> 
#   inner_join(data_for_immune_proportion_relative |> tidybulk::pivot_sample(sample_) )
# 
# res_relative_proportions_sex_tissue |> saveRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/res_relative_proportions_sex_tissue.rds")

res_relative_proportions_sex_tissue = readRDS(glue("{result_directory}/res_relative_proportions_sex_tissue.rds"))


# # Get trend line to be used in the scatter plot of global compositional changes, plot_sex_relative
# line_sex_relative_mean_b_memory_per_tisue =
#   
#   differential_composition_sex_relative |>
#   filter(cell_type_harmonised=="b memory") |>
#   nest(data = -cell_type_harmonised) |>
#   
#   # Add tissue
#   mutate(
#     tissue_harmonised = 
#       list(
#         differential_composition_sex_relative |> 
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
#       (..2 * ..1 |> filter(parameter == glue("{..3}___sex_days")) |> pull(c_effect))   # GROUP-LEVEL SLOPE
#     
#   })) |>
#   dplyr::select(-data) |>
#   unnest(c(x, y)) |>
#   with_groups(x, ~ .x |> mutate(proportion = softmax(y))) |>
#   mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
#   filter(x_corrected |> between(30.0 , 30295.0))
# 
# res_relative_proportions_sex_tissue |>
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
#   # 		line_sex_relative_mean_b_memory_per_tisue 
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

# res_relative_proportions_sex_tissue |> 




# # Fold change for cell type by organ
# # These statistics are used in the paper result section
# differential_composition_sex_relative |> 
#   
#   print_estimate_plus_minus_relative(
#     contrasts_baseline = 
#       differential_composition_sex_relative |>
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
#       differential_composition_sex_relative |>
#       filter(parameter |> str_detect("___sex_days")) |>
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
#       differential_composition_sex_relative |>
#       filter(parameter |> str_detect("___sex_days")) |>
#       distinct(parameter) |>
#       mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
#       tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE)  |> 
#       filter(parameter |> str_detect("adipose", negate = TRUE)) |> 
#       mutate(parameter = glue("__{parameter}")) |> 
#       
#       deframe()
#   )

# Prepare the dataset for drawing the heatmap of changes by tissues and celltypes, Figure 3
df_heatmap_sex_relative_organ_cell_type =
  
  differential_composition_sex_relative |>
  
  # Find stats of random effect with groups
  test_contrasts(
    contrasts =
      differential_composition_sex_relative |>
      filter(parameter |> str_detect("sexmale___")) |>
      distinct(parameter) |>
      mutate(contrast = glue("sexmale + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", ".+___(.+)") |>
      deframe( ),
    #test_composition_above_logit_fold_change = FDR_threshold_1_percent_change_at_20_percent_baseline
  )  |>
  
  filter(cell_type_harmonised != "immune_unclassified") |>
  add_count(cell_type_harmonised) |>
  arrange(parameter, desc(n)) |>
  
  rename(tissue = parameter) |>
  rename(cell_type = cell_type_harmonised) |>
  
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
  
  # Filter for visualisation
  filter(!cell_type %in% c("non_immune", "immune_unclassified")) |>
  
  # Tissue diversity
  with_groups(tissue, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(tissue_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>
  
  # First rank
  with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>
  
  # # Cap
  # mutate(c_effect = c_effect |> pmax(-5) |> pmin(5)) |>
  mutate(Difference = c_effect) |>
  
  rename(`Mean diff` = cell_type_mean_change) |>
  mutate(`Mean diff tissue` = -tissue_mean_change) |>
  mutate(cell_type = cell_type |> str_replace("macrophage", "macro")) |>
  mutate(tissue = tissue |> str_replace_all("_", " ")) |>
  
  # Color
  left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
  left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  )  |>
  
  # Counts
  left_join(
    data_for_immune_proportion_relative |>
      count(tissue_harmonised, name = "count_tissue") |>
      rename(tissue = tissue_harmonised) |>
      mutate(count_tissue = log(count_tissue))
  ) |>
  
  
  # Join_intercept
  left_join(
    
    differential_composition_sex_relative |> 
      select(cell_type_harmonised, count_data) |> 
      unnest(count_data) |> 
      distinct() |> 
      with_groups(sample_, ~ .x |> mutate(proportion = count/sum(count))) |> 
      with_groups(c(tissue_harmonised, cell_type_harmonised), ~ .x |> summarise(mean_proportion = mean(proportion))) |> 
      mutate(mean_proportion = pmax(mean_proportion, 1e-7)) |> 
      mutate(mean_proportion_logit = boot::logit(mean_proportion)) |>
      rename(tissue = tissue_harmonised, cell_type = cell_type_harmonised) |> 
      mutate(mean_proportion_logit = (mean_proportion_logit - min(mean_proportion_logit))/4)
    
  ) |> 
  
  # Shorten names
  mutate(cell_type = cell_type |> 
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
plot_heatmap_sex_relative_organ_cell_type =
  
  df_heatmap_sex_relative_organ_cell_type |>
  
  # Filter
  filter(cell_type != "dnt") |> 
  
  # Heatmap
  heatmap(
    tissue, cell_type, Difference,
    palette_value = circlize::colorRamp2(
      seq(3, -3, length.out = 11),
      RColorBrewer::brewer.pal(11, "RdBu")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 0),
    row_title_gp = gpar(fontsize = 0),
    show_heatmap_legend = FALSE
  ) |>
  
  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  annotation_bar(`Mean diff tissue`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  # annotation_tile(
  #   tissue, show_legend = FALSE,
  #   palette =
  #     df_heatmap_sex_relative_organ_cell_type |>
  #     distinct(tissue, tissue_color) |>
  #     arrange(tissue) |>
  #     deframe(),
  #   size = unit(0.2, "cm")
  # ) |>
  annotation_tile(
    cell_type, show_legend = FALSE,
    palette =
      df_heatmap_sex_relative_organ_cell_type |>
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
  layer_point((c_lower * c_upper)>0, .size = sqrt(mean_proportion  ) * 5 )

# Save heatmap separately
plot_heatmap_sex_relative_organ_cell_type |>
  save_pdf(
    filename = glue("{result_directory}/plot_heatmap_sex_relative_organ_cell_type.pdf"),
    width = 80*1.5, height = 60*1.5, units = "mm"
  )


# DE
library(tidyverse)
library(glue)
library(targets)
result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"
store = glue("{result_directory}/_targets__pseudobulk_non_immune_split3")
# Plot of importance of composition vs transcription

library(furrr)
library(tidybulk)
plan(multisession, workers = 36)
options(future.globals.maxSize = 200000 * 1024^2)


# de_sex_cell_type =
# 	tar_meta( store = store	) |>
# 	dplyr::filter(name |> str_detect("estimates_sex_cell_type_")) |>
# 	filter(!is.na(data)) |>
# 	mutate(se = map(
# 		name,
# 		~ .x |>
# 			tar_read_raw(store=store ) |>
# 			mutate(data = map(data, pivot_transcript)),
# 		.progress=T
# 	))
# de_sex_cell_type |> saveRDS("~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_cell_type.rds")


rank_de_cell_type = 
  readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_cell_type.rds") |>
  select(se) |> 
  unnest(se) |>
  filter(!is.na(cell_type_harmonised)) |> 
  unnest(data) |> 
  count(cell_type_harmonised, P_sex_adjusted < 0.05) |> 
  drop_na() |> 
  spread(`P_sex_adjusted < 0.05`, n) |> 
  mutate(proportion_significant = `TRUE` / (`FALSE` + `TRUE`)) |> 
  arrange(desc(proportion_significant)) |>
  rename(cell_type = cell_type_harmonised) |>
  
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
  
  # TEMPORARY FOR OFFLINE
  #scale_fill_manual(values = cell_type_color) +
  guides(fill="none") +
  xlab("Number of significant genes") +
  
  # TEMPORARY FOR OFFLINE
  #theme_multipanel +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Temporary fix
differential_composition_sex_relative = differential_composition_sex_relative |> tidyseurat:::add_attr(FALSE, "check_outliers")

plot_ranks_cell_type = 
  
  differential_composition_sex_relative |> 
  test_contrasts() |>
  filter(parameter=="sexmale") |>
  arrange(c_FDR) |>
  filter(!cell_type_harmonised %in% c("immune_unclassified", "dnt")) |>
  rowid_to_column("rank_composition") |>
  dplyr::select(cell_type = cell_type_harmonised, rank_composition) |>
  
  full_join(
    
    # DE
    rank_de_cell_type ,
    by = join_by(cell_type)
    
  ) |>
  dplyr::select(cell_type, contains("rank")) |>
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
  
  # TEMPORARY FOR OFFLINE
  # scale_color_manual(values = cell_type_color) +
  
  scale_y_reverse() +
  guides(color = "none") +
  ylab("Rank of tissue with highest significance") + 
  
  # TEMPORARY FOR OFFLINE
  #theme_multipanel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 

# library(furrr)
# library(tidybulk)
# plan(multisession, workers = 18)
# options(future.globals.maxSize = 200000 * 1024^2)
# de_sex_tissue =
# 	tar_meta( store = store	) |>
# 	dplyr::filter(name |> str_detect("estimates_sex_tissue_")) |>
# 	filter(!is.na(data)) |>
# 	mutate(se = future_map(
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
# de_sex_tissue |> saveRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue.rds")


de_sex_tissue = readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue.rds")
de_sex_tissue_non_immune = de_sex_tissue |> filter(!is_immune)
de_sex_tissue = de_sex_tissue |> filter(is_immune)


de_sex_tissue = 
  de_sex_tissue |>
  dplyr::rename(tissue = tissue_harmonised) |>
  mutate(tissue = tissue |> str_remove("data_")) |>
  # Parse cell type
  mutate(tissue = case_when(
    tissue == "node" ~ "lymph node",
    tissue == "large" ~ "intestine large",
    tissue == "small" ~ "intestine small",
    tissue == "gland" ~ "adrenal gland",
    TRUE ~ tissue
  )) |>
  filter(map_lgl(data, ~ "P_sex_adjusted" %in% colnames(.x))) |>
  mutate(data = map(data, ~ .x |> select(.feature, P_sex_adjusted, any_of("P_age_days.sex_adjusted"))))

de_sex_tissue_non_immune =
  de_sex_tissue_non_immune |>
  mutate(tissue = tissue_harmonised) |>
  filter(!is.na(tissue)) |>
  filter(map_lgl(data, ~ "P_sex_adjusted" %in% colnames(.x), .progress = TRUE)) |>
  mutate(data = map(data, ~ .x |> select(.feature, P_sex_adjusted, any_of("P_age_days.sex_adjusted"))))

rank_de_tissue = 
  de_sex_tissue |>
  
  # Drop random effects
  filter(tissue |> str_detect("\\:", negate = TRUE)) |> 
  
  select(tissue, data) |> 
  unnest(data) |> 
  left_join(
    de_sex_tissue_non_immune |>
      select(tissue, data) |> 
      unnest(data) |> 
      dplyr::rename(P_sex_adjusted_non_immune = P_sex_adjusted, P_age_days.sex_adjusted_non_immune = P_age_days.sex_adjusted) 
  ) |>
  
  # Filter out sex chromosomes
  filter(!.feature %in% gene_chr$ID) |> 
  
  # Filter our significant in non immune
  filter(P_sex_adjusted_non_immune > 0.05 | is.na(P_sex_adjusted_non_immune)) |> 
  
  # Summary statistics
  dplyr::count(tissue,  P_sex_adjusted < 0.05) |> 
  
  filter(!is.na(`P_sex_adjusted < 0.05`)) |>
  select(tissue,`P_sex_adjusted < 0.05`, n) |>
  spread(`P_sex_adjusted < 0.05`, n) |> 
  mutate(proportion_significant = `TRUE` / (`FALSE` + `TRUE`)) |> 
  arrange(desc(proportion_significant)) |>
  rowid_to_column("rank_de") 

plot_ranks_tissue_barplot = 
  rank_de_tissue |>
  ggplot(aes(`TRUE`, tissue |> fct_reorder(proportion_significant))) +
  geom_bar(aes(fill = tissue), stat = "identity") +
  #scale_fill_manual(values = tissue_color) +
  guides(fill="none") +
  xlab("Number of significant genes") +
  #theme_multipanel +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_ranks_tissue = 
  
  differential_composition_sex_relative |> 
  sccomp_test() |>
  filter(parameter |> str_detect("sexmale___")) |>
  with_groups(parameter, ~ .x |> summarise(median_FDR= median(c_FDR))) |>
  arrange(median_FDR) |>
  
  # Drop random effects
  filter(parameter |> str_detect("\\:", negate = TRUE)) |> 
  
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
  
  
  #scale_color_manual(values = tissue_color) +
  
  scale_y_reverse() +
  guides(color = "none") +
  theme_multipanel



library(ggVennDiagram)
library(eulerr)

de_sex_tissue_for_venn = 
  de_sex_tissue |>
  
  select(tissue, data) |> 
  unnest(data) |> 
  left_join(
    de_sex_tissue_non_immune |>
      select(tissue, data) |> 
      unnest(data) |> 
      dplyr::rename(P_sex_adjusted_non_immune = P_sex_adjusted, P_age_days.sex_adjusted_non_immune = P_age_days.sex_adjusted) 
  ) |>
  
  # Filter out sex chromosomes
  filter(!.feature %in% gene_chr$ID) 
# |> 
# 	
# 	# Filter our significant in non immune
# 	filter(
# 		(P_sex_adjusted_non_immune > 0.05 | is.na(P_sex_adjusted_non_immune)) &
# 			(P_age_days.sex_adjusted_non_immune > 0.05 | is.na(P_age_days.sex_adjusted_non_immune))
# 	)

venn = 
  list(
    tissue_non_immune_sex = 
      de_sex_tissue_non_immune |> 
      filter(map_lgl(data, ~ "P_sex_adjusted" %in% colnames(.x))) |>
      mutate(data = map(data, ~ .x |> select(P_sex_adjusted, .feature))) |>
      select(data, tissue) |>
      
      unnest(data) |>
      add_count(.feature, name = "n_tissue") |>
      filter(n_tissue >= 5) |>
      dplyr::count(P_sex_adjusted<0.05, .feature, n_tissue) |>
      filter(`P_sex_adjusted < 0.05`) |>
      mutate(proportion = n/n_tissue) |>
      #filter(n > 1) |>
      filter(!.feature %in% gene_chr$ID) |> 
      pull(.feature),
    tissue_sex = 
      de_sex_tissue |> 
      filter(map_lgl(data, ~ "P_sex_adjusted" %in% colnames(.x))) |>
      mutate(data = map(data, ~ .x |> select(P_sex_adjusted, .feature))) |>
      select(data, tissue) |>
      
      unnest(data) |>
      add_count(.feature, name = "n_tissue") |>
      filter(n_tissue >= 5) |>
      dplyr::count(P_sex_adjusted<0.05, .feature, n_tissue) |>
      filter(`P_sex_adjusted < 0.05`) |>
      mutate(proportion = n/n_tissue) |>
      filter(n > 2) |>
      filter(!.feature %in% gene_chr$ID) |> 
      pull(.feature),
    tissue_non_immune_sex_age = 
      de_sex_tissue_non_immune |> 
      filter(map_lgl(data, ~ "P_age_days.sex_adjusted" %in% colnames(.x))) |>
      mutate(data = map(data, ~ .x |> select(P_age_days.sex_adjusted, .feature))) |>
      select(data, tissue) |>
      
      unnest(data) |>
      add_count(.feature, name = "n_tissue") |>
      filter(n_tissue >= 5) |>
      dplyr::count(P_age_days.sex_adjusted<0.05, .feature, n_tissue) |>
      filter(`P_age_days.sex_adjusted < 0.05`) |>
      mutate(proportion = n/n_tissue) |>
      #filter(n > 1) |>
      filter(!.feature %in% gene_chr$ID) |> 
      pull(.feature),
    tissue_sex_age = 
      de_sex_tissue |> 
      filter(map_lgl(data, ~ "P_age_days.sex_adjusted" %in% colnames(.x))) |>
      mutate(data = map(data, ~ .x |> select(P_age_days.sex_adjusted, .feature))) |>
      select(data, tissue) |>
      
      unnest(data) |>
      add_count(.feature, name = "n_tissue") |>
      filter(n_tissue >= 5) |>
      dplyr::count(P_age_days.sex_adjusted<0.05, .feature, n_tissue) |>
      filter(`P_age_days.sex_adjusted < 0.05`) |>
      mutate(proportion = n/n_tissue) |>
      filter(n > 2) |>
      filter(!.feature %in% gene_chr$ID) |> 
      pull(.feature)
  ) |>
  euler(fshape = "ellipse") |>
  plot(quantities = TRUE)

# Rank AGE * SEX
rank_de_tissue = 
  de_sex_tissue |>
  
  # Drop random effects
  filter(tissue |> str_detect("\\:", negate = TRUE)) |> 
  
  select(tissue, data) |> 
  unnest(data) |> 
  left_join(
    de_sex_tissue_non_immune |>
      select(tissue, data) |> 
      unnest(data) |> 
      dplyr::rename(P_sex_adjusted_non_immune = P_sex_adjusted, P_age_days.sex_adjusted_non_immune = P_age_days.sex_adjusted) 
  ) |>
  
  # Filter out sex chromosomes
  filter(!.feature %in% gene_chr$ID) |> 
  
  # Filter our significant in non immune
  filter(P_age_days.sex_adjusted_non_immune > 0.05 | is.na(P_age_days.sex_adjusted_non_immune)) |> 
  
  # Summary statistics
  dplyr::count(tissue, significant = P_age_days.sex_adjusted < 0.05 & P_sex_adjusted > 0.05) |> 
  
  filter(!is.na(significant)) |>
  select(tissue,significant, n) |>
  spread(significant, n) |> 
  mutate(proportion_significant = `TRUE` / (`FALSE` + `TRUE`)) |> 
  arrange(desc(proportion_significant)) |>
  rowid_to_column("rank_de") 

plot_ranks_tissue_barplot = 
  rank_de_tissue |>
  ggplot(aes(`TRUE`, tissue |> fct_reorder(proportion_significant))) +
  geom_bar(aes(fill = tissue), stat = "identity") +
  #scale_fill_manual(values = tissue_color) +
  guides(fill="none") +
  xlab("Number of significant genes") +
  #theme_multipanel +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

plot_ranks_tissue_sex_age = 
  
  differential_composition_sex_relative |> 
  sccomp_test() |>
  filter(parameter |> str_detect("sexmale___")) |>
  with_groups(parameter, ~ .x |> summarise(median_FDR= median(c_FDR))) |>
  arrange(median_FDR) |>
  
  # Drop random effects
  filter(parameter |> str_detect("age_days\\:sexmale")) |> 
  
  rowid_to_column("rank_composition") |>
  mutate(tissue = parameter |> str_remove("age_days\\:sexmale___")) |>
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
  
  
  #scale_color_manual(values = tissue_color) +
  
  scale_y_reverse() +
  guides(color = "none") +
  theme_multipanel



# # Which genes are shared across tissue programs
# 
# readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue.rds") |> 
# 	filter(map_lgl(se, ~ "P_sex_adjusted" %in% colnames(.x))) |>
# 	mutate(se = map(se, ~ .x |> select(P_sex_adjusted, .feature))) |>
# 	select(se, name) |>
# 	unnest(se) |>
# 	nest(tissues = -.feature) |>
# 	mutate(tissues_significant = map(
# 		tissues,
# 		~ .x |> filter(P_sex_adjusted<0.05) 
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



plot_genes_most_altered = 
  de_sex_cell_type |>
  

  readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex.rds") |> 
  filter(map_lgl(se, ~ "P_sex_adjusted" %in% colnames(.x))) |>
  mutate(se = map(se, ~ .x |> select(P_sex_adjusted, .feature))) |>
  select(se, name) |>
  unnest(se) |>
  filter(name |> str_detect("immune_unclassified", negate = TRUE)) |>
  nest(cell_types = -.feature) |>
  
  mutate(cell_types_significant = map(
    cell_types,
    ~ .x |> filter(P_sex_adjusted<0.05) 
  )) |>
  
  mutate(
    n_cell_type = map_int(cell_types, nrow),
    n = map_int(cell_types_significant, nrow)
  ) |>
  filter(n_cell_type >= 5) |>
  mutate(proportion = n/n_cell_type) |>
  #filter(proportion > 0.30) |>
  arrange(desc(n)) |>
  filter(!.feature %in% gene_chr$ID) |>
  mutate(cell_type_names = map_chr(cell_types_significant, ~ .x |> pull(name) |> str_c(collapse = ", ") |> str_remove_all("data_"))) |>
  
  # Plot
  mutate(class = case_when(
    row_number() <= 10 ~ "shared",
    row_number() == 11 ~ "CD4",
    row_number() %in% c(21, 23, 34) ~ "naive",
    row_number() == 41 ~ "Monocyte"
  )) |>
  
  mutate(label = if_else(!is.na(class), .feature, NA)) |>
  mutate(n = as.character(n)) |>
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
# 	
# 	se_adjust =
# 		tar_meta(store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue"), starts_with("data_")) |>
# 		filter(!is.na(data)) |>
# 		
# 		mutate(se = future_map(
# 			name,
# 			~ {
# 				se = tar_read_raw(.x, store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue"))
# 				
# 				if(SummarizedExperiment::ncol(se)==0) return(NULL)
# 				if(!"P_sex_adjusted" %in% colnames(SummarizedExperiment::rowData(se))) return(NULL)
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
# 					tidySummarizedExperiment::filter(P_sex_adjusted<0.05) |>
# 					tidySummarizedExperiment::select(
# 						.feature, .sample, sample_,
# 						contains("counts"), contains("sexmale"),
# 						tissue_harmonised, cell_type_harmonised,
# 						sex, P_sex_adjusted, P_sex
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
# 			readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue_non_immune.rds") |>
# 				mutate(tissue = map_chr(se, ~ .x |> pull(tissue_harmonised) |> unique())) |> 
# 				filter(!is.na(tissue)) |>
# 				filter(map_lgl(se, ~ "P_sex_adjusted" %in% colnames(.x), .progress = TRUE)) |>
# 				mutate(de_genes = map(se, ~ .x |> filter( P_sex_adjusted<0.05) |> pull(.feature))) |>
# 				select(tissue, de_genes),
# 			by = "tissue"
# 		) |>
# 	
# 		mutate(effects = map(
# 			se,
# 			~ .x |>
# 				# Significant
# 				filter(P_sex_adjusted<0.05)  |>
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
# 				nest(data = -c(.feature, contains("mode"), sexmale, P_sex_adjusted)) |>
# 				mutate(cell_type_to_select = map(data, ~ .x |> pull(cell_type_harmonised) |> unique())) |>
# 				select(-data) |>
# 				
# 				
# 				pivot_longer(cols = contains("mode"), names_sep = "_cell_type_harmonised__", names_to = c("cell_type", "stat")) |>
# 				filter(map2_lgl(cell_type, cell_type_to_select, ~ .x %in% .y)) |>
# 				with_groups(c(.feature,P_sex_adjusted, sexmale), ~ .x |> summarise(mode_sd = sd(value), n = n())) |> 
# 				arrange(mode_sd) |>
# 				
# 				# Only filter confident
# 				#filter(n>5) |> 
# 				filter(P_sex_adjusted<0.001) |> 
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
# se_adjust |> qs::qsave("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue_adjusted.qs")

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
# 		tar_meta(store = glue("pseudobulk_0.2.3.4/_targets__sex"), starts_with("data_")) |>
# 		filter(!is.na(data)) |>
# 		
# 		mutate(se = future_map(
# 			name,
# 			~ {
# 				se = tar_read_raw(.x, store = glue("pseudobulk_0.2.3.4/_targets__sex"))
# 				
# 				if(SummarizedExperiment::ncol(se)==0) return(NULL)
# 				if(!"P_sex_adjusted" %in% colnames(SummarizedExperiment::rowData(se))) return(NULL)
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
# 					tidySummarizedExperiment::filter(P_sex_adjusted<0.05) |>
# 					tidySummarizedExperiment::select(
# 						.feature, .sample,
# 						contains("counts"), contains("sexmale"),
# 						tissue_harmonised, cell_type_harmonised,
# 						sex, P_sex_adjusted, P_sex
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
# 				filter(P_sex_adjusted<0.05)  |>
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
# 				nest(data = -c(.feature, contains("mode"), sexmale, P_sex_adjusted)) |>
# 				mutate(cell_type_to_select = map(data, ~ .x |> pull(tissue_harmonised) |> unique())) |>
# 				select(-data) |>
# 				
# 				
# 				pivot_longer(cols = contains("mode"), names_sep = "_tissue_harmonised__", names_to = c("cell_type", "stat")) |>
# 				filter(map2_lgl(cell_type, cell_type_to_select, ~ .x %in% .y)) |>
# 				with_groups(c(.feature,P_sex_adjusted, sexmale), ~ .x |> summarise(mode_sd = sd(value), n = n())) |> 
# 				arrange(mode_sd) |>
# 				
# 				# Only filter confident
# 				#filter(n>5) |> 
# 				filter(P_sex_adjusted<0.001) |> 
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
# se_adjust_cell_type |> qs::qsave("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_cell_type_adjusted.qs")




# Select tissue pathways
xx = 
  qs::qread("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue_adjusted.qs") |> 
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
      plot_sex_absolute_1D | plot_sex_absolute_1D_interaction |
        variability_abundance_plot |
        plot_sex_relative
    ) + 
      plot_layout(width = c(0.8, 0.8, 3, 4))
    ) /
      
      # Row 2
      ((
        wrap_heatmap(
          plot_heatmap_sex_relative_organ_cell_type,
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
  glue("{result_directory}/figure_sex.pdf"),
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
my_cells = 
  # Filter interesting tissues for differentil ageing
  glue("{result_directory}/../sccomp_on_HCA_0.2.3.4/input_relative.rds") |> 
  readRDS() |> 
  filter(tissue_harmonised %in% c("thymus", "lymph node", "heart")) |> 
  tidyr::extract(sample_, "sample_", "([a-zA-Z0-9]+)_.+") |>
  mutate(cell_sample = paste(cell_, sample_) ) 
  

library(CuratedAtlasQueryR)
metadata = 
  get_metadata() |>
  mutate(cell_sample = paste(cell_, sample_) ) |> 
  filter(cell_sample %in% !!my_cells$cell_sample) |>
  

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
    ethnicity_simplified, assay_simplified
  ) 

sce = 
  metadata |>
  get_single_cell_experiment() |>
  mutate(sample_se =
           
           # I need to fix Curated CellAtlas with disease sample, duplication for 
           # file_id=="cc3ff54f-7587-49ea-b197-1515b6d98c4c", cell_type_harmonised=="stromal_cell"
           # for lung
           glue("{sample_}___{disease}___{..2}___{..3}") |>
           str_replace_all(" ", "_") |>
           str_replace_all("/", "__")
  ) 




bone_pseudobulk_samples =
  tar_read_raw("data_bone", store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue")) |>
  distinct(sample_) |> 
  extract(col = sample_, into = "sample_", regex = "([a-zA-Z0-9]+)_*") |>
  pull(sample_)


library(DelayedArray)
library(HDF5Array)

# Genes for bone
xx = readRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue.rds") |>
  rename(tissue = name) |>
  mutate(tissue = tissue |> str_remove("data_")) |>
  # Parse cell type
  mutate(tissue = case_when(
    tissue == "node" ~ "lymph node",
    tissue == "large" ~ "intestine large",
    tissue == "small" ~ "intestine small",
    tissue == "gland" ~ "adrenal gland",
    TRUE ~ tissue
  )) |>
  filter(map_lgl(se, ~ "P_sex_adjusted" %in% colnames(.x))) |>
  
  left_join(
    de_sex_tissue_non_immune |>
      select(tissue, se_non_immune = se) 
  ) |>
  
  mutate(se = map2(
    se,se_non_immune,
    ~ {
      .x = .x |>
        filter(!.feature %in% gene_chr$ID) 
      
      if(!is.null(.y)) 	.x = .x |>	filter(!.feature %in% (.y |> pull(.feature))) 
      .x 
      
    }
  )) |> 
  filter(tissue == "bone") |> 
  unnest(se) |> 
  filter(P_sex_adjusted< 0.05) |>
  
  # Save
  arrange(desc(abs(sexmale))) |> 
  select(.feature, sexmale, P_sex_adjusted, contains("sexmale")) |>
  saveRDS("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/de_sex_tissue_bone_for_alex.rds")


get_metadata() |> 
  #	filter(cell_type_harmonised == "cd4 th17") |>
  filter(tissue_harmonised %in% c("bone")) |> 
  filter(sample_ %in% bone_pseudobulk_samples) |>
  add_count(sample_, name = "count_sample") |> 
  filter(count_sample > 200) |> 
  as_tibble()  |> 
  #with_groups(sample_, ~ .x |> sample_n(min(1000, n() ))) |> 
  get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts") |>
  saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/bone")

# kidney
kidney_pseudobulk_samples =
  tar_read_raw("data_kidney", store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue")) |>
  distinct(sample_) |> 
  extract(col = sample_, into = "sample_", regex = "([a-zA-Z0-9]+)_*") |>
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
  saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/kidney")


# adipose
adipose_pseudobulk_samples =
  tar_read_raw("data_adipose", store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue")) |>
  distinct(sample_) |> 
  extract(col = sample_, into = "sample_", regex = "([a-zA-Z0-9]+)_*") |>
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
  saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/adipose")





do.call(
  cbind,
  tar_meta(store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue"), starts_with("se_filtered")) |>
    filter(!is.na(data)) |>
    mutate(se = map(
      name,
      ~ tar_read_raw(.x, store = glue("pseudobulk_0.2.3.4/_targets__sex_tissue")),
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
  with_groups(tissue_harmonised, ~ .x |> count(cell_type) |> summarise(mean_de = median(n)))
count(cell_type, tissue_harmonised) 

de_table |> 
  count(cell_type, tissue_harmonised) |> 
  arrange(desc(n)) |> 
  pull(n) |> hist()

de_table |> 
  count(symbol) |> 
  arrange(desc(n)) |> 
  pull(n) |> hist()

de_table_tissue = 
  de_table |> 
  nest(data = -c(tissue_harmonised)) |> 
  mutate(n_celltype_in_tissue = map_int(data, ~ .x |> distinct(cell_type) |> nrow())) |> 
  mutate(shared_genes = map2(
    data, n_celltype_in_tissue,
    ~ .x |> 
      count(symbol) |> 
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
      count(symbol) |> 
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
ilc |> saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/ilc")


th17 = 
  get_metadata() |> 
  filter(cell_type_harmonised == "cd4 th17") |>
  filter(tissue_harmonised %in% c("lymph node", "heart", "blood", "kidney", "liver", "lung")) |> 
  add_count(sample_, name = "count_sample") |> 
  filter(count_sample > 200) |> 
  as_tibble()  |> 
  with_groups(sample_, ~ .x |> sample_n(min(1000, n() ))) |> 
  get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated", assays = "counts")

th17 |> saveHDF5SummarizedExperiment("sccomp_on_HCA_0.2.3.7_double_interaction_sex_age/th17")


plot_count_de = 
  de_table |> 
  unnest(de) |> 
  filter(adj.P.Val<0.05) |> 
  count(cell_type, tissue_harmonised) |> 
  
  with_groups(cell_type, ~ .x |> mutate(median_n = median(n))) |> 
  
  ggplot(aes(fct_reorder(cell_type, -median_n), n, color = tissue_harmonised)) +
  geom_point(aes()) +
  scale_y_log10() +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=90))


ggMarginal(plot_count_de, groupColour = TRUE, groupFill = TRUE)
