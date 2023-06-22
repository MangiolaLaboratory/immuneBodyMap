

library(tidyverse)
library(patchwork)



metadata = 
	"/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.rds" |> 
	readRDS() |> 

	# Attach lineage
	left_join(read_csv("~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv") |> replace_na(list(lineage_1 = "other_non_immune"))) |>
	mutate(is_immune = lineage_1 == "immune")  |> 
	filter(is_immune) |> 
	filter(cell_type_harmonised |> str_detect("eryth|epidermal|Lang|thymo", negate = TRUE))

# Load cell type colors
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")
cell_type_color = 
	metadata |> 
	pull(cell_type_harmonised) |> 
	unique() |> 
	get_cell_type_color()

# Theme ggplot
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")


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

plot_count_confidence =
	metadata |>
	filter(cell_type_harmonised != "immune_unclassified") |> 
	mutate( confidence_class = if_else( confidence_class==5, 4,  confidence_class)) |> 
	# mutate(confidence_class = if_else(confidence_class == 4 | is.na(confidence_class), "Low", "High")) |>
	
	ggplot(aes(confidence_class, fill = cell_type_harmonised)) +
	geom_bar() +
	scale_fill_manual(values = cell_type_color) +
	scale_x_reverse() +
	coord_flip() +
	xlab("Confidence class") +
	ylab("Cell count") +
	guides(fill="none") +
	theme_multipanel

plot_proportion_confidence_ethnicity =
	metadata |>
	filter(is_immune) |> 
	filter(confidence_class<5) |> 
	mutate(ethnicity = ethnicity |>
				 	str_remove("or Latin American") |> 
				 	str_replace("African American", "AfroAmerican") |>
				 	str_remove("AfroAmerican or") |> 
				 	str_trim()
				) |> 
	select(`cell_`, ethnicity, cell_type_harmonised, confidence_class) |>
	as_tibble() |>
	dplyr::count(confidence_class, ethnicity) |>
	with_groups(ethnicity, ~ .x |> mutate(proportion = n/(sum(n)))) |>
	
	nest(data = -ethnicity) |>
	mutate(proportion_level_1 = map_dbl(data, ~ .x |> filter(confidence_class ==1) |> pull(proportion) )) |>
	unnest(data) |>

	filter(!is.na(confidence_class)) |>
	ggplot(aes(fct_reorder(ethnicity,1- proportion_level_1), proportion, fill =confidence_class)) +
	geom_bar(stat = "identity") +
	scale_fill_viridis_c(direction=-1) +
	xlab("Ethnicity") +
	guides(fill="none") +
	theme_multipanel +
	theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

plot_proportion_confidence_assay =
	metadata |>
	filter(is_immune) |> 
	filter(confidence_class<5) |> 
	mutate(assay = assay |> str_remove("transcription profiling") |> str_remove("Spatial Gene Expression") |> str_remove("technology") |> str_trim()) |> 
	filter(assay |> str_detect("ATAC", negate = TRUE)) |> 
	select(`cell_`, assay, cell_type_harmonised, confidence_class) |>
	as_tibble() |>
	dplyr::count(confidence_class, assay) |>
	with_groups(assay, ~ .x |> mutate(proportion = n/(sum(n)))) |>
	
	nest(data = -assay) |>
	mutate(proportion_level_1 = map_dbl(data, ~ .x |> filter(confidence_class ==1) |> pull(proportion) )) |>
	unnest(data) |>
	
	filter(!is.na(confidence_class)) |>
	ggplot(aes(fct_reorder(assay,1- proportion_level_1), proportion, fill =confidence_class)) +
	geom_bar(stat = "identity") +
	scale_fill_viridis_c(direction=-1) +
	xlab("Technology") +
	guides(fill="none") +
	theme_multipanel +
	theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


p = 
	plot_count_confidence | (plot_proportion_confidence_assay / plot_proportion_confidence_ethnicity ) +
	plot_layout( guides = 'collect', widths = c(52, 87) ) &
	theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")

ggsave(
	"~/PostDoc/immuneHealthyBodyMap/sccomp_on_HCA_0.2.3/figure_curation.pdf",
	plot = p,
	units = c("mm"),
	width = 135 ,
	height = 50 ,
	limitsize = FALSE
)
