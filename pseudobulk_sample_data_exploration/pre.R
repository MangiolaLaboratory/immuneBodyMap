if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install("HDF5Array",force = TRUE)
BiocManager::install("SummarizedExperiment")

library(HDF5Array)
library(tidyverse)
library(SummarizedExperiment)
library(uwot)

hdf5 = loadHDF5SummarizedExperiment("pseudobulk_sample_cell_type/") 




metadata_list = hdf5 %>% colData()

macrophage_se <- hdf5[, hdf5$cell_type_consensus_harmonised == 'macrophage']

matrix_mac = macrophage_se %>% assay('counts_scaled') %>% as.matrix()

# Run UMAP on the converted matrix
umap_result <- umap(t(matrix_mac))  # Transpose if rows are genes and columns are cells

# Convert UMAP result to a data frame for plotting
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df = merge(umap_df, metadata)




# Plot UMAP
library(ggplot2)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(alpha = 0.7) +
  labs(title = "UMAP of Scaled Counts for Macrophages", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend since we have a single cell type







