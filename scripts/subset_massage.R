library(Seurat)
library(tidyverse)
library(here)
library(SingleR)
library(celldex)

exp_table_subset = read.table(here('data_input/expressions_matrix_subset.tsv'),
                              sep = "\t", row.names = 1, h = T)

metadata_subset = read.table(here('data_input/metadata_subset.tsv'),
                             sep = "\t", row.names = 1, h = T)

seurat_subset = CreateSeuratObject(counts = exp_table_subset, 
                                   meta.data = metadata_subset)


seurat_subset <- NormalizeData(seurat_subset, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

hpca.se <- HumanPrimaryCellAtlasData()
pred.labels<- SingleR(test = as.SingleCellExperiment(seurat_subset), 
                      ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)

pred.labels$anna_labels = pred.labels$pruned.labels
pred.labels$anna_labels[pred.labels$anna_labels %in% c("Macrophage", "Monocyte")] = 
    "Macro/Mono"

pred.labels$anna_labels[! pred.labels$anna_labels %in% 
                            c("Macro/Mono", "B_cell", "NK_cell", "T_cells")] = "Other"

seurat_subset[["CellType"]] = pred.labels$anna_labels

write.table(GetAssayData(seurat_subset, 'count'), 
            file = "data_input/expression_data_normalized.tsv",
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")

write.table(seurat_subset@meta.data, file = "data_input/metadata_subset.tsv",
            quote = FALSE, col.names = NA, row.names = TRUE, sep = "\t")


