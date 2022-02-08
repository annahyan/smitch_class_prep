library(Seurat)
library(celldex)
library(SingleR)
library(tidyverse)
library(here)


cell_info = read.table(here("data_input/GSE120575_patient_ID_single_cells.txt"),
                        skip = 20, nrows = 16291, sep = "\t")

cell_info = cell_info[, 1:7]

colnames(cell_info) = c("Sample", "CellName", "Source", "Organism",
                        "PatientID", "Response", "Therapy")


# exp.table = read.table(here("data_input/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"),
#                     row.names = 1)
                    
exp.table = read.table(here("data_input/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"),
                    row.names = 1, skip = 2)

cell_names = readLines(con =  here("data_input/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"), n = 1)
cell_names = strsplit(cell_names, "\t")
cell_names = cell_names[[1]]

colnames(exp.table) = cell_names[2:length(cell_names)]

genes = rownames(exp.table)

meta.data = cell_info
rownames(meta.data) = meta.data$CellName


seurat_obj = CreateSeuratObject(counts = exp.table, meta.data = meta.data)

saveRDS(seurat_obj, file = "data_input/GSE120575_data_seurat.RDS")

seurat_obj = readRDS(file = "data_input/GSE120575_data_seurat.RDS")


hpca.se <- HumanPrimaryCellAtlasData()


singleR_annots = SingleR(test = as.SingleCellExperiment(seurat_obj), 
                         ref = hpca.se, 
                       assay.type.test = 1, labels = hpca.se$label.main)


seurat_obj$SingleR_CellType = singleR_annots$pruned.labels
saveRDS(seurat_obj, file = "GSE120575_data_seurat.RDS")

write.table(seurat_obj@meta.data, file = here('metadata.tsv'), 
            sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

write.table(as.matrix(GetAssayData(object = seurat_obj, slot = "counts")), 
            file = here("expressions_matrix.tsv"),
            sep = "\t", col.names = NA, row.names = TRUE)

sample(rownames(seurat_obj), size = 1000, )

pretreatment = seurat_obj@meta.data %>% filter(grepl("Pre", PatientID)) %>% pull(CellName)


### Option 1
seurat_subset = subset(seurat_obj, subset = CellName %in% pretreatment)


### Option 2
patients = c("Pre_P1", "Pre_P2", "Pre_P3", "Pre_P4", 
             "Pre_P6", "Pre_P8", "Pre_P24", "Pre_P35")

seurat_subset = subset(seurat_obj, subset = PatientID %in% patients)


seurat_subset <- NormalizeData(seurat_subset, normalization.method = "LogNormalize", 
                               scale.factor = 10000)

# hpca.se <- HumanPrimaryCellAtlasData()
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

