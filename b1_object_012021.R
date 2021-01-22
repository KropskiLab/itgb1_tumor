library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 40000 * 1024^2)
library(RColorBrewer)
library(viridis)
library(hrbrthemes)


#Create Objects
wt_tumor.data <- Read10X_h5("~/Desktop/Haake/WT_tumor/filtered_feature_bc_matrix.h5")
wt_adjacent.data <- Read10X_h5("~/Desktop/Haake/WT_adjacent/filtered_feature_bc_matrix.h5")
B1_tumor.data <- Read10X_h5("~/Desktop/Haake/B1_tumor/filtered_feature_bc_matrix.h5")
B1_adjacent.data <- Read10X_h5("~/Desktop/Haake/B1_adjacent/filtered_feature_bc_matrix.h5")

wt_tumor <- CreateSeuratObject(counts = wt_tumor.data, min.features = 500, project = "WT_Tumor")
wt_adjacent <- CreateSeuratObject(counts = wt_adjacent.data, min.features = 500, project = "WT_Adjacent")
b1_tumor <- CreateSeuratObject(counts = B1_tumor.data, min.features = 500, project = "B1KO_Tumor")
b1_adjacent <- CreateSeuratObject(counts = B1_adjacent.data, min.features = 500, project = "B1KO_Adjacent")

wt_tumor <- PercentageFeatureSet(wt_tumor, pattern = "^mt-", col.name = "percent.mt")
wt_adjacent <- PercentageFeatureSet(wt_adjacent, pattern = "^mt-", col.name = "percent.mt")
b1_tumor <- PercentageFeatureSet(b1_tumor, pattern = "^mt-", col.name = "percent.mt")
b1_adjacent <- PercentageFeatureSet(b1_adjacent, pattern = "^mt-", col.name = "percent.mt")

wt_tumor <- subset(wt_tumor, subset = percent.mt < 10 & percent.mt >0.5)
wt_adjacent <- subset(wt_adjacent, subset = percent.mt < 10 & percent.mt >0.5)  
b1_tumor <- subset(b1_tumor, subset = percent.mt < 10 & percent.mt >0.5)
b1_adjacent <- subset(b1_adjacent, subset = percent.mt < 10 & percent.mt >0.5)

#Merge and Cluster                   
b1_object <- merge(wt_tumor, y = c(wt_adjacent, b1_tumor, b1_adjacent))

Idents(b1_object) <- "orig.ident"
Idents(b1_object, cells = WhichCells(b1_object, idents = c("WT_Tumor", "WT_Adjacent"))) <- "WT"
Idents(b1_object, cells = WhichCells(b1_object, idents = c("B1KO_Tumor", "B1KO_Adjacent"))) <- "ITGB1 KO"
b1_object$genotype <- Idents(b1_object)

Idents(b1_object) <- "orig.ident"
Idents(b1_object, cells = WhichCells(b1_object, idents = c("WT_Tumor", "B1KO_Tumor"))) <- "Tumor"
Idents(b1_object, cells = WhichCells(b1_object, idents = c("WT_Adjacent", "B1KO_Adjacent"))) <- "Adjacent"
b1_object$sample_type <- Idents(b1_object)


b1_object <- SCTransform(b1_object, variable.features.n = 2000, vars.to.regress = c("percent.mt"))
b1_object <- RunPCA(b1_object, verbose = F)
b1_object <- RunUMAP(b1_object, dims = 1:25)
b1_object <- FindNeighbors(b1_object, dims = 1:25)
b1_object <- FindClusters(b1_object, resolution = 0.5)

DimPlot(b1_object, label = T)
DimPlot(b1_object, group.by = "orig.ident")
DotPlot(b1_object, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Krt5', 'Trp63', 'Calca', 'Ascl1', 'Muc5b', 'Mki67'))

#saveRDS(b1_object, file = "~/Desktop/Haake/b1_object.rds")

#Epithelial_subset
epithelial <- subset(b1_object, idents = c(2,19,21,25))
epithelial <- SCTransform(epithelial, vars.to.regress = c("percent.mt"))
epithelial <- RunPCA(epithelial, verbose = F)
epithelial <- RunUMAP(epithelial, dims = 1:12)
epithelial <- FindNeighbors(epithelial, dims = 1:12)
epithelial <- FindClusters(epithelial, resolution = 0.5)
DimPlot(epithelial, label = T)
DimPlot(epithelial, split.by = "orig.ident", label = T)

DotPlot(epithelial, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", 'Sftpc', 'Foxj1', 'Scgb1a1', 'Hopx', 'Krt5', 'Trp63', 'Calca', 'Ascl1', 'Muc5b', 'Mki67'))

#exclude doublets (cluster 6,9)
epithelial <- subset(epithelial, idents = c(0:5,7,8))
epithelial <- SCTransform(epithelial, vars.to.regress = c("percent.mt"))
epithelial <- RunPCA(epithelial, verbose = F)
epithelial <- RunUMAP(epithelial, dims = 1:12)
epithelial <- FindNeighbors(epithelial, dims = 1:12)
epithelial <- FindClusters(epithelial, resolution = 0.5)
DimPlot(epithelial, label = T)
DimPlot(epithelial, split.by = "orig.ident", label = T)

DotPlot(epithelial, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", 'Sftpc', 'Abca3', 'Lamp3', 'Foxj1', 'Scgb1a1', 'Hopx', 'Krt5', 'Trp63', 'Calca', 'Ascl1', 'Muc5b', 'Mki67', "Ctgf", 'Cdkn1a', 'Cldn4'))

epi_markers <- FindAllMarkers(epithelial)

#exclude doublets (cluster 4)
epithelial <- subset(epithelial, idents = c(0:3, 5:9))
epithelial <- SCTransform(epithelial, vars.to.regress = c("percent.mt"))
epithelial <- RunPCA(epithelial, verbose = F)
epithelial <- RunUMAP(epithelial, dims = 1:12)
epithelial <- FindNeighbors(epithelial, dims = 1:12)
epithelial <- FindClusters(epithelial, resolution = 0.8)
DimPlot(epithelial, label = T)
DimPlot(epithelial, split.by = "orig.ident", label = T)

DotPlot(epithelial, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", 'Sftpc', 'Abca3', 'Lamp3', 'Foxj1', 'Scgb1a1', 'Hopx', 'Krt5', 'Trp63', 'Calca', 'Ascl1', 'Muc5b', 'Mki67', "Ctgf", 'Cdkn1a', 'Cldn4'))

#annotated population and cell types
Idents(epithelial) <- "seurat_clusters"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(2,6,9,11,12))) <- "Epithelial"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(0,1,3,4,5,7,8,10))) <- "Tumor"
epithelial$population <- Idents(epithelial)

Idents(epithelial) <- "seurat_clusters"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(12))) <- "Proliferating epithelial"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(11))) <- "Transitional"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(9))) <- "Secretory/Ciliated"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(6))) <- "AT1"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(2))) <- "AT2"
Idents(epithelial, cells = WhichCells(epithelial, idents = c(0,1,3,4,5,7,8,10))) <- "Tumor"
epithelial$celltype <- Idents(epithelial)

saveRDS(epithelial, file = '~/Desktop/Haake/epithelial_annotated_012021.rds')

#stromal cells
DotPlot(b1_object, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", "Cspg4", 'Acta2', 'Mki67'))

#Stromal_object
stromal <- subset(b1_object, idents = c(5,15,30,9,17,26))
stromal <- SCTransform(stromal, vars.to.regress = c("percent.mt"))
stromal <- RunPCA(stromal, verbose = F)
stromal <- RunUMAP(stromal, dims = 1:15)
stromal <- FindNeighbors(stromal, dims = 1:15)
stromal <- FindClusters(stromal, resolution = 0.5)
DimPlot(stromal, label = T)
DimPlot(stromal, split.by = "orig.ident", label = T)

DotPlot(stromal, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", "Cspg4", 'Acta2', 'Wnt2', 'Dcn', 'Wt1', 'Plvap', "Ackr1", 'Hey1', 'Apln', 'Aplnr', 'Car4', 'Ccl21a', 'Mki67'))

#exclude doublets
stromal <- subset(stromal, idents = c(0:4,6,7,9,12,13,16,18))
stromal <- SCTransform(stromal, vars.to.regress = c("percent.mt"))
stromal <- RunPCA(stromal, verbose = F)
stromal <- RunUMAP(stromal, dims = 1:15)
stromal <- FindNeighbors(stromal, dims = 1:15)
stromal <- FindClusters(stromal, resolution = 0.5)
DimPlot(stromal, label = T)
DimPlot(stromal, split.by = "orig.ident", label = T)

DotPlot(stromal, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Col1a1", "Pdgfra", "Cspg4", 'Acta2', 'Wnt2', 'Dcn', 'Pi16', 'Wnt5a', 'Cthrc1', 'Wt1', 'Plvap', "Ackr1", 'Hey1', 'Apln', 'Aplnr', 'Car4', 'Ccl21a', 'Mki67'))

stromal_markers <- FindAllMarkers(stromal)

#annotate cell populations and types
Idents(stromal) <- "seurat_clusters"
Idents(stromal, cells = WhichCells(stromal, idents = c(1,2,3,4,6,8))) <- "Endothelial"
Idents(stromal, cells = WhichCells(stromal, idents = c(0,5,7,9,10,11,12))) <- "Mesenchyme"
stromal$population <- Idents(stromal)

Idents(stromal) <- "seurat_clusters"
Idents(stromal, cells = WhichCells(stromal, idents = c(12))) <- "Lgr5+"
Idents(stromal, cells = WhichCells(stromal, idents = c(11))) <- "Mesothelial"
Idents(stromal, cells = WhichCells(stromal, idents = c(10))) <- "SMC"
Idents(stromal, cells = WhichCells(stromal, idents = c(9))) <- "MyoFB"
Idents(stromal, cells = WhichCells(stromal, idents = c(8))) <- "Venule"
Idents(stromal, cells = WhichCells(stromal, idents = c(7))) <- "Pericyte"
Idents(stromal, cells = WhichCells(stromal, idents = c(5))) <- "Adventitial FB"
Idents(stromal, cells = WhichCells(stromal, idents = c(4))) <- "Capillary - tumor"
Idents(stromal, cells = WhichCells(stromal, idents = c(3))) <- "Capillary - Ca4+"
Idents(stromal, cells = WhichCells(stromal, idents = c(2,6))) <- "Capillary"
Idents(stromal, cells = WhichCells(stromal, idents = c(1))) <- "Arteriole"
Idents(stromal, cells = WhichCells(stromal, idents = c(0))) <- "FB - Wnt2+"
stromal$celltype <- Idents(stromal)

DimPlot(stromal, label = T)

saveRDS(stromal, file = '~/Desktop/Haake/stromal_annotated_012021.rds')

#immune subset
DotPlot(b1_object, features = c("Ptprc", "Epcam", "Pecam1", "Lyz1", "Cd3e", "Ms4a1", "Cd19", "Irf7", 'Jchain', 'Mki67'))

immune <- subset(b1_object, idents = c(0,1,3,4,6,7,8,10,11,12,13,14,16,17,18,20,22,23,24,27,28,29))
immune <- SCTransform(immune, vars.to.regress = c("percent.mt"))
immune <- RunPCA(immune, verbose = F)
immune <- RunUMAP(immune, dims = 1:25)
immune <- FindNeighbors(immune, dims = 1:25)
immune <- FindClusters(immune, resolution = 0.5)
DimPlot(immune, label = T)
DimPlot(immune, split.by = "orig.ident", label = T)

DotPlot(immune, features = c("Ptprc", "Epcam", 'Sftpc', 'Krt19', 'Col1a1', 'Acta2', "Pecam1", "Lyz1", 'Lyz2', 'Cd68', 'Cd14', 'Pparg', 'Cd86', 'Cd3e', 'Cd4', 'Cd8a', 'Il7r', 'Nkg7', 'Gzmb', 'Pdcd1', 'Foxp3', "Ms4a1", "Cd19", "Irf7", 'Jchain', 'Mki67', 'Cpa3' ))

immune_markers <- FindAllMarkers(immune)

#exclude doublets
immune <- subset(immune, idents = c(0:6,8:14,17,18,22,23,24,25))
immune <- SCTransform(immune, vars.to.regress = c("percent.mt"))
immune <- RunPCA(immune, verbose = F)
immune <- RunUMAP(immune, dims = 1:25)
immune <- FindNeighbors(immune, dims = 1:25)
immune <- FindClusters(immune, resolution = 0.5)
DimPlot(immune, label = T)
DimPlot(immune, split.by = "orig.ident", label = T)

DotPlot(immune, features = c("Ptprc", "Epcam", 'Sftpc', 'Krt19', 'Col1a1', 'Acta2', "Pecam1", "Lyz1", 'Lyz2', 'Itgam', 'Itgax', 'Cd68', 'Cd14', 'Pparg', 'Cd86', 'Cd3e', 'Cd4', 'Cd8a', 'Il7r', 'Nkg7', 'Gzmb', 'Pdcd1', 'Foxp3', "Ms4a1", "Cd19", "Irf7", 'Jchain', 'Mki67', 'Cpa3', 'Ly6g' ))

immune_markers2 <- FindAllMarkers(immune)


#annotate cell populations and types
Idents(immune) <- "seurat_clusters"
Idents(immune, cells = WhichCells(immune, idents = c(1,4,5,8,9,12,13,15,16,20,21,22))) <- "Myeloid"
Idents(immune, cells = WhichCells(immune, idents = c(2,3,6,7,10,11,17,18))) <- "T/NK"
Idents(immune, cells = WhichCells(immune, idents = c(0,14,19))) <- "B/Plasma"
immune$population <- Idents(immune)

Idents(immune) <- "seurat_clusters"
Idents(immune, cells = WhichCells(immune, idents = c(22))) <- "Mast"
Idents(immune, cells = WhichCells(immune, idents = c(21))) <- "pDC"
Idents(immune, cells = WhichCells(immune, idents = c(20))) <- "DC - Ccl22+"
Idents(immune, cells = WhichCells(immune, idents = c(19))) <- "Plasma"
Idents(immune, cells = WhichCells(immune, idents = c(17))) <- "CD4 - Trbv19+"
Idents(immune, cells = WhichCells(immune, idents = c(16))) <- "Proliferating immune"
Idents(immune, cells = WhichCells(immune, idents = c(15))) <- "DC - Cd207+"
Idents(immune, cells = WhichCells(immune, idents = c(13))) <- "MDM - Ace+"
Idents(immune, cells = WhichCells(immune, idents = c(11))) <- "CD4 - Pd1+"
Idents(immune, cells = WhichCells(immune, idents = c(10))) <- "T-cell gamma/delta"
Idents(immune, cells = WhichCells(immune, idents = c(9))) <- "Neutrophil"
Idents(immune, cells = WhichCells(immune, idents = c(8))) <- "MDM"
Idents(immune, cells = WhichCells(immune, idents = c(7))) <- "CD8"
Idents(immune, cells = WhichCells(immune, idents = c(6))) <- "NK"
Idents(immune, cells = WhichCells(immune, idents = c(1,5,12))) <- "Macrophage"
Idents(immune, cells = WhichCells(immune, idents = c(4))) <- "Macrophage - Dcstamp+"
Idents(immune, cells = WhichCells(immune, idents = c(3))) <- "CD4 - Icos+/Ctla4+"
Idents(immune, cells = WhichCells(immune, idents = c(2,18))) <- "CD4"
Idents(immune, cells = WhichCells(immune, idents = c(0,14))) <- "B cells"
immune$celltype <- Idents(immune)

DimPlot(immune, label = T, repel = T)

saveRDS(immune, file = '~/Desktop/Haake/immune_annotated_012021.rds')

#Merge Annotated Objects
b1_merged <- merge(epithelial, y = c(immune, stromal))
b1_merged <- SCTransform(b1_merged, vars.to.regress = c("percent.mt"), return.only.var.genes = FALSE)
b1_merged <- RunPCA(b1_merged, verbose = F)
b1_merged <- RunUMAP(b1_merged, dims = 1:35)
DimPlot(b1_merged, label = T)
DimPlot(b1_merged, split.by = "orig.ident")


Idents(b1_merged) <- 'celltype'

DimPlot(b1_merged, split.by = "sample_type")
DimPlot(b1_merged, split.by = "genotype")
DimPlot(b1_merged, group.by = "population")


VlnPlot(b1_merged, features = 'Itgb1', split.by = 'genotype', pt.size = 0.1)

saveRDS(b1_merged, file = '~/Desktop/Haake/b1_merged_annotated_012021.rds')

#Quantify Cell Type Proportions
write.csv((prop.table(table(Idents(b1_merged), b1_merged$orig.ident), margin = 2)), file = "~/Desktop/Haake/cell_type_proportions.csv")
write.csv((prop.table(table(Idents(b1_merged), b1_merged$genotype), margin = 2)), file = "~/Desktop/Haake/cell_type_proportions_by_genotype.csv")


#Define Cell Type Markers
all_markers <- FindAllMarkers(b1_merged)
write.csv(all_markers, file = "~/Desktop/Haake/b1_merged_all_markers.csv")


#Cell Type Differential Expression
tumor <- subset(b1_merged, idents = c("Tumor"))
Idents(tumor) <- "genotype"
tumor_DE <- FindMarkers(tumor, ident.1 = "ITGB1 KO")
write.csv(tumor_DE, file = "~/Desktop/Haake/tumor_DE.csv")

#Cell Type Proportions
write.csv((prop.table(table(Idents(epithelial), epithelial$orig.ident), margin = 2)), file = "~/Desktop/Haake/epithelial_cell_type_proportions.csv")
write.csv((prop.table(table(Idents(stromal), stromal$orig.ident), margin = 2)), file = "~/Desktop/Haake/stromal_cell_type_proportions.csv")
write.csv((prop.table(table(Idents(immune), immune$orig.ident), margin = 2)), file = "~/Desktop/Haake/immune_cell_type_proportions.csv")
write.csv((prop.table(table(Idents(immune), immune$genotype), margin = 2)), file = "~/Desktop/Haake/immune_cell_type_proportions_by_genotype.csv")

write.csv((prop.table(table(Idents(stromal), stromal$genotype), margin = 2)), file = "~/Desktop/Haake/stromal_cell_type_proportions_by_genotype.csv")
write.csv((prop.table(table(Idents(stromal), stromal$orig.ident), margin = 2)), file = "~/Desktop/Haake/stromal_cell_type_proportions_by_origident.csv")


#calculate Itgb1 differential expression
airway <- subset(epithelial, idents = c('Secretory/Ciliated'))
Idents(airway) <- 'genotype'
FindMarkers(airway, features = 'Itgb1', ident.1 = 'ITGB1 KO', test.use = 'negbinom', logfc.threshold = 0.01)

at1 <- subset(epithelial, idents = c('AT1'))
Idents(at1) <- 'genotype'
FindMarkers(at1, features = 'Itgb1', ident.1 = 'ITGB1 KO', test.use = 'negbinom', logfc.threshold = 0.01)

at2 <- subset(epithelial, idents = c('AT2'))
Idents(at2) <- 'genotype'
FindMarkers(at2, features = 'Itgb1', ident.1 = 'ITGB1 KO', test.use = 'negbinom', logfc.threshold = 0.01)

tumor <- subset(epithelial, idents = c('Tumor'))
Idents(tumor) <- 'genotype'
FindMarkers(tumor, features = 'Itgb1', ident.1 = 'ITGB1 KO', test.use = 'negbinom', logfc.threshold = 0.01)

#export as h5ad
library(SeuratDisk)
setwd("~/Desktop/Haake/")

b1_merged <- UpdateSeuratObject(b1_merged)
SaveH5Seurat(b1_merged, filename = "b1_merged.h5Seurat", overwrite = T)
Convert("b1_merged.h5Seurat", dest = "h5ad")

