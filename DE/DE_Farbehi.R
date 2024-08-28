library(Seurat)
library(SeuratObject)
library(harmony)
library(dplyr)
library(patchwork)
library(cowplot)
library(sctransform)
library(ggplot2)
library(DoubletFinder)
library(fields)
library(parallel)
library(remotes)
library(SeuratDisk)

Clusters <- readRDS('RP_Clusters.rds')

levels(Clusters)

options(repr.plot.height = 12, repr.plot.width = 13)
DimPlot(Clusters, reduction = "umap", label.size = 4, pt.size =0.5)

Clusters$celltype <- Idents(Clusters)

Clusters <- SetIdent(Clusters, value =Clusters$celltype )

table(Clusters@meta.data$celltype)

Clusters$celltype_orig.ident <- paste(Clusters$orig.ident, sep = '_', Idents(Clusters))

Idents(Clusters)<- Clusters$celltype_orig.ident

table(Clusters@meta.data$celltype_orig.ident)

levels(Clusters)

RP_B_Cells_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_B_Cells','MI_7TIP_Sham_B_Cells'),
                                    ident.1 = c('MI_3TIP_B_Cells','MI_3GFP_B_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(RP_B_Cells_DEG_D3)

RP_B_SEG_D3<-subset(RP_B_Cells_DEG_D3, p_val_adj<0.05 )

RP_B_SEG_D3

write.csv(RP_B_SEG_D3, 'RP_B_SEG_D3.csv')

write.csv(RP_B_Cells_DEG_D3, 'RP_B_Cells_DEG_D3.csv')

RP_B_Cells_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_B_Cells','MI_7TIP_Sham_B_Cells'),
                                    ident.1 = c('MI_7TIP_B_Cells','MI_7GFP_B_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(RP_B_Cells_DEG_D7)

RP_B_SEG_D7<-subset(RP_B_Cells_DEG_D7, p_val_adj<0.05 )

write.csv(RP_B_SEG_D7, 'RP_B_SEG_D7.csv')

write.csv(RP_B_Cells_DEG_D7, 'RP_B_Cells_DEG_D7.csv')





RP_Dendritic_Cells_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_Dendritic_Cells','MI_7TIP_Sham_Dendritic_Cells'),
                                    ident.1 = c('MI_3TIP_Dendritic_Cells','MI_3GFP_Dendritic_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(RP_Dendritic_Cells_DEG_D3)

RP_Dendritic_Cells_SEG_D3<-subset(RP_Dendritic_Cells_DEG_D3, p_val_adj<0.05 )

write.csv(RP_Dendritic_Cells_SEG_D3, 'RP_Dendritic_Cells_SEG_D3.csv')

write.csv(RP_Dendritic_Cells_DEG_D3, 'RP_Dendritic_Cells_DEG_D3.csv')

RP_Dendritic_Cells_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_Dendritic_Cells','MI_7TIP_Sham_Dendritic_Cells'),
                                    ident.1 = c('MI_7TIP_Dendritic_Cells','MI_7GFP_Dendritic_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(RP_Dendritic_Cells_DEG_D7)

RP_Dendritic_Cells_SEG_D7<-subset(RP_Dendritic_Cells_DEG_D7, p_val_adj<0.05 )

write.csv(RP_Dendritic_Cells_SEG_D7, 'RP_Dendritic_Cells_SEG_D7.csv')

write.csv(RP_Dendritic_Cells_DEG_D7, 'RP_Dendritic_Cells_DEG_D7.csv')

RP_Monocytes_Macrophages_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_Monocytes&Macrophages','MI_7TIP_Sham_Monocytes&Macrophages'),
                                    ident.1 = c('MI_3TIP_Monocytes&Macrophages','MI_3GFP_Monocytes&Macrophages'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)

head(RP_Monocytes_Macrophages_DEG_D3)

RP_Monocytes_Macrophages_SEG_D3<-subset(RP_Monocytes_Macrophages_DEG_D3, p_val_adj<0.05 )

write.csv(RP_Monocytes_Macrophages_SEG_D3, 'RP_Monocytes_Macrophages_SEG_D3.csv')

write.csv(RP_Monocytes_Macrophages_DEG_D3, 'RP_Monocytes_Macrophages_DEG_D3.csv')

RP_Monocytes_Macrophages_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_Monocytes&Macrophages','MI_7TIP_Sham_Monocytes&Macrophages'),
                                    ident.1 = c('MI_7TIP_Monocytes&Macrophages','MI_7GFP_Monocytes&Macrophages'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(RP_Monocytes_Macrophages_DEG_D7)

RP_Monocytes_Macrophages_SEG_D7<-subset(RP_Monocytes_Macrophages_DEG_D7, p_val_adj<0.05 )

write.csv(RP_Monocytes_Macrophages_SEG_D7, 'RP_Monocytes_Macrophages_SEG_D7.csv')

write.csv(RP_Monocytes_Macrophages_DEG_D7, 'RP_Monocytes_Macrophages_DEG_D7.csv')

RP_Granulocytes_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_Granulocytes','MI_7TIP_Sham_Granulocytes'),
                                    ident.1 = c('MI_3TIP_Granulocytes','MI_3GFP_Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)

head(RP_Granulocytes_DEG_D3)

RP_Granulocytes_SEG_D3<-subset(RP_Granulocytes_DEG_D3, p_val_adj<0.05 )

write.csv(RP_Granulocytes_SEG_D3, 'RP_Granulocytes_SEG_D3.csv')

write.csv(RP_Granulocytes_DEG_D3, 'RP_Granulocytes_DEG_D3.csv')

RP_Granulocytes_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('MI_7GFP_Sham_Granulocytes','MI_7TIP_Sham_Granulocytes'),
                                    ident.1 = c('MI_7GFP_Granulocytes','MI_7TIP_Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)

head(RP_Granulocytes_DEG_D7)

write.csv(RP_Granulocytes_DEG_D7, 'RP_Granulocytes_DEG_D7.csv')

RP_Granulocytes_SEG_D7<-subset(RP_Granulocytes_DEG_D7, p_val_adj<0.05 )

write.csv(RP_Granulocytes_SEG_D7, 'RP_Granulocytes_SEG_D7.csv')
