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

Forte_Clusters <- readRDS('Forte_Clusters.rds')

levels(Forte_Clusters)

head(Forte_Clusters@meta.data)

Forte_Clusters$celltype <- Idents(Forte_Clusters)

Forte_Clusters <- SetIdent(Forte_Clusters, value =Forte_Clusters$celltype )

table(Forte_Clusters@meta.data$celltype)

head(Forte_Clusters@meta.data)

Forte_Clusters$celltype_orig.ident <- paste(Forte_Clusters$orig.ident, sep = '_', Idents(Forte_Clusters))

Idents(Forte_Clusters)<- Forte_Clusters$celltype_orig.ident

levels(Forte_Clusters)

table(Forte_Clusters@meta.data$celltype_orig.ident)


Monocytes_Granulocytes_D1 <- FindMarkers(Forte_Clusters,
                                    ident.1 = 'day_1_Monocytes&Granulocytes',
                                    ident.2 = c('Sham_day_7_1_a_Monocytes&Granulocytes','Sham_day_7_1_b_Monocytes&Granulocytes','Sham_day_7_2_a_Monocytes&Granulocytes','Sham_day_7_2_b_Monocytes&Granulocytes','Sham_day_7_3_Monocytes&Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monocytes_Granulocytes_D1)



Monocytes_Granulocytes_SEG_D1<-subset(Monocytes_Granulocytes_D1,p_val_adj<0.05 )

write.csv(Monocytes_Granulocytes_SEG_D1, 'Monocytes&Granulocytes/Monocytes_Granulocytes_SEG_D1.csv')

write.csv(Monocytes_Granulocytes_D1, 'Monocytes_Granulocytes_D1.csv')

Monocytes_Granulocytes_D3 <- FindMarkers(Forte_Clusters,
                                    ident.1 = 'day_3_Monocytes&Granulocytes',
                                    ident.2 = c('Sham_day_7_1_a_Monocytes&Granulocytes','Sham_day_7_1_b_Monocytes&Granulocytes','Sham_day_7_2_a_Monocytes&Granulocytes','Sham_day_7_2_b_Monocytes&Granulocytes','Sham_day_7_3_Monocytes&Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monocytes_Granulocytes_D3)

Monocytes_Granulocytes_SEG_D3<-subset(Monocytes_Granulocytes_D3,p_val_adj<0.05 )

write.csv(Monocytes_Granulocytes_SEG_D3, 'Monocytes&Granulocytes/Monocytes_Granulocytes_SEG_D3.csv')

write.csv(Monocytes_Granulocytes_D3, 'Monocytes&Granulocytes/Monocytes_Granulocytes_D3.csv')

Monocytes_Granulocytes_D5 <- FindMarkers(Forte_Clusters,
                                    ident.1 = c('Sham_day_7_1_a_Monocytes&Granulocytes','Sham_day_7_1_b_Monocytes&Granulocytes','Sham_day_7_2_a_Monocytes&Granulocytes','Sham_day_7_2_b_Monocytes&Granulocytes','Sham_day_7_3_Monocytes&Granulocytes'),
                                    ident.2 = 'day_5_Monocytes&Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monocytes_Granulocytes_D5)

Monocytes_Granulocytes_SEG_D5<-subset(Monocytes_Granulocytes_D5,p_val_adj<0.05 )

write.csv(Monocytes_Granulocytes_SEG_D5, 'Monocytes&Granulocytes/Monocytes_Granulocytes_SEG_D5.csv')

write.csv(Monocytes_Granulocytes_D5, 'Monocytes&Granulocytes/Monocytes_Granulocytes_D5.csv')

Monocytes_Granulocytes_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes&Granulocytes','Sham_day_7_1_b_Monocytes&Granulocytes','Sham_day_7_2_a_Monocytes&Granulocytes','Sham_day_7_2_b_Monocytes&Granulocytes','Sham_day_7_3_Monocytes&Granulocytes'),
                                    ident.1 = 'day_7_Monocytes&Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monocytes_Granulocytes_D7)

Monocytes_Granulocytes_SEG_D7<-subset(Monocytes_Granulocytes_D7,p_val_adj<0.05 )

write.csv(Monocytes_Granulocytes_SEG_D7, 'Monocytes&Granulocytes/Monocytes_Granulocytes_SEG_D7.csv')

write.csv(Monocytes_Granulocytes_D7, 'Monocytes&Granulocytes/Monocytes_Granulocytes_D7.csv')

Monocytes_Granulocytes_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes&Granulocytes','Sham_day_7_1_b_Monocytes&Granulocytes','Sham_day_7_2_a_Monocytes&Granulocytes','Sham_day_7_2_b_Monocytes&Granulocytes','Sham_day_7_3_Monocytes&Granulocytes'),
                                    ident.1 = 'day_14_Monocytes&Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monocytes_Granulocytes_D14)

Monocytes_Granulocytes_SEG_D14<-subset(Monocytes_Granulocytes_D14,p_val_adj<0.05 )

write.csv(Monocytes_Granulocytes_SEG_D14, 'Monocytes&Granulocytes/Monocytes_Granulocytes_SEG_D14.csv')

write.csv(Monocytes_Granulocytes_D14, 'Monocytes&Granulocytes/Monocytes_Granulocytes_D14.csv')

Monocytes_Granulocytes_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes&Granulocytes','Sham_day_7_1_b_Monocytes&Granulocytes','Sham_day_7_2_a_Monocytes&Granulocytes','Sham_day_7_2_b_Monocytes&Granulocytes','Sham_day_7_3_Monocytes&Granulocytes'),
                                    ident.1 = 'day_28_Monocytes&Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monocytes_Granulocytes_D28)

Monocytes_Granulocytes_SEG_D28<-subset(Monocytes_Granulocytes_D28,p_val_adj<0.05 )

write.csv(Monocytes_Granulocytes_SEG_D28, 'Monocytes&Granulocytes/Monocytes_Granulocytes_SEG_D28.csv')

write.csv(Monocytes_Granulocytes_D28, 'Monocytes&Granulocytes/Monocytes_Granulocytes_D28.csv')



Macro_DEG_D1 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Macrophages','Sham_day_7_1_b_Macrophages','Sham_day_7_2_a_Macrophages','Sham_day_7_2_b_Macrophages','Sham_day_7_3_Macrophages'),
                                    ident.1 = 'day_1_Macrophages',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Macro_DEG_D1)

write.csv(Macro_DEG_D1, 'Macrophages/Macro_DEG_D1.csv')

Macrophages_SEG_D1<-subset(Macro_DEG_D1,p_val_adj<0.05 )

write.csv(Macrophages_SEG_D1, 'Macrophages/Macrophages_SEG_D1.csv')



Macro_DEG_D3 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Macrophages','Sham_day_7_1_b_Macrophages','Sham_day_7_2_a_Macrophages','Sham_day_7_2_b_Macrophages','Sham_day_7_3_Macrophages'),
                                    ident.1 = 'day_3_Macrophages',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Macro_DEG_D3)

Macrophages_SEG_D3<-subset(Macro_DEG_D3, p_val_adj<0.05 )

write.csv(Macrophages_SEG_D1, 'Macrophages/Macrophages_SEG_D3.csv')

write.csv(Macro_DEG_D3, 'Macrophages/Macro_DEG_D3.csv')



Macro_DEG_D5 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Macrophages','Sham_day_7_1_b_Macrophages','Sham_day_7_2_a_Macrophages','Sham_day_7_2_b_Macrophages','Sham_day_7_3_Macrophages'),
                                    ident.1 = 'day_5_Macrophages',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Macro_DEG_D5)

Macrophages_SEG_D5<-subset(Macro_DEG_D5, p_val_adj<0.05 )

write.csv(Macrophages_SEG_D5, 'Macrophages/Macrophages_SEG_D5.csv')

write.csv(Macro_DEG_D5, 'Macrophages/Macro_DEG_D5.csv')



Macro_DEG_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Macrophages','Sham_day_7_1_b_Macrophages','Sham_day_7_2_a_Macrophages','Sham_day_7_2_b_Macrophages','Sham_day_7_3_Macrophages'),
                                    ident.1 = 'day_7_Macrophages',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Macro_DEG_D7)

Macrophages_SEG_D7<-subset(Macro_DEG_D7, p_val_adj<0.05 )

write.csv(Macrophages_SEG_D7, 'Macrophages/Macrophages_SEG_D7.csv')

write.csv(Macro_DEG_D7, 'Macrophages/Macro_DEG_D7.csv')



Macro_DEG_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Macrophages','Sham_day_7_1_b_Macrophages','Sham_day_7_2_a_Macrophages','Sham_day_7_2_b_Macrophages','Sham_day_7_3_Macrophages'),
                                    ident.1 = 'day_14_Macrophages',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Macro_DEG_D14)

Macrophages_SEG_D14<-subset(Macro_DEG_D14, p_val_adj<0.05 )

write.csv(Macrophages_SEG_D14, 'Macrophages/Macrophages_SEG_D14.csv')

write.csv(Macro_DEG_D14, 'Macro_DEG_D14.csv')

Macro_DEG_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Macrophages','Sham_day_7_1_b_Macrophages','Sham_day_7_2_a_Macrophages','Sham_day_7_2_b_Macrophages','Sham_day_7_3_Macrophages'),
                                    ident.1 = 'day_28_Macrophages',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Macro_DEG_D28)

write.csv(Macro_DEG_D28, 'Macrophages/Macro_DEG_D28.csv')

Macrophages_SEG_D28<-subset(Macro_DEG_D28, p_val_adj<0.05 )

write.csv(Macrophages_SEG_D28, 'Macrophages/Macrophages_SEG_D28.csv')

Monnocyte_DEG_D1 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes','Sham_day_7_1_b_Monocytes','Sham_day_7_2_a_Monocytes','Sham_day_7_2_b_Monocytes','Sham_day_7_3_Monocytes'),
                                    ident.1 = 'day_1_Monocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monnocyte_DEG_D1)

write.csv(Monnocyte_DEG_D1, 'Monocyte/Monnocyte_DEG_D1.csv')

Monocyte_SEG_D1<-subset(Monnocyte_DEG_D1, p_val_adj<0.05 )

write.csv(Monocyte_SEG_D1, 'Monocyte/Monocyte_SEG_D1.csv')

Monnocyte_DEG_D3 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes','Sham_day_7_1_b_Monocytes','Sham_day_7_2_a_Monocytes','Sham_day_7_2_b_Monocytes','Sham_day_7_3_Monocytes'),
                                    ident.1 = 'day_3_Monocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monnocyte_DEG_D3)

write.csv(Monnocyte_DEG_D3, 'Monocyte/Monnocyte_DEG_D3.csv')

Monocyte_SEG_D3<-subset(Monnocyte_DEG_D3, p_val_adj<0.05 )

write.csv(Monocyte_SEG_D3, 'Monocyte/Monocyte_SEG_D3.csv')

Monnocyte_DEG_D5 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes','Sham_day_7_1_b_Monocytes','Sham_day_7_2_a_Monocytes','Sham_day_7_2_b_Monocytes','Sham_day_7_3_Monocytes'),
                                    ident.1 = 'day_5_Monocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monnocyte_DEG_D5)

write.csv(Monnocyte_DEG_D5, 'Monnocyte_DEG_D5.csv')

Monocyte_SEG_D5<-subset(Monnocyte_DEG_D5, p_val_adj<0.05 )

write.csv(Monocyte_SEG_D5, 'Monocyte/Monocyte_SEG_D5.csv')

Monnocyte_DEG_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes','Sham_day_7_1_b_Monocytes','Sham_day_7_2_a_Monocytes','Sham_day_7_2_b_Monocytes','Sham_day_7_3_Monocytes'),
                                    ident.1 = 'day_7_Monocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monnocyte_DEG_D7)

write.csv(Monnocyte_DEG_D7, 'Monnocyte_DEG_D7.csv')

Monocyte_SEG_D7<-subset(Monnocyte_DEG_D7, p_val_adj<0.05 )

write.csv(Monocyte_SEG_D7, 'Monocyte/Monocyte_SEG_D7.csv')

Monnocyte_DEG_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes','Sham_day_7_1_b_Monocytes','Sham_day_7_2_a_Monocytes','Sham_day_7_2_b_Monocytes','Sham_day_7_3_Monocytes'),
                                    ident.1 = 'day_14_Monocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monnocyte_DEG_D14)

write.csv(Monnocyte_DEG_D14, 'Monnocyte_DEG_D14.csv')

Monocyte_SEG_D14<-subset(Monnocyte_DEG_D14, p_val_adj<0.05 )

write.csv(Monocyte_SEG_D14, 'Monocyte/Monocyte_SEG_D14.csv')

Monnocyte_DEG_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Monocytes','Sham_day_7_1_b_Monocytes','Sham_day_7_2_a_Monocytes','Sham_day_7_2_b_Monocytes','Sham_day_7_3_Monocytes'),
                                    ident.1 = 'day_28_Monocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Monnocyte_DEG_D28)

write.csv(Monnocyte_DEG_D28, 'Monnocyte_DEG_D28.csv')

Monocyte_SEG_D28<-subset(Monnocyte_DEG_D28, p_val_adj<0.05 )

write.csv(Monocyte_SEG_D28, 'Monocyte/Monocyte_SEG_D28.csv')

Granulocyte_DEG_D1 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Granulocytes','Sham_day_7_1_b_Granulocytes','Sham_day_7_2_a_Granulocytes','Sham_day_7_2_b_Granulocytes','Sham_day_7_3_Granulocytes'),
                                    ident.1 = 'day_1_Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Granulocyte_DEG_D1)

write.csv(Granulocyte_DEG_D1, 'Granulocyte/Granulocyte_DEG_D1.csv')

Granulocyte_SEG_D1<-subset(Granulocyte_DEG_D1, p_val_adj<0.05 )

write.csv(Granulocyte_SEG_D1, 'Granulocyte/Granulocyte_SEG_D1.csv')

Granulocyte_DEG_D3 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Granulocytes','Sham_day_7_1_b_Granulocytes','Sham_day_7_2_a_Granulocytes','Sham_day_7_2_b_Granulocytes','Sham_day_7_3_Granulocytes'),
                                    ident.1 = 'day_3_Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Granulocyte_DEG_D3)

Granulocyte_SEG_D3<-subset(Granulocyte_DEG_D3, p_val_adj<0.05 )

write.csv(Granulocyte_SEG_D3, 'Granulocyte/Granulocyte_SEG_D3.csv')

write.csv(Granulocyte_DEG_D3, 'Granulocyte/Granulocyte_DEG_D3.csv')

Granulocyte_DEG_D5 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Granulocytes','Sham_day_7_1_b_Granulocytes','Sham_day_7_2_a_Granulocytes','Sham_day_7_2_b_Granulocytes','Sham_day_7_3_Granulocytes'),
                                    ident.1 = 'day_5_Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Granulocyte_DEG_D5)

Granulocyte_SEG_D5<-subset(Granulocyte_DEG_D5, p_val_adj<0.05 )

write.csv(Granulocyte_SEG_D5, 'Granulocyte/Granulocyte_SEG_D5.csv')

write.csv(Granulocyte_DEG_D5, 'Granulocyte/Granulocyte_DEG_D5.csv')

Granulocyte_DEG_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Granulocytes','Sham_day_7_1_b_Granulocytes','Sham_day_7_2_a_Granulocytes','Sham_day_7_2_b_Granulocytes','Sham_day_7_3_Granulocytes'),
                                    ident.1 = 'day_7_Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Granulocyte_DEG_D7)

Granulocyte_SEG_D7<-subset(Granulocyte_DEG_D7, p_val_adj<0.05 )

write.csv(Granulocyte_SEG_D7, 'Granulocyte/Granulocyte_SEG_D7.csv')

write.csv(Granulocyte_DEG_D7, 'Granulocyte/Granulocyte_DEG_D7.csv')

Granulocyte_DEG_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Granulocytes','Sham_day_7_1_b_Granulocytes','Sham_day_7_2_a_Granulocytes','Sham_day_7_2_b_Granulocytes','Sham_day_7_3_Granulocytes'),
                                    ident.1 = 'day_14_Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Granulocyte_DEG_D14)

write.csv(Granulocyte_DEG_D14, 'Granulocyte/Granulocyte_DEG_D14.csv')

Granulocyte_SEG_D14<-subset(Granulocyte_DEG_D14, p_val_adj<0.05 )

write.csv(Granulocyte_SEG_D14, 'Granulocyte/Granulocyte_SEG_D14.csv')

Granulocyte_DEG_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_Granulocytes','Sham_day_7_1_b_Granulocytes','Sham_day_7_2_a_Granulocytes','Sham_day_7_2_b_Granulocytes','Sham_day_7_3_Granulocytes'),
                                    ident.1 = 'day_28_Granulocytes',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Granulocyte_DEG_D28)

Granulocyte_SEG_D28<-subset(Granulocyte_DEG_D28, p_val_adj<0.05 )

write.csv(Granulocyte_SEG_D28, 'Granulocyte/Granulocyte_SEG_D28.csv')

write.csv(Granulocyte_DEG_D28, 'Granulocyte/Granulocyte_DEG_D28.csv')

T_NK_Cells_DEG_D1 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_T_NK_Cells','Sham_day_7_1_b_T_NK_Cells','Sham_day_7_2_a_T_NK_Cells','Sham_day_7_2_b_T_NK_Cells','Sham_day_7_3_T_NK_Cells'),
                                    ident.1 = 'day_1_T_NK_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(T_NK_Cells_DEG_D1)

write.csv(T_NK_Cells_DEG_D1, 'T_NK_Cells/T_NK_Cells_DEG_D1.csv')

T_NK_Cells_SEG_D1<-subset(T_NK_Cells_DEG_D1, p_val_adj<0.05 )

write.csv(T_NK_Cells_SEG_D1, 'T_NK_Cells/T_NK_Cells_SEG_D1.csv')

T_NK_Cells_DEG_D3 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_T_NK_Cells','Sham_day_7_1_b_T_NK_Cells','Sham_day_7_2_a_T_NK_Cells','Sham_day_7_2_b_T_NK_Cells','Sham_day_7_3_T_NK_Cells'),
                                    ident.1 = 'day_3_T_NK_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(T_NK_Cells_DEG_D3)

write.csv(T_NK_Cells_DEG_D3, 'T_NK_Cells/T_NK_Cells_DEG_D3.csv')

T_NK_Cells_SEG_D3<-subset(T_NK_Cells_DEG_D3, p_val_adj<0.05 )

write.csv(T_NK_Cells_SEG_D3, 'T_NK_Cells/T_NK_Cells_SEG_D3.csv')

T_NK_Cells_DEG_D5 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_T_NK_Cells','Sham_day_7_1_b_T_NK_Cells','Sham_day_7_2_a_T_NK_Cells','Sham_day_7_2_b_T_NK_Cells','Sham_day_7_3_T_NK_Cells'),
                                    ident.1 = 'day_5_T_NK_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(T_NK_Cells_DEG_D5)

write.csv(T_NK_Cells_DEG_D5, 'T_NK_Cells/T_NK_Cells_DEG_D5.csv')

T_NK_Cells_SEG_D5<-subset(T_NK_Cells_DEG_D5, p_val_adj<0.05 )

write.csv(T_NK_Cells_SEG_D5, 'T_NK_Cells/T_NK_Cells_SEG_D5.csv')

T_NK_Cells_DEG_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_T_NK_Cells','Sham_day_7_1_b_T_NK_Cells','Sham_day_7_2_a_T_NK_Cells','Sham_day_7_2_b_T_NK_Cells','Sham_day_7_3_T_NK_Cells'),
                                    ident.1 = 'day_7_T_NK_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(T_NK_Cells_DEG_D7)

write.csv(T_NK_Cells_DEG_D7, 'T_NK_Cells/T_NK_Cells_DEG_D7.csv')

T_NK_Cells_SEG_D7<-subset(T_NK_Cells_DEG_D7, p_val_adj<0.05 )

write.csv(T_NK_Cells_SEG_D7, 'T_NK_Cells/T_NK_Cells_SEG_D7.csv')

T_NK_Cells_DEG_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_T_NK_Cells','Sham_day_7_1_b_T_NK_Cells','Sham_day_7_2_a_T_NK_Cells','Sham_day_7_2_b_T_NK_Cells','Sham_day_7_3_T_NK_Cells'),
                                    ident.1 = 'day_14_T_NK_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(T_NK_Cells_DEG_D14)

write.csv(T_NK_Cells_DEG_D14, 'T_NK_Cells_DEG_D14.csv')

T_NK_Cells_SEG_D14<-subset(T_NK_Cells_DEG_D14, p_val_adj<0.05 )

write.csv(T_NK_Cells_SEG_D14, 'T_NK_Cells/T_NK_Cells_SEG_D14.csv')

T_NK_Cells_DEG_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_T_NK_Cells','Sham_day_7_1_b_T_NK_Cells','Sham_day_7_2_a_T_NK_Cells','Sham_day_7_2_b_T_NK_Cells','Sham_day_7_3_T_NK_Cells'),
                                    ident.1 = 'day_28_T_NK_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(T_NK_Cells_DEG_D28)

write.csv(T_NK_Cells_DEG_D28, 'T_NK_Cells_DEG_D28.csv')

T_NK_Cells_SEG_D28<-subset(T_NK_Cells_DEG_D28, p_val_adj<0.05 )

write.csv(T_NK_Cells_SEG_D28, 'T_NK_Cells/T_NK_Cells_SEG_D28.csv')

B_Cells_DEG_D1 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_B_Cells','Sham_day_7_1_b_B_Cells','Sham_day_7_2_a_B_Cells','Sham_day_7_2_b_B_Cells','Sham_day_7_3_B_Cells'),
                                    ident.1 = 'day_1_B_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(B_Cells_DEG_D1)

B_D1_SG <- subset(B_Cells_DEG_D1, p_val_adj < 0.05)

write.csv(B_D1_SG, 'B_Cells/B_D1_SG.csv')

write.csv(B_Cells_DEG_D1, 'B_Cells/B_Cells_DEG_D1.csv')

B_Cells_DEG_D3 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_B_Cells','Sham_day_7_1_b_B_Cells','Sham_day_7_2_a_B_Cells','Sham_day_7_2_b_B_Cells','Sham_day_7_3_B_Cells'),
                                    ident.1 = 'day_3_B_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(B_Cells_DEG_D3)

write.csv(B_Cells_DEG_D3, 'B_Cells/B_Cells_DEG_D3.csv')

B_D3_SG <- subset(B_Cells_DEG_D3, p_val_adj < 0.05)

write.csv(B_D3_SG, 'B_Cells/B_D3_SG.csv')



B_Cells_DEG_D5 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_B_Cells','Sham_day_7_1_b_B_Cells','Sham_day_7_2_a_B_Cells','Sham_day_7_2_b_B_Cells','Sham_day_7_3_B_Cells'),
                                    ident.1 = 'day_5_B_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(B_Cells_DEG_D5)

B_D5_SG <- subset(B_Cells_DEG_D5, p_val_adj < 0.05)

write.csv(B_D5_SG, 'B_Cells/B_D5_SG.csv')

write.csv(B_Cells_DEG_D5, 'B_Cells/B_Cells_DEG_D5.csv')



B_Cells_DEG_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_B_Cells','Sham_day_7_1_b_B_Cells','Sham_day_7_2_a_B_Cells','Sham_day_7_2_b_B_Cells','Sham_day_7_3_B_Cells'),
                                    ident.1 = 'day_7_B_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(B_Cells_DEG_D7)

B_D7_SG <- subset(B_Cells_DEG_D7, p_val_adj < 0.05)

write.csv(B_D7_SG, 'B_Cells/B_D7_SG.csv')

write.csv(B_Cells_DEG_D7, 'B_Cells/B_Cells_DEG_D7.csv')



B_Cells_DEG_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_B_Cells','Sham_day_7_1_b_B_Cells','Sham_day_7_2_a_B_Cells','Sham_day_7_2_b_B_Cells','Sham_day_7_3_B_Cells'),
                                    ident.1 = 'day_14_B_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(B_Cells_DEG_D14)

B_D14_SG <- subset(B_Cells_DEG_D14, p_val_adj < 0.05)

write.csv(B_D14_SG, 'B_Cells/B_D14_SG.csv')

write.csv(B_Cells_DEG_D14, 'B_Cells/B_Cells_DEG_D14.csv')



B_Cells_DEG_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 = c('Sham_day_7_1_a_B_Cells','Sham_day_7_1_b_B_Cells','Sham_day_7_2_a_B_Cells','Sham_day_7_2_b_B_Cells','Sham_day_7_3_B_Cells'),
                                    ident.1 = 'day_28_B_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(B_Cells_DEG_D28)

B_D28_SG <- subset(B_Cells_DEG_D28, p_val_adj < 0.05)

write.csv(B_D28_SG, 'B_Cells/B_D28_SG.csv')

write.csv(B_Cells_DEG_D28, 'B_Cells/B_Cells_DEG_D28.csv')

Dendritic_Cells_DEG_D1 <- FindMarkers(Forte_Clusters,
                                    ident.2 =  c('Sham_day_7_1_a_Dendritic_Cells','Sham_day_7_1_b_Dendritic_Cells','Sham_day_7_2_a_Dendritic_Cells','Sham_day_7_2_b_Dendritic_Cells','Sham_day_7_3_Dendritic_Cells'),
                                    ident.1 = 'day_1_Dendritic_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Dendritic_Cells_DEG_D1)

write.csv(Dendritic_Cells_DEG_D1, 'Dendritic_Cells/Dendritic_Cells_DEG_D1.csv')

Dendritic_Cells_SEG_D1 <- subset(Dendritic_Cells_DEG_D1, p_val_adj < 0.05)

write.csv(Dendritic_Cells_SEG_D1, 'Dendritic_Cells/Dendritic_Cells_SEG_D1.csv')

Dendritic_Cells_DEG_D3 <- FindMarkers(Forte_Clusters,
                                    ident.2 =  c('Sham_day_7_1_a_Dendritic_Cells','Sham_day_7_1_b_Dendritic_Cells','Sham_day_7_2_a_Dendritic_Cells','Sham_day_7_2_b_Dendritic_Cells','Sham_day_7_3_Dendritic_Cells'),
                                    ident.1 = 'day_3_Dendritic_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Dendritic_Cells_DEG_D3)

Dendritic_Cells_SEG_D3 <- subset(Dendritic_Cells_DEG_D3, p_val_adj < 0.05)

write.csv(Dendritic_Cells_SEG_D3, 'Dendritic_Cells/Dendritic_Cells_SEG_D3.csv')

write.csv(Dendritic_Cells_DEG_D3, 'Dendritic_Cells/Dendritic_Cells_DEG_D3.csv')

Dendritic_Cells_DEG_D5 <- FindMarkers(Forte_Clusters,
                                    ident.2 =  c('Sham_day_7_1_a_Dendritic_Cells','Sham_day_7_1_b_Dendritic_Cells','Sham_day_7_2_a_Dendritic_Cells','Sham_day_7_2_b_Dendritic_Cells','Sham_day_7_3_Dendritic_Cells'),
                                    ident.1 = 'day_5_Dendritic_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Dendritic_Cells_DEG_D5)

Dendritic_Cells_SEG_D5 <- subset(Dendritic_Cells_DEG_D5, p_val_adj < 0.05)

write.csv(Dendritic_Cells_SEG_D5, 'Dendritic_Cells/Dendritic_Cells_SEG_D5.csv')

write.csv(Dendritic_Cells_DEG_D5, 'Dendritic_Cells/Dendritic_Cells_DEG_5.csv')

Dendritic_Cells_DEG_D7 <- FindMarkers(Forte_Clusters,
                                    ident.2 =  c('Sham_day_7_1_a_Dendritic_Cells','Sham_day_7_1_b_Dendritic_Cells','Sham_day_7_2_a_Dendritic_Cells','Sham_day_7_2_b_Dendritic_Cells','Sham_day_7_3_Dendritic_Cells'),
                                    ident.1 = 'day_7_Dendritic_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Dendritic_Cells_DEG_D7)

Dendritic_Cells_SEG_D7 <- subset(Dendritic_Cells_DEG_D7, p_val_adj < 0.05)

write.csv(Dendritic_Cells_SEG_D7, 'Dendritic_Cells/Dendritic_Cells_SEG_D7.csv')

write.csv(Dendritic_Cells_DEG_D7, 'Dendritic_Cells/Dendritic_Cells_DEG_7.csv')

Dendritic_Cells_DEG_D14 <- FindMarkers(Forte_Clusters,
                                    ident.2 =  c('Sham_day_7_1_a_Dendritic_Cells','Sham_day_7_1_b_Dendritic_Cells','Sham_day_7_2_a_Dendritic_Cells','Sham_day_7_2_b_Dendritic_Cells','Sham_day_7_3_Dendritic_Cells'),
                                    ident.1 = 'day_14_Dendritic_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Dendritic_Cells_DEG_D14)

Dendritic_Cells_SEG_D14 <- subset(Dendritic_Cells_DEG_D14, p_val_adj < 0.05)

write.csv(Dendritic_Cells_SEG_D14, 'Dendritic_Cells/Dendritic_Cells_SEG_D14.csv')

write.csv(Dendritic_Cells_DEG_D14, 'Dendritic_Cells/Dendritic_Cells_DEG_14.csv')

Dendritic_Cells_DEG_D28 <- FindMarkers(Forte_Clusters,
                                    ident.2 =  c('Sham_day_7_1_a_Dendritic_Cells','Sham_day_7_1_b_Dendritic_Cells','Sham_day_7_2_a_Dendritic_Cells','Sham_day_7_2_b_Dendritic_Cells','Sham_day_7_3_Dendritic_Cells'),
                                    ident.1 = 'day_28_Dendritic_Cells',
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(Dendritic_Cells_DEG_D28)

Dendritic_Cells_SEG_D28 <- subset(Dendritic_Cells_DEG_D28, p_val_adj < 0.05)

write.csv(Dendritic_Cells_SEG_D28, 'Dendritic_Cells/Dendritic_Cells_SEG_D28.csv')

write.csv(Dendritic_Cells_DEG_D2, 'Dendritic_Cells/Dendritic_Cells_DEG_28.csv')
