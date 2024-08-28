
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

Clusters <- readRDS('KJ_Clusters.rds')

levels(Clusters)

options(repr.plot.height = 10, repr.plot.width = 10)
DimPlot(Clusters, reduction = "umap",label.size = 5, label.color = "black", label = 'FALSE')

Clusters$celltype <- Idents(Clusters)

Clusters <- SetIdent(Clusters, value =Clusters$celltype )

table(Clusters@meta.data$celltype)

Clusters$celltype_orig.ident <- paste(Clusters$orig.ident, sep = '_', Idents(Clusters))

Idents(Clusters)<- Clusters$celltype_orig.ident

levels(Clusters)

table(Clusters$celltype_orig.ident)

head(Clusters)

library(ggplot2)

celltype_counts <- table(Clusters$celltype_orig.ident)  # Count occurrences of each cell type
celltype_data <- as.data.frame(celltype_counts)  # Convert to data frame

# Rename the columns for clarity
colnames(celltype_data) <- c("celltype", "count")

# View the first few rows of the data frame
head(celltype_data)

celltype_data

celltype_data <- data.frame(celltype = Clusters$celltype_orig.ident, timepoint = Clusters$orig.ident)

# Aggregate counts for each cell type by time point
aggregated_data <- celltype_data %>%
    group_by(timepoint, celltype) %>%
    summarize(count = n(), .groups = 'drop')  # Count instances

# Calculate total counts per time point
total_counts <- aggregated_data %>%
    group_by(timepoint) %>%
    summarize(total_count = sum(count), .groups = 'drop')

# Join to get total counts for proportion calculation
aggregated_data <- aggregated_data %>%
    left_join(total_counts, by = "timepoint") %>%
    mutate(proportion = count / total_count)  # Calculate proportion

# Create the proportion bar plot
ggplot(aggregated_data, aes(x = reorder(celltype, -proportion), y = proportion, fill = timepoint)) +
    geom_bar(stat = "identity", position = "dodge") +  # Use position = "dodge" for side-by-side bars
    labs(title = "Proportion of Immune Cell Types by Time Point", x = "Cell Type", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Adjust text size and angle for clarity
          axis.text.y = element_text(size = 8))


KJ_Monocytes_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Monocytes','Sham_2_Monocytes','Sham_3_Monocytes','Sham_4_Monocytes'),
                                    ident.1 = c('MI_Day_3_1_Monocytes','MI_Day_3_2_Monocytes','MI_Day_3_3_Monocytes','MI_Day_3_4_Monocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Monocytes_DEG_D3)

KJ_Monocytes_SEG_D3<-subset(KJ_Monocytes_DEG_D3, p_val_adj<0.05)

write.csv(KJ_Monocytes_SEG_D3,"KJ_Monocytes_SEG_D3.csv")

write.csv(KJ_Monocytes_DEG_D3, 'KJ_Monocytes_DEG_D3.csv')

KJ_Monocytes_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Monocytes','Sham_2_Monocytes','Sham_3_Monocytes','Sham_4_Monocytes'),
                                    ident.1 = c('MI_Day_7_1_Monocytes','MI_Day_7_2_Monocytes','MI_Day_7_3_Monocytes','MI_Day_7_4_Monocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Monocytes_DEG_D7)

KJ_Monocytes_SEG_D7<-subset(KJ_Monocytes_DEG_D7, p_val_adj<0.05)

write.csv(KJ_Monocytes_SEG_D7,"KJ_Monocytes_SEG_D7.csv")



write.csv(KJ_Monocytes_DEG_D7, 'KJ_Monocytes_DEG_D7.csv')

KJ_Monocytes_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Monocytes','Sham_2_Monocytes','Sham_3_Monocytes','Sham_4_Monocytes'),
                                    ident.1 = c('MI_Day_14_1_Monocytes','MI_Day_14_2_Monocytes','MI_Day_14_3_Monocytes','MI_Day_14_4_Monocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Monocytes_DEG_D14)

KJ_Monocytes_SEG_D14<-subset(KJ_Monocytes_DEG_D14, p_val_adj<0.05)

write.csv(KJ_Monocytes_SEG_D14, 'KJ_Monocytes_SEG_D14.csv')

write.csv(KJ_Monocytes_DEG_D14, 'KJ_Monocytes_DEG_D14.csv')

KJ_T_NK_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_T_NK_Cells','Sham_2_T_NK_Cells','Sham_3_T_NK_Cells','Sham_4_T_NK_Cells'),
                                    ident.1 = c('MI_Day_3_1_T_NK_Cells','MI_Day_3_2_T_NK_Cells','MI_Day_3_3_T_NK_Cells','MI_Day_3_4_T_NK_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_T_NK_DEG_D3)

write.csv(KJ_T_NK_DEG_D3, 'KJ_T_NK_DEG_D3.csv')

KJ_T_NK_SEG_D3<-subset(KJ_T_NK_DEG_D3, p_val_adj<0.05)

write.csv(KJ_T_NK_SEG_D3, 'KJ_T_NK_SEG_D3.csv')

KJ_T_NK_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_T_NK_Cells','Sham_2_T_NK_Cells','Sham_3_T_NK_Cells','Sham_4_T_NK_Cells'),
                                    ident.1 = c('MI_Day_7_1_T_NK_Cells','MI_Day_7_2_T_NK_Cells','MI_Day_7_3_T_NK_Cells','MI_Day_7_4_T_NK_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_T_NK_DEG_D7)

write.csv(KJ_T_NK_DEG_D7, 'KJ_T_NK_DEG_D7.csv')

KJ_T_NK_SEG_D7<-subset(KJ_T_NK_DEG_D7, p_val_adj<0.05)

write.csv(KJ_T_NK_SEG_D7, 'KJ_T_NK_SEG_D7.csv')

KJ_T_NK_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_T_NK_Cells','Sham_2_T_NK_Cells','Sham_3_T_NK_Cells','Sham_4_T_NK_Cells'),
                                    ident.1 = c('MI_Day_14_1_T_NK_Cells','MI_Day_14_2_T_NK_Cells','MI_Day_14_3_T_NK_Cells','MI_Day_14_4_T_NK_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='T_NK_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_T_NK_DEG_D14)

write.csv(KJ_T_NK_DEG_D14, 'KJ_T_NK_DEG_D14.csv')

KJ_T_NK_SEG_D14<-subset(KJ_T_NK_DEG_D14, p_val_adj<0.05)

write.csv(KJ_T_NK_SEG_D14, 'KJ_T_NK_SEG_D14.csv')

KJ_Macrophages_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Macrophages','Sham_2_Macrophages','Sham_3_Macrophages','Sham_4_Macrophages'),
                                    ident.1 = c('MI_Day_3_1_Macrophages','MI_Day_3_2_Macrophages','MI_Day_3_3_Macrophages','MI_Day_3_4_Macrophages'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Macrophages_DEG_D3)

KJ_Macrophages_SEG_D3<-subset(KJ_Macrophages_DEG_D3,p_val_adj<0.05)

write.csv(KJ_Macrophages_SEG_D3, 'KJ_Macrophages_SEG_D3.csv')

write.csv(KJ_Macrophages_DEG_D3, 'KJ_Macrophages_DEG_D3.csv')

KJ_Macrophages_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Macrophages','Sham_2_Macrophages','Sham_3_Macrophages','Sham_4_Macrophages'),
                                    ident.1 = c('MI_Day_7_1_Macrophages','MI_Day_7_2_Macrophages','MI_Day_7_3_Macrophages','MI_Day_7_4_Macrophages'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Macrophages_DEG_D7)

KJ_Macrophages_SEG_D7<-subset(KJ_Macrophages_DEG_D7,p_val_adj<0.05)

write.csv(KJ_Macrophages_SEG_D7, 'KJ_Macrophages_SEG_D7.csv')

write.csv(KJ_Macrophages_DEG_D7, 'KJ_Macrophages_DEG_D7.csv')

KJ_Macrophages_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Macrophages','Sham_2_Macrophages','Sham_3_Macrophages','Sham_4_Macrophages'),
                                    ident.1 = c('MI_Day_14_1_Macrophages','MI_Day_14_2_Macrophages','MI_Day_14_3_Macrophages','MI_Day_14_4_Macrophages'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Macrophages_DEG_D14)

KJ_Macrophages_SEG_D14<-subset(KJ_Macrophages_DEG_D14, p_val_adj<0.05)

write.csv(KJ_Macrophages_SEG_D14, 'KJ_Macrophages_SEG_D14.csv')

write.csv(KJ_Macrophages_DEG_D14, 'KJ_Macrophages_DEG_D14.csv')

KJ_Dendritic_Cells_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Dendritic_Cells','Sham_2_Dendritic_Cells','Sham_3_Dendritic_Cells','Sham_4_Dendritic_Cells'),
                                    ident.1 = c('MI_Day_3_1_Dendritic_Cells','MI_Day_3_2_Dendritic_Cells','MI_Day_3_3_Dendritic_Cells','MI_Day_3_4_Dendritic_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Dendritic_Cells_DEG_D3)

KJ_Dendritic_Cells_SEG_D3<-subset(KJ_Dendritic_Cells_DEG_D3, p_val_adj<0.05)

write.csv(KJ_Dendritic_Cells_SEG_D3, 'KJ_Dendritic_Cells_SEG_D3.csv')

write.csv(KJ_Dendritic_Cells_DEG_D3, 'KJ_Dendritic_Cells_DEG_D3.csv')

KJ_Dendritic_Cells_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Dendritic_Cells','Sham_2_Dendritic_Cells','Sham_3_Dendritic_Cells','Sham_4_Dendritic_Cells'),
                                    ident.1 = c('MI_Day_7_1_Dendritic_Cells','MI_Day_7_2_Dendritic_Cells','MI_Day_7_3_Dendritic_Cells','MI_Day_7_4_Dendritic_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Dendritic_Cells_DEG_D7)

KJ_Dendritic_Cells_SEG_D7<-subset(KJ_Dendritic_Cells_DEG_D7, p_val_adj<0.05)

write.csv(KJ_Dendritic_Cells_SEG_D7, 'KJ_Dendritic_Cells_SEG_D7.csv')

write.csv(KJ_Dendritic_Cells_DEG_D7, 'KJ_Dendritic_Cells_DEG_D7.csv')

KJ_Dendritic_Cells_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Dendritic_Cells','Sham_2_Dendritic_Cells','Sham_3_Dendritic_Cells','Sham_4_Dendritic_Cells'),
                                    ident.1 = c('MI_Day_14_1_Dendritic_Cells','MI_Day_14_2_Dendritic_Cells','MI_Day_14_3_Dendritic_Cells','MI_Day_14_4_Dendritic_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Dendritic_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Dendritic_Cells_DEG_D14)

KJ_Dendritic_Cells_SEG_D14<-subset(KJ_Dendritic_Cells_DEG_D14, p_val_adj<0.05)

write.csv(KJ_Dendritic_Cells_SEG_D14, 'KJ_Dendritic_Cells_SEG_D14.csv')

write.csv(KJ_Dendritic_Cells_DEG_D14, 'KJ_Dendritic_Cells_DEG_D14.csv')

KJ_B_Cell_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_B_Cells','Sham_2_B_Cells','Sham_3_B_Cells','Sham_4_B_Cells'),
                                    ident.1 = c('MI_Day_3_1_B_Cells','MI_Day_3_2_B_Cells','MI_Day_3_3_B_Cells','MI_Day_3_4_B_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_B_Cell_DEG_D3)

KJ_B_Cell_SEG_D3<-subset(KJ_B_Cell_DEG_D3, p_val_adj<0.05)

write.csv(KJ_B_Cell_SEG_D3, 'KJ_B_Cell_SEG_D3.csv')

write.csv(KJ_B_Cell_DEG_D3, 'KJ_B_Cell_DEG_D3.csv')

KJ_B_Cell_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_B_Cells','Sham_2_B_Cells','Sham_3_B_Cells','Sham_4_B_Cells'),
                                    ident.1 = c('MI_Day_7_1_B_Cells','MI_Day_7_2_B_Cells','MI_Day_7_3_B_Cells','MI_Day_7_4_B_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_B_Cell_DEG_D7)

KJ_B_Cell_SEG_D7<-subset(KJ_B_Cell_DEG_D7, p_val_adj<0.05)

write.csv(KJ_B_Cell_SEG_D7, 'KJ_B_Cell_SEG_D7.csv')

write.csv(KJ_B_Cell_DEG_D7, 'KJ_B_Cell_DEG_D7.csv')

KJ_B_Cell_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_B_Cells','Sham_2_B_Cells','Sham_3_B_Cells','Sham_4_B_Cells'),
                                    ident.1 = c('MI_Day_14_1_B_Cells','MI_Day_14_2_B_Cells','MI_Day_14_3_B_Cells','MI_Day_14_4_B_Cells'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='B_Cells',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_B_Cell_DEG_D14)

KJ_B_Cell_SEG_D14<-subset(KJ_B_Cell_DEG_D14, p_val_adj<0.05)

write.csv(KJ_B_Cell_SEG_D14, 'KJ_B_Cell_SEG_D14.csv')

write.csv(KJ_B_Cell_DEG_D14, 'KJ_B_Cell_DEG_D14.csv')

KJ_Monocytes_Granulocytes_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Monocytes&Granulocytes','Sham_2_Monocytes&Granulocytes','Sham_3_Monocytes&Granulocytes','Sham_4_Monocytes&Granulocytes'),
                                    ident.1 = c('MI_Day_3_1_Monocytes&Granulocytes','MI_Day_3_2_Monocytes&Granulocytes','MI_Day_3_3_Monocytes&Granulocytes','MI_Day_3_4_Monocytes&Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Monocytes_Granulocytes_DEG_D3)

KJ_Monocytes_Granulocytes_SEG_D3<-subset(KJ_Monocytes_Granulocytes_DEG_D3,p_val_adj<0.05 )

write.csv(KJ_Monocytes_Granulocytes_SEG_D3, 'KJ_Monocytes_Granulocytes_SEG_D3.csv')

write.csv(KJ_Monocytes_Granulocytes_DEG_D3, 'KJ_Monocytes_Granulocytes_DEG_D3.csv')

KJ_Monocytes_Granulocytes_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Monocytes&Granulocytes','Sham_2_Monocytes&Granulocytes','Sham_3_Monocytes&Granulocytes','Sham_4_Monocytes&Granulocytes'),
                                    ident.1 = c('MI_Day_7_1_Monocytes&Granulocytes','MI_Day_7_2_Monocytes&Granulocytes','MI_Day_7_3_Monocytes&Granulocytes','MI_Day_7_4_Monocytes&Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Monocytes_Granulocytes_DEG_D7)

KJ_Monocytes_Granulocytes_SEG_D7<-subset(KJ_Monocytes_Granulocytes_DEG_D7, p_val_adj<0.05)

write.csv(KJ_Monocytes_Granulocytes_SEG_D7, 'KJ_Monocytes_Granulocytes_SEG_D7.csv')

write.csv(KJ_Monocytes_Granulocytes_DEG_D7, 'KJ_Monocytes_Granulocytes_DEG_D7.csv')

KJ_Monocytes_Granulocytes_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Monocytes&Granulocytes','Sham_2_Monocytes&Granulocytes','Sham_3_Monocytes&Granulocytes','Sham_4_Monocytes&Granulocytes'),
                                    ident.1 = c('MI_Day_14_1_Monocytes&Granulocytes','MI_Day_14_2_Monocytes&Granulocytes','MI_Day_14_3_Monocytes&Granulocytes','MI_Day_14_4_Monocytes&Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Monocytes&Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Monocytes_Granulocytes_DEG_D14)

KJ_Monocytes_Granulocytes_SEG_D14<-subset(KJ_Monocytes_Granulocytes_DEG_D14, p_val_adj<0.05)

write.csv(KJ_Monocytes_Granulocytes_SEG_D14, 'KJ_Monocytes_Granulocytes_SEG_D14.csv')

write.csv(KJ_Monocytes_Granulocytes_DEG_D14, 'KJ_Monocytes_Granulocytes_DEG_D14.csv')

KJ_Macrophages_Monocytes_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Macrophages&Monocytes','Sham_2_Macrophages&Monocytes','Sham_3_Macrophages&Monocytes','Sham_4_Macrophages&Monocytes'),
                                    ident.1 = c('MI_Day_3_1_Macrophages&Monocytes','MI_Day_3_2_Macrophages&Monocytes','MI_Day_3_3_Macrophages&Monocytes','MI_Day_3_4_Macrophages&Monocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages&Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Macrophages_Monocytes_DEG_D3)

KJ_Macrophages_Monocytes_SEG_D3<-subset(KJ_Macrophages_Monocytes_DEG_D3,p_val_adj<0.05 )

write.csv(KJ_Macrophages_Monocytes_SEG_D3,'KJ_Macrophages_Monocytes_SEG_D3.csv')

write.csv(KJ_Macrophages_Monocytes_DEG_D3, 'KJ_Macrophages_Monocytes_DEG_D3.csv')

KJ_Macrophages_Monocytes_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Macrophages&Monocytes','Sham_2_Macrophages&Monocytes','Sham_3_Macrophages&Monocytes','Sham_4_Macrophages&Monocytes'),
                                    ident.1 = c('MI_Day_7_1_Macrophages&Monocytes','MI_Day_7_2_Macrophages&Monocytes','MI_Day_7_3_Macrophages&Monocytes','MI_Day_7_4_Macrophages&Monocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages&Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Macrophages_Monocytes_DEG_D7)

KJ_Macrophages_Monocytes_SEG_D7<-subset(KJ_Macrophages_Monocytes_DEG_D7,p_val_adj<0.05 )

write.csv(KJ_Macrophages_Monocytes_SEG_D7,'KJ_Macrophages_Monocytes_SEG_D7.csv')

write.csv(KJ_Macrophages_Monocytes_DEG_D7, 'KJ_Macrophages_Monocytes_DEG_D7.csv')

KJ_Macrophages_Monocytes_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Macrophages&Monocytes','Sham_2_Macrophages&Monocytes','Sham_3_Macrophages&Monocytes','Sham_4_Macrophages&Monocytes'),
                                    ident.1 = c('MI_Day_14_1_Macrophages&Monocytes','MI_Day_14_2_Macrophages&Monocytes','MI_Day_14_3_Macrophages&Monocytes','MI_Day_14_4_Macrophages&Monocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Macrophages&Monocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Macrophages_Monocytes_DEG_D14)

KJ_Macrophages_Monocytes_SEG_D14<-subset(KJ_Macrophages_Monocytes_DEG_D14,p_val_adj<0.05 )

write.csv(KJ_Macrophages_Monocytes_SEG_D14, 'KJ_Macrophages_Monocytes_SEG_D14.csv')

write.csv(KJ_Macrophages_Monocytes_DEG_D14, 'KJ_Macrophages_Monocytes_DEG_D14.csv')

KJ_Granulocytes_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Granulocytes','Sham_2_Granulocytes','Sham_3_Granulocytes','Sham_1_Granulocytes'),
                                    ident.1 = c('MI_Day_3_1_Granulocytes','MI_Day_3_2_Granulocytes','MI_Day_3_3_Granulocytes','MI_Day_3_4_Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Granulocytes_DEG_D3)

KJ_Granulocytes_SEG_D3<-subset(KJ_Granulocytes_DEG_D3, p_val_adj<0.05)

write.csv(KJ_Granulocytes_SEG_D3,'KJ_Granulocytes_SEG_D3.csv')

write.csv(KJ_Granulocytes_DEG_D3, 'KJ_Granulocytes_DEG_D3.csv')

KJ_Granulocytes_DEG_D7 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Granulocytes','Sham_2_Granulocytes','Sham_3_Granulocytes','Sham_1_Granulocytes'),
                                    ident.1 = c('MI_Day_7_1_Granulocytes','MI_Day_7_2_Granulocytes','MI_Day_7_3_Granulocytes','MI_Day_7_4_Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Granulocytes_DEG_D7)

KJ_Granulocytes_SEG_D7<-subset(KJ_Granulocytes_DEG_D7, p_val_adj<0.05)

write.csv(KJ_Granulocytes_SEG_D7,'KJ_Granulocytes_SEG_D7.csv')

write.csv(KJ_Granulocytes_DEG_D7, 'KJ_Granulocytes_DEG_D7.csv')

KJ_Granulocytes_DEG_D14 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Granulocytes','Sham_2_Granulocytes','Sham_3_Granulocytes','Sham_1_Granulocytes'),
                                    ident.1 = c('MI_Day_14_1_Granulocytes','MI_Day_14_2_Granulocytes','MI_Day_14_3_Granulocytes','MI_Day_14_4_Granulocytes'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Granulocytes',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Granulocytes_DEG_D14)

KJ_Granulocytes_SEG_D14<-subset(KJ_Granulocytes_DEG_D14, p_val_adj<0.05)

write.csv(KJ_Granulocytes_SEG_D14, 'KJ_Granulocytes_SEG_D14.csv')

write.csv(KJ_Granulocytes_DEG_D14, 'KJ_Granulocytes_DEG_D14.csv')



KJ_Neutrophil_DEG_D3 <- FindMarkers(Clusters,
                                    ident.2 = c('Sham_1_Neutrophils','Sham_2_Neutrophils','Sham_3_Neutrophils','Sham_4_Neutrophils'),
                                    ident.1 = c('MI_Day_3_1_Neutrophils','MI_Day_3_2_Neutrophils','MI_Day_3_3_Neutrophils','MI_Day_3_4_Neutrophils'),
                                    group.by = 'celltype_orig.ident',
                                    celltype ='Neutrophils',
                                    only.pos = FALSE,
                                    logfc.threshold = 0.25)
head(KJ_Neutrophil_DEG_D3)

KJ_Neutrophil_SEG_D3<-subset(KJ_Neutrophil_DEG_D3, p_val_adj<0.05)
