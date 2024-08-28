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

Kaiyu_Jin_Combined_data <- readRDS('Kaiyu_Jin_Combined_data.rds')

Kaiyu_Jin_Combined_data

Combined_data_transformed<- SCTransform(Kaiyu_Jin_Combined_data)
Combined_data_transformed <- RunPCA(object = Combined_data_transformed)
ElbowPlot(Combined_data_transformed)
Combined_data_transformed <- RunUMAP(Combined_data_transformed, dims = 1:20, reduction = 'pca')

DimPlot(Combined_data_transformed, reduction='pca', group.by = 'orig.ident')

DimPlot(Combined_data_transformed, reduction = 'umap', group.by = 'orig.ident')

# Harmony
Combined_transformed.harmony <- Combined_data_transformed %>%
  RunHarmony(group.by.vars = 'orig.ident',
             plot_convergence = FALSE)
Combined_transformed.harmony@reductions

Combined_transformed.harmony.embed <- Embeddings(Combined_transformed.harmony, "harmony")

Combined_transformed.harmony <-Combined_transformed.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

DimPlot(Combined_transformed.harmony, reduction = "harmony", group.by = 'orig.ident')

Combined_transformed.harmony <- PrepSCTFindMarkers(Combined_transformed.harmony)

cluster.markers <- FindAllMarkers(Combined_transformed.harmony)



# B cells 1,22,9,13
#CellTypist
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Ms4a1","Cd79a","Cd19", "Hbg2", "Hba2", "Ptma", "Malat1", "Hba1","Cd74", "C1qa", "Cd68","Trem2","C1qa", "Ftl", "RNAse1", "Apoe", "Alf1", "Apoc1", "C1qa", "Ptprc","Igfbp1")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

# T_NK Cells #9,11,12,22
#CellTypist, Forte
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Nkg7","Gnly","Cd8a", "Prdm16", "Trdv2", "Snx29P1", "Linc01013", "Ms4a4b", "Gzma")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

#Macrophages
#CellTypist
options(repr.plot.height = 10, repr.plot.width = 15)
markergenes=c("C1qa", "Cd68","Trem2","C1qa", "Ftl", "RNAse1", "Apoe", "Alf1", "Apoc1", "C1qa", "Ptprc","Igfbp1")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

# DC
#CellTypist, Forte
options(repr.plot.height = 10, repr.plot.width = 15)
markergenes=c("Cd1c", "Fcer1a","Clec10a","Fcer1a", "Gpx1", "Rgs1", "Foxh1", "Cap2", "H2ab1", "Cd74")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

#Monocyte
#CellTypist, Forte
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Kit", "Cpa3","Slc45a3","Ctag2", "Adcyap1", "Siglec17p", "Rps26p8", "Page5", "Chil3", "Plac8", "Saa3", "Arg1", "Cux1", "Ctsd", "Add3", "Ly6c2", "Plac8", "Ccr2")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

# Granulocyte
#CellTypist, Forte
options(repr.plot.height = 20, repr.plot.width = 20)
markergenes=c("Lyz", "Fcn1","Tyrobp","Ftl", "Igkc", "Trac", "Neat1", "Sat1", "Lyz", "Cotl1", "Chil3", "Plac8", "Saa3", "Arg1", "S100a8", "S100a9")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

# Neutrophils
#CellTypist
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Lcn2", "Orm1","Mmp8","Page5", "Sele", "Sycp1", "Syn3", "Xcr1", "Defa3", "Ctag2")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

# Endothelial 
#CellTypist
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Cldn5", "Plvap","Sparcl1","Ramp3", "Sele", "Ackr1", "Adgrl4", "Plvap", "Cfi", "Mmrn2", "Meox2", "Fabp4", "Pecam1", "Lyve1", "Cldn5")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

# Fibroblasts 14
#CellTypist
options(repr.plot.height = 30, repr.plot.width = 15)
markergenes=c("Col1a1", "Col1a2","Dcn","Dpt", "Smoc2", "Sfrp1", "Prrx1", "Mxra5", "Ebf2", "Ntrk2", "Olfml1", "Srpx", "Glt8d2", "Gsn", "Dcn", "Wif1", "Dkk3", "Mt2", "Timp1")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

#SMC 19
#Forte
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Rgs5", "Vtn", "Kcnj8", "Myl9")
FeaturePlot(object=Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Combined_transformed.harmony, features=markergenes,ncol = 2 )

new.cluster.ids <- c("Monocytes",'B_Cells','Monocytes&Granulocytes','Macrophages&Monocytes','Macrophages&Monocytes','Monocytes&Granulocytes','Macrophages','Dendritic_Cells',
                     'Macrophages','T_NK_Cells','Macrophages','T_NK_Cells','T_NK_Cells','Granulocytes','Fibroblasts','Granulocytes','Dendritic_Cells','Macrophages','Dendritic_Cells','SMC','Macrophages','Endothelial','Granulocytes','Unknown','Dendritic_Cells')
names(new.cluster.ids) <- levels(Combined_transformed.harmony)
Clusters <- RenameIdents(Combined_transformed.harmony, new.cluster.ids)

saveRDS(Clusters, 'KJ_Clusters.rds')
