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

Combined_Object <- readRDS('Combined_Object.rds')

Combined_Object

Combined_data_transformed<- SCTransform(Combined_Object)
Combined_data_transformed <- RunPCA(object = Combined_data_transformed)
ElbowPlot(Combined_data_transformed)
Combined_data_transformed <- RunUMAP(Combined_data_transformed, dims = 1:20, reduction = 'pca')

DimPlot(Combined_data_transformed, reduction = 'pca', group.by = 'orig.ident')

DimPlot(Combined_data_transformed, reduction = 'umap', group.by = 'orig.ident')



# Harmony
Forte_Combined_transformed.harmony <- Combined_data_transformed %>%
  RunHarmony(group.by.vars = 'orig.ident',
             plot_convergence = FALSE)
Forte_Combined_transformed.harmony@reductions

Forte_Combined_transformed.harmony.embed <- Embeddings(Forte_Combined_transformed.harmony, "harmony")

Forte_Combined_transformed.harmony <- Forte_Combined_transformed.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

DimPlot(Forte_Combined_transformed.harmony, reduction = 'umap', group.by = 'orig.ident')

DimPlot(Forte_Combined_transformed.harmony, reduction = "harmony", group.by = 'orig.ident')

Forte_Combined_transformed.harmony <- PrepSCTFindMarkers(Forte_Combined_transformed.harmony)

Forte_Combined_transformed.harmony@meta.data

cluster.markers <- FindAllMarkers(Forte_Combined_transformed.harmony)

# B cells #9
#CellTypist
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Ms4a1","Cd79a","Cd19", "Hbg2", "Hba2", "Ptma", "Malat1", "Hba1","Cd74", "C1qa", "Cd68","Trem2","C1qa", "Ftl", "RNAse1", "Apoe", "Alf1", "Apoc1", "C1qa", "Ptprc","Igfbp1")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2 )

# T_NK Cells #12
#CellTypist, Forte
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Nkg7","Gnly","Cd8a", "Prdm16", "Trdv2", "Snx29P1", "Linc01013", "Ms4a4b", "Gzma")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2 )

#Macrophages #5,13
#CellTypist
options(repr.plot.height = 10, repr.plot.width = 15)
markergenes=c("C1qa", "Cd68","Trem2","C1qa", "Ftl", "RNAse1", "Apoe", "Alf1", "Apoc1", "C1qa", "Ptprc","Igfbp1")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2 )

#Monocyte #8
#CellTypist, Forte
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Kit", "Cpa3","Slc45a3","Ctag2", "Adcyap1", "Siglec17p", "Rps26p8", "Page5", "Chil3", "Plac8", "Saa3", "Arg1", "Cux1", "Ctsd", "Add3", "Ly6c2", "Plac8", "Ccr2")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

# DC #17
#CellTypist, Forte
options(repr.plot.height = 10, repr.plot.width = 15)
markergenes=c("Cd1c", "Fcer1a","Clec10a","Fcer1a", "Gpx1", "Rgs1", "Foxh1", "Cap2", "H2ab1", "Cd74")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

# Granulocyte #11
#CellTypist, Forte
options(repr.plot.height = 20, repr.plot.width = 20)
markergenes=c("Lyz", "Fcn1","Tyrobp","Ftl", "Igkc", "Trac", "Neat1", "Sat1", "Lyz", "Cotl1", "Chil3", "Plac8", "Saa3", "Arg1", "S100a8", "S100a9")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

# Neutrophils
#CellTypist
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Lcn2", "Orm1","Mmp8","Page5", "Sele", "Sycp1", "Syn3", "Xcr1", "Defa3", "Ctag2")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

# Endothelial #4,14
#CellTypist
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Cldn5", "Plvap","Sparcl1","Ramp3", "Sele", "Ackr1", "Adgrl4", "Plvap", "Cfi", "Mmrn2", "Meox2", "Fabp4", "Pecam1", "Lyve1", "Cldn5")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

# Fibroblasts #0,1,2,6,7,10,15,16,19
#CellTypist
options(repr.plot.height = 30, repr.plot.width = 15)
markergenes=c("Col1a1", "Col1a2","Dcn","Dpt", "Smoc2", "Sfrp1", "Prrx1", "Mxra5", "Ebf2", "Ntrk2", "Olfml1", "Srpx", "Glt8d2", "Gsn", "Dcn", "Wif1", "Dkk3", "Mt2", "Timp1")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

#SMC #18
#Forte
options(repr.plot.height = 15, repr.plot.width = 15)
markergenes=c("Rgs5", "Vtn", "Kcnj8", "Myl9")
FeaturePlot(object=Forte_Combined_transformed.harmony, reduction = "umap", label = TRUE, features=markergenes, order=TRUE, ncol = 2)
VlnPlot(object=Forte_Combined_transformed.harmony, features=markergenes,ncol = 2)

new.cluster.ids <- c("Fibroblasts",'Fibroblasts','Fibroblasts','Monocytes&Granulocytes','Endothelial','Macrophages','Fibroblasts','Fibroblasts',
                     'Monocytes','B_Cells','Fibroblasts','Granulocytes','T_NK_Cells','Macrophages&Monocytes','Endothelial','Fibroblasts','Fibroblasts','Dendritic_Cells','SMC','Fibroblasts')
names(new.cluster.ids) <- levels(Forte_Combined_transformed.harmony)
Forte_Clusters <- RenameIdents(Forte_Combined_transformed.harmony, new.cluster.ids)

DimPlot(Forte_Clusters, reduction = "umap", label = TRUE, pt.size = 0.5)

saveRDS(Forte_Clusters, 'Forte_Clusters.rds')

