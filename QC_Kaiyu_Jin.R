
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(cowplot)
library(sctransform)
library(ggplot2)
library(SeuratObject)
library(DoubletFinder)
library(fields)
library(parallel)
library(remotes)
library(SeuratDisk)

setwd("/rds/general/user/mg2523/home/apps/OUTS/Kaiyu_Jin")

# Read data
Sham_1 <- Read10X_h5('D_13292569/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_2 <- Read10X_h5('D_13292570/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_3 <- Read10X_h5('D_13292571/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_4 <- Read10X_h5('D_13292572/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_3_1 <- Read10X_h5('D_13292573/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_3_2 <- Read10X_h5('D_13292574/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_3_3 <- Read10X_h5('D_13292575/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_3_4 <- Read10X_h5('D_13292576/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_7_1 <- Read10X_h5('D_13292577/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_7_2 <- Read10X_h5('D_13292578/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_7_3 <- Read10X_h5('D_13292579/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_7_4 <- Read10X_h5('D_13292580/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_14_1 <- Read10X_h5('D_13292581/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_14_2 <- Read10X_h5('D_13292582/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_14_3 <- Read10X_h5('D_13292583/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_Day_14_4 <- Read10X_h5('D_13292584/outs/Cellbender_outs/output_file_seurat.h5', use.names = TRUE, unique.features = TRUE)

Sham_1 <- CreateSeuratObject(counts = Sham_1, project = "Sham_1", min.cells = 3, min.features = 200)
Sham_2 <- CreateSeuratObject(counts = Sham_2, project = "Sham_2", min.cells = 3, min.features = 200)
Sham_3 <- CreateSeuratObject(counts = Sham_3, project = "Sham_3", min.cells = 3, min.features = 200)
Sham_4 <- CreateSeuratObject(counts = Sham_4, project = "Sham_4", min.cells = 3, min.features = 200)
MI_Day_3_1 <- CreateSeuratObject(counts = MI_Day_3_1, project = "MI_Day_3_1", min.cells = 3, min.features = 200)
MI_Day_3_2 <- CreateSeuratObject(counts = MI_Day_3_2, project = "MI_Day_3_2", min.cells = 3, min.features = 200)
MI_Day_3_3 <- CreateSeuratObject(counts = MI_Day_3_3, project = "MI_Day_3_3", min.cells = 3, min.features = 200)
MI_Day_3_4 <- CreateSeuratObject(counts = MI_Day_3_4, project = "MI_Day_3_4", min.cells = 3, min.features = 200)
MI_Day_7_1 <- CreateSeuratObject(counts = MI_Day_7_1, project = "MI_Day_7_1", min.cells = 3, min.features = 200)
MI_Day_7_2 <- CreateSeuratObject(counts = MI_Day_7_2, project = "MI_Day_7_2", min.cells = 3, min.features = 200)
MI_Day_7_3 <- CreateSeuratObject(counts = MI_Day_7_3, project = "MI_Day_7_3", min.cells = 3, min.features = 200)
MI_Day_7_4 <- CreateSeuratObject(counts = MI_Day_7_4, project = "MI_Day_7_4", min.cells = 3, min.features = 200)
MI_Day_14_1 <- CreateSeuratObject(counts = MI_Day_14_1, project = "MI_Day_14_1", min.cells = 3, min.features = 200)
MI_Day_14_2 <- CreateSeuratObject(counts = MI_Day_14_2, project = "MI_Day_14_2", min.cells = 3, min.features = 200)
MI_Day_14_3 <- CreateSeuratObject(counts = MI_Day_14_3, project = "MI_Day_14_3", min.cells = 3, min.features = 200)
MI_Day_14_4 <- CreateSeuratObject(counts = MI_Day_14_4, project = "MI_Day_14_4", min.cells = 3, min.features = 200)


Sham_1 
Sham_2
Sham_3 
Sham_4 
MI_Day_3_1
MI_Day_3_2
MI_Day_3_3 
MI_Day_3_4
MI_Day_7_1
MI_Day_7_2
MI_Day_7_3
MI_Day_7_4
MI_Day_14_1
MI_Day_14_2
MI_Day_14_3
MI_Day_14_4



QC_function=function(x)
    {
        obj=x
        obj$mito.percent <- PercentageFeatureSet(obj, pattern = '^mt-')
        obj <-subset(obj, subset = nCount_RNA > 800 & 
                                nFeature_RNA > 200 &
                                nFeature_RNA < 5000 &
                                mito.percent < 5)
        obj<-NormalizeData(object=obj, normalization.method="LogNormalize", verbose=TRUE, scale.factor=10000)
        obj<-FindVariableFeatures(obj, selection.method="vst", nfeatures=2000)
        obj<-ScaleData(object=obj)
        obj<-RunPCA(object=obj, verbose=TRUE)
        obj<-RunUMAP(obj, dims=1:10, verbose=TRUE)

        return(obj)    
    }

MI_day3_1_QC<-QC_function(MI_Day_3_1)
MI_day3_2_QC<-QC_function(MI_Day_3_2)
MI_day3_3_QC<-QC_function(MI_Day_3_3)
MI_day3_4_QC<-QC_function(MI_Day_3_4)

MI_day7_1_QC<-QC_function(MI_Day_7_1)
MI_day7_2_QC<-QC_function(MI_Day_7_2)
MI_day7_3_QC<-QC_function(MI_Day_7_3)
MI_day7_4_QC<-QC_function(MI_Day_7_4)

MI_day14_1_QC<-QC_function(MI_Day_14_1)
MI_day14_2_QC<-QC_function(MI_Day_14_2)
MI_day14_3_QC<-QC_function(MI_Day_14_3)
MI_day14_4_QC<-QC_function(MI_Day_14_4)

Sham_1_QC<-QC_function(Sham_1)
Sham_2_QC<-QC_function(Sham_2)
Sham_3_QC<-QC_function(Sham_3)
Sham_4_QC<-QC_function(Sham_4)

MI_day3_1_QC
MI_day3_2_QC
MI_day3_3_QC
MI_day3_4_QC

MI_day7_1_QC
MI_day7_2_QC
MI_day7_3_QC
MI_day7_4_QC

MI_day14_1_QC
MI_day14_2_QC
MI_day14_3_QC
MI_day14_4_QC

Sham_1_QC
Sham_2_QC
Sham_3_QC
Sham_4_QC



DF_Function=function(x)
    {
        obj=x
        sweep.res.list <- paramSweep(obj, PCs = 1:20, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
          filter(BCmetric == max(BCmetric)) %>%
          select(pK) 
        pK <- as.numeric(as.character(pK[[1]]))

        ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
        annotations <- obj@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
        nExp_poi <- round(0.076*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        # run doubletFinder 
        obj <- doubletFinder(obj, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)
        
        return(obj)    
    }

MI_day3_1_doublet<-DF_Function(MI_day3_1_QC)
MI_day3_2_doublet<-DF_Function(MI_day3_2_QC)
MI_day3_3_doublet<-DF_Function(MI_day3_3_QC)
MI_day3_4_doublet<-DF_Function(MI_day3_4_QC)

MI_day7_1_doublet<-DF_Function(MI_day7_1_QC)
MI_day7_2_doublet<-DF_Function(MI_day7_2_QC)
MI_day7_3_doublet<-DF_Function(MI_day7_3_QC)
MI_day7_4_doublet<-DF_Function(MI_day7_4_QC)

MI_day14_1_doublet<-DF_Function(MI_day14_1_QC)
MI_day14_2_doublet<-DF_Function(MI_day14_2_QC)
MI_day14_3_doublet<-DF_Function(MI_day14_3_QC)
MI_day14_4_doublet<-DF_Function(MI_day14_4_QC)

Sham_1_doublet<-DF_Function(Sham_1_QC)
Sham_2_doublet<-DF_Function(Sham_2_QC)
Sham_3_doublet<-DF_Function(Sham_3_QC)
Sham_4_doublet<-DF_Function(Sham_4_QC)


DimPlot(MI_day3_1_doublet, reduction='pca', group.by='DF.classifications_0.25_0.25_448', pt.size=1)

head(MI_day3_1_doublet@meta.data)

table(MI_day3_1_doublet@meta.data$DF.classifications_0.25_0.25_448)
MI_day3_1_Object <- subset(MI_day3_1_doublet, DF.classifications_0.25_0.25_448 == 'Singlet')
MI_day3_1_Object

head(MI_day3_2_doublet@meta.data)

DimPlot(MI_day3_2_doublet, reduction='pca', group.by='DF.classifications_0.25_0.17_435', pt.size=1)

table(MI_day3_2_doublet@meta.data$DF.classifications_0.25_0.17_435)
MI_day3_2_Object <- subset(MI_day3_2_doublet, DF.classifications_0.25_0.17_435 == 'Singlet')
MI_day3_2_Object

head(MI_day3_3_doublet@meta.data)

DimPlot(MI_day3_3_doublet, reduction='pca', group.by='DF.classifications_0.25_0.29_451', pt.size=1)

table(MI_day3_3_doublet@meta.data$DF.classifications_0.25_0.29_451)
MI_day3_3_Object <- subset(MI_day3_3_doublet, DF.classifications_0.25_0.29_451 == 'Singlet')
MI_day3_3_Object

head(MI_day3_4_doublet@meta.data)

DimPlot(MI_day3_4_doublet, reduction='pca', group.by='DF.classifications_0.25_0.24_445', pt.size=1)

table(MI_day3_4_doublet@meta.data$DF.classifications_0.25_0.24_445)
MI_day3_4_Object <- subset(MI_day3_4_doublet, DF.classifications_0.25_0.24_445 == 'Singlet')
MI_day3_4_Object

head(MI_day7_1_doublet@meta.data)

DimPlot(MI_day7_1_doublet, reduction='pca', group.by='DF.classifications_0.25_0.005_357', pt.size=1)

table(MI_day7_1_doublet@meta.data$DF.classifications_0.25_0.005_357)
MI_day7_1_Object <- subset(MI_day7_1_doublet, DF.classifications_0.25_0.005_357 == 'Singlet')
MI_day7_1_Object

head(MI_day7_2_doublet@meta.data)

DimPlot(MI_day7_2_doublet, reduction='pca', group.by='DF.classifications_0.25_0.14_348', pt.size=1)

table(MI_day7_2_doublet@meta.data$DF.classifications_0.25_0.14_348)
MI_day7_2_Object <- subset(MI_day7_2_doublet, DF.classifications_0.25_0.14_348 == 'Singlet')
MI_day7_2_Object

head(MI_day7_3_doublet@meta.data)

DimPlot(MI_day7_3_doublet, reduction='pca', group.by='DF.classifications_0.25_0.2_359', pt.size=1)

table(MI_day7_3_doublet@meta.data$DF.classifications_0.25_0.2_359)
MI_day7_3_Object <- subset(MI_day7_3_doublet, DF.classifications_0.25_0.2_359 == 'Singlet')
MI_day7_3_Object

head(MI_day7_4_doublet@meta.data)

DimPlot(MI_day7_4_doublet, reduction='pca', group.by='DF.classifications_0.25_0.005_366', pt.size=1)

table(MI_day7_4_doublet@meta.data$DF.classifications_0.25_0.005_366)
MI_day7_4_Object <- subset(MI_day7_4_doublet, DF.classifications_0.25_0.005_366 == 'Singlet')
MI_day7_4_Object

head(MI_day14_1_doublet@meta.data)

DimPlot(MI_day14_1_doublet, reduction='pca', group.by='DF.classifications_0.25_0.3_209', pt.size=1)

table(MI_day14_1_doublet@meta.data$DF.classifications_0.25_0.3_209)
MI_day14_1_Object <- subset(MI_day14_1_doublet, DF.classifications_0.25_0.3_209 == 'Singlet')
MI_day14_1_Object

head(MI_day14_2_doublet@meta.data)

DimPlot(MI_day14_2_doublet, reduction='pca', group.by='DF.classifications_0.25_0.15_208', pt.size=1)

table(MI_day14_2_doublet@meta.data$DF.classifications_0.25_0.15_208)
MI_day14_2_Object <- subset(MI_day14_2_doublet, DF.classifications_0.25_0.15_208 == 'Singlet')
MI_day14_2_Object

head(MI_day14_3_doublet@meta.data)

DimPlot(MI_day14_3_doublet, reduction='pca', group.by='DF.classifications_0.25_0.27_208', pt.size=1)

table(MI_day14_3_doublet@meta.data$DF.classifications_0.25_0.27_208)
MI_day14_3_Object <- subset(MI_day14_3_doublet, DF.classifications_0.25_0.27_208 == 'Singlet')
MI_day14_3_Object

head(MI_day14_4_doublet@meta.data)

DimPlot(MI_day14_4_doublet, reduction='pca', group.by='DF.classifications_0.25_0.22_205', pt.size=1)

table(MI_day14_4_doublet@meta.data$DF.classifications_0.25_0.22_205)
MI_day14_4_Object <- subset(MI_day14_4_doublet, DF.classifications_0.25_0.22_205 == 'Singlet')
MI_day14_4_Object

head(Sham_1_doublet@meta.data)

DimPlot(Sham_1_doublet, reduction='pca', group.by='DF.classifications_0.25_0.09_167', pt.size=1)

table(Sham_1_doublet@meta.data$DF.classifications_0.25_0.09_167)
Sham_1_Object <- subset(Sham_1_doublet, DF.classifications_0.25_0.09_167 == 'Singlet')
Sham_1_Object

head(Sham_2_doublet@meta.data)

DimPlot(Sham_2_doublet, reduction='pca', group.by='DF.classifications_0.25_0.08_168', pt.size=1)

table(Sham_2_doublet@meta.data$DF.classifications_0.25_0.08_168)
Sham_2_Object <- subset(Sham_2_doublet, DF.classifications_0.25_0.08_168 == 'Singlet')
Sham_2_Object

head(Sham_3_doublet@meta.data)

DimPlot(Sham_3_doublet, reduction='pca', group.by='DF.classifications_0.25_0.07_166', pt.size=1)

table(Sham_3_doublet@meta.data$DF.classifications_0.25_0.07_166)
Sham_3_Object <- subset(Sham_3_doublet, DF.classifications_0.25_0.07_166 == 'Singlet')
Sham_3_Object

head(Sham_4_doublet@meta.data)

DimPlot(Sham_4_doublet, reduction='pca', group.by='DF.classifications_0.25_0.09_166', pt.size=1)

table(Sham_4_doublet@meta.data$DF.classifications_0.25_0.09_166)
Sham_4_Object <- subset(Sham_4_doublet, DF.classifications_0.25_0.09_166 == 'Singlet')
Sham_4_Object



Combined_data=merge(Sham_1_Object, y=c(Sham_2_Object, Sham_3_Object, Sham_4_Object, MI_day3_1_Object, MI_day3_2_Object,MI_day3_3_Object, MI_day3_4_Object, MI_day7_1_Object, MI_day7_2_Object, MI_day7_3_Object, MI_day7_4_Object, MI_day14_1_Object, MI_day14_2_Object, MI_day14_3_Object, MI_day14_4_Object),
                    add.cell.ids=c("Sham_1", "Sham_2", "Sham_3", "Sham_4", "MI_Day_3_1", "MI_Day_3_2","MI_Day_3_3", "MI_Day_3_4","MI_Day_7_1","MI_Day_7_2","MI_Day_7_3","MI_Day_7_4","MI_Day_14_1","MI_Day_14_2","MI_Day_14_3","MI_Day_14_4"), merge.data = TRUE, project="Kaiyu_Jin")

Combined_data

options(repr.plot.height = 5, repr.plot.width = 15)
VlnPlot(Combined_data, features = c("nFeature_RNA", "nCount_RNA", "mito.percent"), group.by="orig.ident")

saveRDS(Combined_data, 'Kaiyu_Jin_Combined_data.rds')

