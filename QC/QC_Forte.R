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

#Read in the data
MI_day1_data <- Read10X_h5('cellranger_outs_12040/outs/Cellbender_outs/output_12040_filtered_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_day3_data <- Read10X_h5('cellranger_outs_12041/outs/Cellbender_outs/output_12041_filtered.seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_day5_data <- Read10X_h5('cellranger_outs_12042/outs/Cellbender_outs/output_12042_filtered.seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_day7_data <- Read10X_h5('cellranger_outs_12043/outs/Cellbender_outs/output_12043_filtered.seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_day14_data <- Read10X_h5('cellranger_outs_15642/outs/Cellbender_outs/output_15642_filtered.seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_day28_data <- Read10X_h5('cellranger_outs_15643/outs/Cellbender_outs/output_15643_filtered.seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_day7_data_1_a <- Read10X_h5('cellranger_outs_10999/outs/Cellbender_outs/output_10999_filtered_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_day7_data_1_b<- Read10X_h5('cellranger_outs_13609/outs/Cellbender_outs/output_13609_filtered_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_day7_data_2_a<- Read10X_h5('cellranger_outs_11003/outs/Cellbender_outs/11003_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_day7_data_2_b <- Read10X_h5('cellranger_outs_13613/outs/Cellbender_outs/output_13613_filtered_seurat.h5', use.names = TRUE, unique.features = TRUE)
Sham_day7_data_3 <- Read10X_h5('cellranger_outs_03499/outs/Cellbender_outs/output_03499_seurat.h5', use.names = TRUE, unique.features = TRUE)

#Creating a Seurat Objects for each read in data
MI_day1_data <- CreateSeuratObject(counts = MI_day1_data, project = "day_1", min.cells = 3, min.features = 200)
MI_day3_data <- CreateSeuratObject(counts = MI_day3_data, project = "day_3", min.cells = 3, min.features = 200)
MI_day5_data <- CreateSeuratObject(counts = MI_day5_data, project = "day_5", min.cells = 3, min.features = 200)
MI_day7_data <- CreateSeuratObject(counts = MI_day7_data, project = "day_7", min.cells = 3, min.features = 200)
MI_day14_data <- CreateSeuratObject(counts = MI_day14_data, project = "day_14", min.cells = 3, min.features = 200)
MI_day28_data <- CreateSeuratObject(counts = MI_day28_data, project = "day_28", min.cells = 3, min.features = 200)

MI_day1_data
MI_day3_data
MI_day5_data
MI_day7_data
MI_day14_data
MI_day28_data

Sham_day7_data_1_a <- CreateSeuratObject(counts = Sham_day7_data_1_a, project = "Sham_day_7_1_a", min.cells = 3, min.features = 200)
Sham_day7_data_1_b <- CreateSeuratObject(counts = Sham_day7_data_1_b, project = "Sham_day_7_1_b", min.cells = 3, min.features = 200)
Sham_day7_data_2_a <- CreateSeuratObject(counts = Sham_day7_data_2_a, project = "Sham_day_7_2_a", min.cells = 3, min.features = 200)
Sham_day7_data_2_b <- CreateSeuratObject(counts = Sham_day7_data_2_b, project = "Sham_day_7_2_b", min.cells = 3, min.features = 200)
Sham_day7_data_3 <- CreateSeuratObject(counts = Sham_day7_data_3, project = "Sham_day_7_3", min.cells = 3, min.features = 200)

Sham_day7_data_1_a
Sham_day7_data_1_b
Sham_day7_data_2_a
Sham_day7_data_2_b
Sham_day7_data_3

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

MI_day1_QC<-QC_function(MI_day1_data)
MI_day3_QC<-QC_function(MI_day3_data)
MI_day5_QC<-QC_function(MI_day5_data)
MI_day7_QC<-QC_function(MI_day7_data)
MI_day14_QC<-QC_function(MI_day14_data)
MI_day28_QC<-QC_function(MI_day28_data)
Sham_day7_data_1_a_QC<-QC_function(Sham_day7_data_1_a)
Sham_day7_data_1_b_QC<-QC_function(Sham_day7_data_1_b)
Sham_day7_data_2_a_QC<-QC_function(Sham_day7_data_2_a)
Sham_day7_data_2_b_QC<-QC_function(Sham_day7_data_2_b)
Sham_day7_data_3_QC<-QC_function(Sham_day7_data_3)


MI_day1_QC
MI_day3_QC
MI_day5_QC
MI_day7_QC
MI_day14_QC
MI_day28_QC
Sham_day7_data_1_a_QC
Sham_day7_data_1_b_QC
Sham_day7_data_2_a_QC
Sham_day7_data_2_b_QC
Sham_day7_data_3_QC

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

MI_day1_doublet<-DF_Function(MI_day1_QC)
MI_d1_Object<-subset(MI_day1_doublet, DF.classifications_0.25_0.03_363=='Singlet')
MI_d1_Object

MI_day3_doublet<-DF_Function(MI_day3_QC)
MI_d3_Object<-subset(MI_day3_doublet, DF.classifications_0.25_0.06_414=='Singlet')
MI_d3_Object

MI_day5_doublet<-DF_Function(MI_day5_QC)
MI_d5_Object<-subset(MI_day5_doublet, DF.classifications_0.25_0.04_532=='Singlet')
MI_d5_Object

MI_day7_doublet<-DF_Function(MI_day7_QC)
MI_d7_Object<-subset(MI_day7_doublet, DF.classifications_0.25_0.13_469=='Singlet')
MI_d7_Object

MI_day14_doublet<-DF_Function(MI_day14_QC)
MI_d14_Object<-subset(MI_day14_doublet, DF.classifications_0.25_0.005_495=='Singlet')
MI_d14_Object

MI_day28_doublet<-DF_Function(MI_day28_QC)
MI_d28_Object<-subset(MI_day28_doublet, DF.classifications_0.25_0.02_403=='Singlet')
MI_d28_Object

Sham_day7_data_1_a_doublet<-DF_Function(Sham_day7_data_1_a_QC)
Sham_day7_data_1_a_Object<-subset(Sham_day7_data_1_a_doublet, DF.classifications_0.25_0.01_137=='Singlet')
Sham_day7_data_1_a_Object

Sham_day7_data_1_b_doublet<-DF_Function(Sham_day7_data_1_b_QC)
Sham_day7_data_1_b_Object<-subset(Sham_day7_data_1_b_doublet, DF.classifications_0.25_0.01_138=='Singlet')
Sham_day7_data_1_b_Object

Sham_day7_data_2_a_doublet<-DF_Function(Sham_day7_data_2_a_QC)
Sham_day7_data_2_a_Object<-subset(Sham_day7_data_2_a_doublet, DF.classifications_0.25_0.04_230=='Singlet')
Sham_day7_data_2_a_Object

Sham_day7_data_2_b_doublet<-DF_Function(Sham_day7_data_2_b_QC)
Sham_day7_data_2_b_Object<-subset(Sham_day7_data_2_b_doublet, DF.classifications_0.25_0.01_234=='Singlet')
Sham_day7_data_2_b_Object

Sham_day7_data_3_doublet<-DF_Function(Sham_day7_data_3_QC)
Sham_day7_data_3_Object<-subset(Sham_day7_data_3_doublet, DF.classifications_0.25_0.15_343=='Singlet')
Sham_day7_data_3_Object

Forte_Combined_Object=merge(Control_Object, y=c(MI_d1_Object, MI_d3_Object, MI_d5_Object, MI_d7_Object, MI_d14_Object, MI_d28_Object, Sham_day7_data_1_a_Object, Sham_day7_data_1_b_Object, Sham_day7_data_2_a_Object, Sham_day7_data_2_b_Object, Sham_day7_data_3_Object  ), add.cell.ids=c("Control", "Day_1", "Day_3", "Day_5", "Day_7", "Day_14","Day
                                                                                                                                           _28", "Sham_day7_1_a", "Sham_day7_1_b", "Sham_day7_2_a", "Sham_day7_2_b", "Sham_day7_3"), merge.data = TRUE, project="Forte")
Forte_Combined_Object

options(repr.plot.height = 5, repr.plot.width = 15)
VlnPlot(Forte_Combined_Object, features = c("nFeature_RNA", "nCount_RNA", "mito.percent"), group.by="orig.ident")

saveRDS(Forte_Combined_Object, 'Combined_Object.rds')
