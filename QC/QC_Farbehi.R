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

# Read data
MI_3GFP <- Read10X_h5('Cellranger_Rp_d3GFP/outs/CB/output_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_3TIP <- Read10X_h5('Cellranger_Rp_d3TIP/outs/CB/output_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_7GFP <- Read10X_h5('Cellranger_Rp_d7GFP/outs/CB/output_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_7TIP <- Read10X_h5('Cellranger_Rp_d7TIP/outs/CB/output_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_3GFP_Sham <- Read10X_h5('Cellranger_Sham_d3GFP/outs/output_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_7GFP_Sham <- Read10X_h5('Cellranger_Sham_d7GFP/outs/output_seurat.h5', use.names = TRUE, unique.features = TRUE)
MI_7TIP_Sham <- Read10X_h5('Cellranger_Sham_d7TIP/outs/output_seurat.h5', use.names = TRUE, unique.features = TRUE)

#Create Seurat Object
MI_3GFP_data<- CreateSeuratObject(counts = MI_3GFP, project = "MI_3GFP", min.cells = 3, min.features = 200)
MI_7GFP_data<- CreateSeuratObject(counts = MI_7GFP, project = "MI_7GFP", min.cells = 3, min.features = 200)
MI_3TIP_data<- CreateSeuratObject(counts = MI_3TIP, project = "MI_3TIP", min.cells = 3, min.features = 200)
MI_7TIP_data<- CreateSeuratObject(counts = MI_7TIP, project = "MI_7TIP", min.cells = 3, min.features = 200)
MI_3GFP_Sham_data<- CreateSeuratObject(counts = MI_3GFP_Sham, project = "MI_3GFP_Sham", min.cells = 3, min.features = 200)
MI_7GFP_Sham_data<- CreateSeuratObject(counts = MI_7GFP_Sham, project = "MI_7GFP_Sham", min.cells = 3, min.features = 200)
MI_7TIP_Sham_data<- CreateSeuratObject(counts = MI_7TIP_Sham, project = "MI_7TIP_Sham", min.cells = 3, min.features = 200)



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


MI_3GFP_QC<-QC_function(MI_3GFP_data)
MI_7GFP_QC<-QC_function(MI_7GFP_data)
MI_3TIP_QC<-QC_function(MI_3TIP_data)
MI_7TIP_QC<-QC_function(MI_7TIP_data)
MI_3GFP_Sham_QC<-QC_function(MI_3GFP_Sham_data)
MI_7GFP_Sham_QC<-QC_function(MI_7GFP_Sham_data)
MI_7TIP_Sham_QC<-QC_function(MI_7TIP_Sham_data)

MI_3GFP_QC
MI_7GFP_QC
MI_3TIP_QC
MI_7TIP_QC
MI_3GFP_Sham_QC
MI_7GFP_Sham_QC
MI_7TIP_Sham_QC




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


MI_3GFP_doublet<-DF_Function(MI_3GFP_QC)
MI_7GFP_doublet<-DF_Function(MI_7GFP_QC)
MI_3TIP_doublet<-DF_Function(MI_3TIP_QC)
MI_7TIP_doublet<-DF_Function(MI_7TIP_QC)
MI_3GFP_Sham_doublet<-DF_Function(MI_3GFP_Sham_QC)
MI_7GFP_Sham_doublet<-DF_Function(MI_7GFP_Sham_QC)
MI_7TIP_Sham_doublet<-DF_Function(MI_7TIP_Sham_QC)


head(MI_3GFP_doublet@meta.data)
DimPlot(MI_3GFP_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.13_261', pt.size=1)
table(MI_3GFP_doublet@meta.data$DF.classifications_0.25_0.13_261)
MI_3GFP_Object<- subset(MI_3GFP_doublet, DF.classifications_0.25_0.13_261 == 'Singlet')
MI_3GFP_Object

head(MI_7GFP_doublet@meta.data)

DimPlot(MI_7GFP_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.21_452', pt.size=1)

table(MI_7GFP_doublet@meta.data$DF.classifications_0.25_0.21_452)
MI_7GFP_Object<- subset(MI_7GFP_doublet, DF.classifications_0.25_0.21_452 == 'Singlet')
MI_7GFP_Object

head(MI_3TIP_doublet@meta.data)

DimPlot(MI_3TIP_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.06_487',, pt.size=1)

table(MI_3TIP_doublet@meta.data$DF.classifications_0.25_0.06_487)
MI_3TIP_Object<- subset(MI_3TIP_doublet, DF.classifications_0.25_0.06_487 == 'Singlet')
MI_3TIP_Object

head(MI_7TIP_doublet@meta.data)

DimPlot(MI_7TIP_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.03_417', pt.size=1)

table(MI_7TIP_doublet@meta.data$DF.classifications_0.25_0.03_417)
MI_7TIP_Object<- subset(MI_7TIP_doublet, DF.classifications_0.25_0.03_417 == 'Singlet')
MI_7TIP_Object

head(MI_3GFP_Sham_doublet@meta.data)

DimPlot(MI_3GFP_Sham_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.29_173', pt.size=1)

table(MI_3GFP_Sham_doublet@meta.data$DF.classifications_0.25_0.29_173)
MI_3GFP_Sham_Object<- subset(MI_3GFP_Sham_doublet, DF.classifications_0.25_0.29_173=='Singlet')
MI_3GFP_Sham_Object

head(MI_7GFP_Sham_doublet@meta.data)

DimPlot(MI_7GFP_Sham_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.16_421', pt.size=1)

table(MI_7GFP_Sham_doublet@meta.data$DF.classifications_0.25_0.16_421)
MI_7GFP_Sham_Object<- subset(MI_7GFP_Sham_doublet, DF.classifications_0.25_0.16_421=='Singlet')
MI_7GFP_Sham_Object

head(MI_7TIP_Sham_doublet@meta.data)

DimPlot(MI_7TIP_Sham_doublet, reduction = 'pca', group.by='DF.classifications_0.25_0.3_603', pt.size=1)

table(MI_7TIP_Sham_doublet@meta.data$DF.classifications_0.25_0.3_603)
MI_7TIP_Sham_Object<- subset(MI_7TIP_Sham_doublet, DF.classifications_0.25_0.3_603 =='Singlet')
MI_7TIP_Sham_Object

Combined_data=merge(MI_3GFP_Sham_Object, y=c(MI_7GFP_Sham_Object, MI_7TIP_Sham_Object, MI_3GFP_Object, MI_7GFP_Object, MI_3TIP_Object,MI_7TIP_Object),
                    add.cell.ids=c("MI_3GFP_Sham_Object", "MI_7GFP_Sham_Object", "MI_7TIP_Sham_Object", "MI_3GFP_Object", "MI_7GFP_Object", "MI_3TIP_Object","MI_7TIP_Object"), merge.data = TRUE, project="RP")

options(repr.plot.height = 5, repr.plot.width = 15)
VlnPlot(Combined_data, features = c("nFeature_RNA", "nCount_RNA", "mito.percent"), group.by="orig.ident")

saveRDS(Combined_data, 'RP_Combined_data.rds')
