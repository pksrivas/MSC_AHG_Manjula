
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

library(metap)

B3<-read.csv('Fishers_B - D3 (2).csv')

as.matrix(B3[1,])

BD3_Result<-apply(as.matrix(B3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

BD3_Result

write.csv(BD3_Result,'BD3_Result.csv')

B_D7<-read.csv('Fishers_B - D7 (1).csv', row.names = 1)

as.matrix(B_D7[1,])

BD7_Result<-apply(as.matrix(B_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(BD7_Result,'BD7_Result.csv')

B_D14<-read.csv('Fishers_B - D14.csv',row.names = 1)

as.matrix(B_D14[1,])

BD14_Result<-apply(as.matrix(B_D14),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(BD14_Result,'BD14_Result.csv')

G_D3<-read.csv('Granulocyte - 3 (2).csv',row.names = 1)

as.matrix(G_D3[1,])

G_D3_Result<-apply(as.matrix(G_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(G_D3_Result,'G_D3_Result.csv')



G_D7<-read.csv('Granulocyte - 7.csv',row.names = 1)

as.matrix(G_D7[1,])

G_D7_Result<-apply(as.matrix(G_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(G_D7_Result,'G_D7_Result.csv')



G_D14<-read.csv('Granulocyte - 14 (2).csv',row.names = 1)

G_D14

as.matrix(G_D14[1,])

G_D14_Result<-apply(as.matrix(G_D14),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(G_D14_Result,'G_D14_Result.csv')



DC_D3<-read.csv('DC_Common  - DC3 (1).csv',row.names = 1)

as.matrix(DC_D3[1,])

DC_D3

DC_D3_Result<-apply(as.matrix(DC_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(DC_D3_Result,'DC_D3_Result.csv')

DC_D7<-read.csv('DC_Common  - D7.csv',row.names = 1)

Macrophages_D3<-read.csv('Macrophages - D3.csv',row.names = 1)

as.matrix(Macrophages_D3[1,])

Macrophages_D3_Result<-apply(as.matrix(Macrophages_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Macrophages_D3_Result,'CombinedP/Macrophages_D3_Result.csv')

Macrophages_D7<-read.csv('Macrophages_FishersTest - D7 (1).csv',row.names = 1)

Macrophages_D7

as.matrix(Macrophages_D7[1,])

Macrophages_D7_Result<-apply(as.matrix(Macrophages_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Macrophages_D7_Result,'CombinedP/Macrophages_D7_Result.csv')



Macrophages_D14<-read.csv('Macrophages_FishersTest - D14 (1).csv',row.names = 1)

Macrophages_D14_Result<-apply(as.matrix(Macrophages_D14),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Macrophages_D14_Result,'Macrophages_D14_Result.csv')



Monocyte_D3<-read.csv('Monocyte - D3 (2).csv',row.names = 1)

Monocyte_D3

Monocyte_D3_Result<-apply(as.matrix(Monocyte_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Monocyte_D3_Result,'Monocyte_D3_Result.csv')

Monocyte_D7<-read.csv('Monocyte - D7 (1).csv',row.names = 1)

Monocyte_D3

Monocyte_D7_Result<-apply(as.matrix(Monocyte_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Monocyte_D7_Result,'Monocyte_D7_Result.csv')



Monocyte_D14<-read.csv('Monocyte - D14 (2).csv',row.names = 1)

Monocyte_D14_Result<-apply(as.matrix(Monocyte_D14),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Monocyte_D14_Result,'Monocyte_D14_Result.csv')



TNK_D3<-read.csv('TNK - D3 (2).csv',row.names = 1)

TNK_D3

TNK_D3_Result<-apply(as.matrix(TNK_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(TNK_D3_Result, 'CombinedP/TNK_D3_Result.csv')

TNK_D7<-read.csv('TNK - D7 (1).csv',row.names = 1)

TNK_D14<-read.csv('TNK - D14 (1).csv',row.names = 1)

TNK_D7_Result<-apply(as.matrix(TNK_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

TNK_D14_Result<-apply(as.matrix(TNK_D14),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(TNK_D7_Result, 'TNK_D7_Result.csv')

write.csv(TNK_D14_Result, 'TNK_D14_Result.csv')

Macro_Mono_D3

Macro_Mono_D3<-read.csv('Mono&Macro - D3 (3).csv',row.names = 1)

Macro_Mono_D3_Result<-apply(as.matrix(Macro_Mono_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Macro_Mono_D3_Result, 'Macro_Mono_D3_Result.csv')



Macro_Mono_D7<-read.csv('Mono&Macro - D7 (1).csv',row.names = 1)

Macro_Mono_D7_Result<-apply(as.matrix(Macro_Mono_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Macro_Mono_D7_Result, 'Macro_Mono_D7_Result.csv')



Mono_Granulocyte_D3<-read.csv('Mono&Granulocyte - D3 (1).csv',row.names = 1)

Mono_Granulocyte_D3_Result<-apply(as.matrix(Mono_Granulocyte_D3),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Mono_Granulocyte_D3_Result, 'Mono_Granulocyte_D3_Result.csv')



Mono_Granulocyte_D7<-read.csv('Mono&Granulocyte - D7 (1).csv',row.names = 1)

Mono_Granulocyte_D7_Result<-apply(as.matrix(Mono_Granulocyte_D7),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Mono_Granulocyte_D7_Result, 'Mono_Granulocyte_D7_Result.csv')



Mono_Granulocyte_D14<-read.csv('Mono&Granulocyte - D14 (1).csv',row.names = 1)

Mono_Granulocyte_D14_Result<-apply(as.matrix(Mono_Granulocyte_D14),1,function(x){sumlog(as.numeric(x[which(x!="#N/A")]))$p})

write.csv(Mono_Granulocyte_D14_Result, 'Mono_Granulocyte_D14_Result.csv')
