library(clusterProfiler)

library(org.Mm.eg.db)

library(AnnotationDbi)

library(ggplot2)

getwd()

B_D3<- read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/BD3_Result.csv',row.names=1)

B_D3

genes_to_test<-rownames(B_D3)

GO_results_B3<-enrichGO(gene=genes_to_test, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")

GO_results_B3

fit3<-plot(dotplot(GO_results_B3, showCategory = 20))

options(repr.plot.height = 13, repr.plot.width = 13)
cnetplot(GO_results_B3)

B_D7<- read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/BD7_Result.csv',row.names=1)

B_D7_genes<-rownames(B_D7)

GO_results<-enrichGO(gene=B_D7_genes, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")

fit7<-plot(dotplot(GO_results, showCategory = 20))



B_D14<- read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/BD14_Result.csv',row.names=1)

Genes_B_D14<-rownames(B_D14)

Genes_B_D14_GO_results<-enrichGO(gene=Genes_B_D14, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")

options(repr.plot.height = 11, repr.plot.width = 8)
fit14<-plot(dotplot(Genes_B_D14_GO_results, showCategory = 20))



DC_D3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/DC_D3_Result.csv',row.names=1)

DC_D7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/DC_D7_Result.csv',row.names=1)



Genes_DC_D3<-rownames(DC_D3)

Genes_DC_D7<-rownames(DC_D7)



Genes_DC_D3_GO_results<-enrichGO(gene=Genes_DC_D3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")

Genes_DC_D7_GO_results<-enrichGO(gene=Genes_DC_D7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")



options(repr.plot.height = 8, repr.plot.width = 8)
fit<-plot(dotplot(Genes_DC_D3_GO_results, showCategory = 20))

options(repr.plot.height = 8, repr.plot.width = 8)
fit<-plot(dotplot(Genes_DC_D7_GO_results, showCategory = 20))

getwd()

Granulocyte_D3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/G_D3_Result.csv',row.names=1)

Genes_Granulocyte_D3<-rownames(Granulocyte_D3)
Genes_Granulocyte_D3GO_results<-enrichGO(gene=Genes_Granulocyte_D3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")

options(repr.plot.height = 9, repr.plot.width = 8)
fit<-plot(dotplot(Genes_Granulocyte_D3GO_results, showCategory = 20))

Granulocyte_D7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/G_D7_Result.csv',row.names=1)
Genes_Granulocyte_D7<-rownames(Granulocyte_D7)
Genes_Granulocyte_D7GO_results<-enrichGO(gene=Genes_Granulocyte_D7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
fit<-plot(dotplot(Genes_Granulocyte_D7GO_results, showCategory = 20))

Granulocyte_D14<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/G_D14_Result.csv',row.names=1)
Genes_Granulocyte_D14<-rownames(Granulocyte_D14)
Genes_Granulocyte_D14GO_results<-enrichGO(gene=Genes_Granulocyte_D14, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
fit<-plot(dotplot(Genes_Granulocyte_D14GO_results, showCategory = 20))

Macro_Mono_D3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Macro_Mono_D3_Result.csv',row.names=1)
Macro_Mono_D3<-rownames(Macro_Mono_D3)
Macro_Mono_D3GO_results<-enrichGO(gene=Macro_Mono_D3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 11, repr.plot.width = 8)
fit<-plot(dotplot(Macro_Mono_D3GO_results, showCategory = 20))

Macro_Mono_D7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Macro_Mono_D7_Result.csv',row.names=1)
Macro_Mono_D7<-rownames(Macro_Mono_D7)
Macro_Mono_D7GO_results<-enrichGO(gene=Macro_Mono_D7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 9, repr.plot.width = 8)
fit<-plot(dotplot(Macro_Mono_D7GO_results, showCategory = 20))

MacroD3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Macrophages_D3_Result.csv',row.names=1)
MacroD3<-rownames(MacroD3)
MacroD3_results<-enrichGO(gene=MacroD3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 9, repr.plot.width = 8)
fit<-plot(dotplot(MacroD3_results, showCategory = 20))

MacroD7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Macrophages_D3_Result.csv',row.names=1)
MacroD7<-rownames(MacroD7)
MacroD7_results<-enrichGO(gene=MacroD7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 9, repr.plot.width = 8)
fit<-plot(dotplot(MacroD7_results, showCategory = 20))

MacroD14<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Macrophages_D14_Result.csv',row.names=1)
MacroD14<-rownames(MacroD14)
MacroD14_results<-enrichGO(gene=MacroD14, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 9, repr.plot.width = 8)
fit<-plot(dotplot(MacroD14_results, showCategory = 20))

Mono_D3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Monocyte_D3_Result.csv',row.names=1)
Genes_Mono_D3<-rownames(Mono_D3)
Genes_Mono_D3GO_results<-enrichGO(gene=Genes_Mono_D3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 11, repr.plot.width = 8)
fit<-plot(dotplot(Genes_Mono_D3GO_results, showCategory = 20))

Mono_D7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Monocyte_D7_Result.csv',row.names=1)
Genes_Mono_D7<-rownames(Mono_D3)
Genes_Mono_D7GO_results<-enrichGO(gene=Genes_Mono_D7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 11, repr.plot.width = 8)
fit<-plot(dotplot(Genes_Mono_D7GO_results, showCategory = 20))

Mono_D14<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Monocyte_D14_Result.csv',row.names=1)
Genes_Mono_D14<-rownames(Mono_D14)
Genes_Mono_D14GO_results<-enrichGO(gene=Genes_Mono_D14, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 8, repr.plot.width = 8)
fit<-plot(dotplot(Genes_Mono_D14GO_results, showCategory = 20))

Mono_GranulocyteD3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Mono_Granulocyte_D3_Result.csv',row.names=1)

Genes_Mono_GranulocyteD3<-rownames(Mono_GranulocyteD3)
Genes_Mono_GranulocyteD3GO_results<-enrichGO(gene=Genes_Mono_GranulocyteD3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")

options(repr.plot.height = 11, repr.plot.width = 8)
fit<-plot(dotplot(Genes_Mono_GranulocyteD3GO_results, showCategory = 20))

Mono_GranulocyteD7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Mono_Granulocyte_D7_Result.csv',row.names=1)
Genes_Mono_GranulocyteD7<-rownames(Mono_GranulocyteD7)
Mono_GranulocyteD7GO_results<-enrichGO(gene=Genes_Mono_GranulocyteD7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 8, repr.plot.width = 8)
fit<-plot(dotplot(Mono_GranulocyteD7GO_results, showCategory = 20))

Mono_GranulocyteD14<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/Mono_Granulocyte_D14_Result.csv',row.names=1)
Genes_Mono_GranulocyteD14<-rownames(Mono_GranulocyteD14)
Mono_GranulocyteD14GO_results<-enrichGO(gene=Genes_Mono_GranulocyteD14, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 10, repr.plot.width = 8)
fit<-plot(dotplot(Mono_GranulocyteD14GO_results, showCategory = 20))

TNK_D3<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/TNK_D3_Result.csv',row.names=1)
Genes_TNK_D3<-rownames(TNK_D3)
Genes_TNK_D3GO_results<-enrichGO(gene=Genes_TNK_D3, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 10, repr.plot.width = 8)
fit<-plot(dotplot(Genes_TNK_D3GO_results, showCategory = 20))

TNK_D7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/TNK_D7_Result.csv',row.names=1)
Genes_TNK_D7<-rownames(TNK_D7)
Genes_TNK_D7GO_results<-enrichGO(gene=Genes_TNK_D7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 11, repr.plot.width = 8)
fit<-plot(dotplot(Genes_TNK_D7GO_results, showCategory = 20))

TNK_D14<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/TNK_D14_Result.csv',row.names=1)
Genes_TNK_D14<-rownames(TNK_D14)
Genes_TNK_D14GO_results<-enrichGO(gene=Genes_TNK_D14, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 8, repr.plot.width = 8)
fit<-plot(dotplot(Genes_TNK_D14GO_results, showCategory = 20))

DC_D7<-read.csv('/rds/general/user/mg2523/home/apps/OUTS/CombinedP/DC_D7_Result.csv',row.names=1)
Genes_DC_D7<-rownames(DC_D7)
Genes_DC_D7GO_results<-enrichGO(gene=Genes_DC_D7, OrgDb= "org.Mm.eg.db", keyType= "SYMBOL", ont="BP")
options(repr.plot.height = 11, repr.plot.width = 8)
fit<-plot(dotplot(Genes_DC_D7GO_results, showCategory = 20))

head(DC_D7)


