rm(list=ls())

install.packages('pak')
pak::pkg_install("r-lib/rlang")
devtools::install_github("igordot/msigdbr")
BiocManager::install("fgsea")

library(data.table)
library(fgsea)
library(ggplot2)
library(BiocParallel)
library(msigdbr)
library(xlsx)
library(dplyr)


setwd("F:/Ko/67ESCC_PCN_integration_07-19-22/immune_cell_T4_M4/fgsea")
######################################################################################################
####################### fgsea using Wilcoxon gene rank list from Scanpy ##############################
######################################################################################################

###################### B2T3 to others ################################################################

### KEGG ###

DEG <- read.csv("T_cell_BtoT_cluster_fgesa.csv", sep = ',')
DEG_1 <- DEG %>% select(B2T3_names, B2T3_score, B2T3_pvals)
DEG_1 %>% filter(B2T3_pvals < 0.05)
DEG_1$gene <- DEG_1$B2T3_names
row.names(DEG_1) <- DEG_1$gene

##### rank with log2FC may not correct ####
#### I used new rank (avg_log2FC * -log10(p.adj)) ####calculate and make new column in Excel ####
ranks_1 <- DEG_1$B2T3_score
#ranks_1 <- DEG_1$rank
names(ranks_1) <- row.names(DEG_1)
head(ranks_1)
fgsea.p1<-barplot(sort(ranks_1, decreasing = TRUE))

msigdbr_show_species()
genesets_1 = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
msigdbr_list = split(x = genesets_1$gene_symbol, f = genesets_1$gs_name)
fgseaRes_1 <- fgsea(msigdbr_list, ranks_1, minSize=15, maxSize = 500)
head(fgseaRes_1[order(padj, -abs(NES)), ], n=10)

topUp <- fgseaRes_1 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown <- fgseaRes_1 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
fgsea.p2<-plotGseaTable(msigdbr_list[topPathways$pathway], ranks_1, fgseaRes_1, gseaParam = 0.5)

fgseaResTidy <- fgseaRes_1 %>%  as_tibble() %>% arrange(desc(NES))

fwrite(fgseaResTidy, file="Tcell_B2T3_GSEA_results_KEGG.csv", sep=",", sep2=c("", " ", ""))


ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

#### REACTOME  ####

genesets_2 = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
msigdbr_list_2 = split(x = genesets_2$gene_symbol, f = genesets_2$gs_name)
msigdbr_list_2
fgseaRes_2 <- fgsea(pathways=msigdbr_list_2, ranks_1, minSize=5, maxSize = 500)
head(fgseaRes_2[order(padj, -abs(NES)), ], n=10)


topUp_2 <- fgseaRes_2 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown_2 <- fgseaRes_2 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways_2 <- bind_rows(topUp_2, topDown_2) %>% arrange(-ES)

## To remove previoud plots, use this code ## 
##dev.off() 
fgsea.p3<-plotGseaTable(msigdbr_list_2[topPathways_2$pathway], ranks_1, fgseaRes_2, gseaParam = 0.5)



fgseaResTidy_2 <- fgseaRes_2 %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy_2 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fwrite(fgseaResTidy_2, file="Tcell_B2T3_GSEA_results_REACTOME.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_2 = read.csv('Tcell_B2T3_GSEA_results_REACTOME.csv', sep=",")


ggplot(fgseaResTidy_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="REACTOME Hallmark pathways") + 
  theme_minimal()

plotEnrichment(msigdbr_list_2[["REACTOME_PD_1_SIGNALING"]],
               ranks_1) + labs(title="REACTOME_PD_1_SIGNALING")
plotEnrichment(msigdbr_list_2[["REACTOME_SIGNALING_BY_INTERLEUKINS"]],
               ranks_1) + labs(title="REACTOME_SIGNALING_BY_INTERLEUKINS")
plotEnrichment(msigdbr_list_2[["REACTOME_TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"]],
               ranks_1) + labs(title="REACTOME_TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION")
plotEnrichment(msigdbr_list_2[["REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT"]],
               ranks_1) + labs(title="REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT")
plotEnrichment(msigdbr_list_2[["REACTOME_FCGR_ACTIVATION"]],
               ranks_1) + labs(title="REACTOME_FCGR_ACTIVATION")



####### GOBP
genesets_5 = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
msigdbr_list_5 = split(x = genesets_5$gene_symbol, f = genesets_5$gs_name)
msigdbr_list_5
fgseaRes_5 <- fgsea(pathways=msigdbr_list_5, ranks_1, minSize=5, maxSize = 500)
head(fgseaRes_5[order(padj, -abs(NES)), ], n=10)


topUp_5 <- fgseaRes_5 %>% filter(ES > 0) %>% top_n(5, wt=-padj)
topDown_5 <- fgseaRes_5 %>% filter(ES < 0) %>% top_n(5, wt=-padj)
topPathways_5 <- bind_rows(topUp_5, topDown_5) %>% arrange(-ES)

## To remove previoud plots, use this code ## 
dev.off() 
fgsea.p5<-plotGseaTable(msigdbr_list_5[topPathways_5$pathway], ranks_DEG, fgseaRes_5, gseaParam = 0.5)



fgseaResTidy_5 <- fgseaRes_5 %>%  as_tibble() %>% arrange(desc(NES))
fgseaResTidy_5 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fwrite(fgseaResTidy_5, file="Tcell_B2T3_GSEA_results_GOBP.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_5 = read.csv("Tcell_B2T3_GSEA_results_GOBP.csv", sep = ",")

plotEnrichment(msigdbr_list_5[["GOBP_NEGATIVE_REGULATION_OF_LYMPHOCYTE_ACTIVATION"]],
               ranks_1) + labs(title="GOBP_Neg._Reg._Lymphocyte_Act._SIGNALING")

ggplot(fgseaResTidy_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="B2T3vs.others GOBP pathways") + 
  theme_minimal()
