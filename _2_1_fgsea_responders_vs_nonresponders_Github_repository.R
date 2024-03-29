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


setwd("F:/Ko/PD-1_responders/fgsea")
######################################################################################################
####################### fgsea using Wilcoxon gene rank list from Scanpy ##############################
######################################################################################################

###################### Responders to non-responders ################################################################

### KEGG ###

DEG <- read.csv("T_cell_response_fgesa.csv", sep = ',')
DEG_1 <- DEG %>% select(Responders_names, Responders_score, Responders_pvals)
DEG_1 %>% filter(Responders_pvals < 0.05)
DEG_1$gene <- DEG_1$Responders_names
row.names(DEG_1) <- DEG_1$gene

##### rank with log2FC may not correct ####
#### I used new rank (avg_log2FC * -log10(p.adj)) ####calculate and make new column in Excel ####
ranks_1 <- DEG_1$Responders_score
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

fwrite(fgseaResTidy, file="T_cell_Responders_GSEA_results_KEGG.csv", sep=",", sep2=c("", " ", ""))


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

fwrite(fgseaResTidy_2, file="T_cell_Responders_GSEA_results_REACTOME.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_2 = read.csv('T_cell_Responders_GSEA_results_REACTOME.csv', sep=",")


ggplot(fgseaResTidy_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="REACTOME Hallmark pathways") + 
  theme_minimal()


plotEnrichment(msigdbr_list_2[["REACTOME_PD_1_SIGNALING"]],
               ranks_1) + labs(title="REACTOME_PD_1_SIGNALING")
plotEnrichment(msigdbr_list_2[["REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA"]],
               ranks_1) + labs(title="REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA")
plotEnrichment(msigdbr_list_2[["REACTOME_INTERFERON_SIGNALING"]],
               ranks_1) + labs(title="REACTOME_INTERFERON_SIGNALING")
plotEnrichment(msigdbr_list_2[["REACTOME_INTERLEUKIN_1_SIGNALING"]],
               ranks_1) + labs(title="REACTOME_INTERLEUKIN_1_SIGNALING")
plotEnrichment(msigdbr_list_2[["REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING"]],
               ranks_1) + labs(title="REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING")
plotEnrichment(msigdbr_list_2[["REACTOME_SIGNALING_BY_INTERLEUKINS"]],
               ranks_1) + labs(title="REACTOME_SIGNALING_BY_INTERLEUKINS")
plotEnrichment(msigdbr_list_2[["REACTOME_TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"]],
               ranks_1) + labs(title="REACTOME_TGF_BETA_RECEPTOR_SIGNALING_IN_EMT_EPITHELIAL_TO_MESENCHYMAL_TRANSITION")




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

fwrite(fgseaResTidy_5, file="T_cell_Responders_GSEA_results_GOBP.csv", sep=",", sep2=c("", " ", ""))
fgseaResTidy_5 = read.csv("T_cell_Responders_GSEA_results_GOBP.csv", sep = ",")


ggplot(fgseaResTidy_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="F1T4vs.others GOBP pathways") + 
  theme_minimal()

plotEnrichment(msigdbr_list_5[["GOBP_RESPONSE_TO_INTERFERON_ALPHA"]],
               ranks_1) + labs(title="GOBP_Response_to_IFN_alpha_SIGNALING")
plotEnrichment(msigdbr_list_5[["GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION"]],
               ranks_1) + labs(title="GOBP_pos.regulation_of_Type_I_IFN_production")
plotEnrichment(msigdbr_list_5[["GOBP_INTERFERON_BETA_PRODUCTION"]],
               ranks_1) + labs(title="GOBP_INTERFERON_BETA_PRODUCTION")
plotEnrichment(msigdbr_list_5[["GOBP_TYPE_I_INTERFERON_PRODUCTION"]],
               ranks_1) + labs(title="GOBP_TYPE_I_INTERFERON_PRODUCTION")

plotEnrichment(msigdbr_list_5[["GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL"]],
               ranks_1) + labs(title="GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL")
plotEnrichment(msigdbr_list_5[["GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL"]],
               ranks_1) + labs(title="GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL")

