rm(list = ls())
##########################################################################################
####################merge M2T3 and M1T4
##########################################################################################
library(reticulate)
reticulate::py_install(packages = 'umap-learn')
reticulate::py_install(packages= "numpy")

library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)

setwd('F:/Ko/67ESCC_PCN_integration_07-19-22/immune_cell_T4_M4')

M2T3 <- readRDS('M2T3.celltype_detail_after_cellchat.rds') 
M1T4 <- readRDS('M1T4.celltype_detail_after_cellchat.rds')



levels(M2T3@idents)
levels(M1T4@idents)


#Lift up CellChat object and merge together
#Since there are additional two populations (i.e., dermal DC and pericytes) specific to E14.5 
#compared to E13.5, we lift up cellchat.E13 by lifting up the cell groups to the same cell 
#labels as E14.5. liftCellChat will only update the slot related to the cell-cell 
#communication network, including slots object@net, object@netP and object@idents.

# Define the cell labels to lift up
group.new = levels(M2T3@idents)
M2T3 <- liftCellChat(M2T3, group.new)
M1T4 <- liftCellChat(M1T4, group.new)
GroupAM1_T3 <- liftCellChat(GroupAM1_T3, group.new)

saveRDS(M2T3, "M2T3_for_merge.rds")
saveRDS(M1T4, "M1T4_for_merge.rds")
saveRDS(GroupAM1_T3, "GroupAM1_T3_for_merge.rds")


#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...

object.list_12 <- list(M2T3 = M2T3, M1T4 = M1T4)
cellchat_12 <- mergeCellChat(object.list_12, add.names = names(object.list_12), cell.prefix = TRUE)
cellchat_12
saveRDS(cellchat_12, 'cellchat.M2T3_M1T4_merged.rds')



###############################################################
#Compare the total number of interactions and interaction strength
#To answer on question on whether the cell-cell communication is enhanced or not, 
#CellChat compares the the total number of interactions and interaction strength of 
#the inferred cell-cell communication networks from different biological conditions.
gg1 <- compareInteractions(cellchat_12, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_12, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#To identify the interaction between which cell populations showing significant changes, 
#CellChat compares the number of interactions and interaction strength among different cell populations.

#The differential number of interactions or interaction strength in the cell-cell communication network 
#between two datasets can be visualized using circle plot, where red (or blue) colored edges represent 
#increased (or decreased) signaling in the second dataset compared to the first one.

# red: increased in M1T4, blue: decreased in M1T4
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_12, weight.scale = T)
netVisual_diffInteraction(cellchat_12, weight.scale = T, measure = "weight")


#We can also show differential number of interactions or interaction strength in a greater details using 
#a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap 
#(incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). 
#In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared 
#to the first one.
######## The differential network analysis only works for pairwise datasets. 

# red: increased in M1T4, blue: decreased in M1T4
gg1 <- netVisual_heatmap(cellchat_12)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_12, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


#To better control the node size and edge weights of the inferred networks across different datasets, we compute 
#the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) 
#across all datasets.
weight.max <- getMaxWeight(object.list_12, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_12)) {
  netVisual_circle(object.list_12[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list_12)[i]))
}

#Identify signaling groups based on their functional similarit
library(umap)
cellchat_12 <- computeNetSimilarityPairwise(cellchat_12, type = "functional")
cellchat_12 <- netEmbedding(cellchat_12, type = "functional")
cellchat_12 <- netClustering(cellchat_12, type = "functional")
netVisual_embeddingPairwise(cellchat_12, type = "functional", label.size = 3.5)
#> Compute the distance of signaling networks between datasets 1 2
gg1 <- rankNet(cellchat_12, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_12, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


table(cellchat_12@idents$M2T3)
table(cellchat_12@idents$M1T4)


#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
group.new
netVisual_bubble(cellchat_12, sources.use = 11, targets.use = c(5:7,9),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat_12, sources.use = 11, targets.use = c(2,3,10,11,13), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat_12, sources.use = c(5,6,7,9), targets.use = c(2,3,10,12), comparison = c(1, 2), angle.x = 45)

#> Comparing communications on a merged object

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "M1T4"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat_12_1 <- identifyOverExpressedGenes(cellchat_12, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_12_1, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_12_1, net = net, datasets = "M1T4",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_12_1, net = net, datasets = "M2T3",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_12_1)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_12_1)

#We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram

#   [1] "B cell"           "CD4 T cell"       "CD8 T cell"       "Fibroblast"       "M1 Macrophage"    "M2 Macrophage"    "Macrophage"       "Mast cell"       
#   [9] "Monocyte"         "T cell"           "effector T cell"  "epithelial"       "exhausted T cell"

#c(12), targets.use = c(5:7,9)
# upregulated signaling
pairLR.use.up = net.up[, "interaction_name", drop = F]
# epi to Myeloid
gg1 <- netVisual_bubble(cellchat_12_1, pairLR.use = pairLR.use.up, sources.use = c(11), targets.use = c(5:7,9), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("upregulated in ", names(object.list_12)[1]))
gg1
# T to Myeloid  
gg1 <- netVisual_bubble(cellchat_12_1, pairLR.use = pairLR.use.up, sources.use = c(2,3,10,12), targets.use = c(5:7,9,11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("upregulated in ", names(object.list_12)[1]))
gg1
# Myeloid to T
gg1 <- netVisual_bubble(cellchat_12_1, pairLR.use = pairLR.use.up, sources.use = c(5:7,9), targets.use = c(2,3,10,11,12), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("upregulated in ", names(object.list_12)[1]))
gg1




























