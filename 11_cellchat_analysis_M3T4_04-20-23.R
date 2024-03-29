rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ade4)
library(CellChat)

setwd('F:/Ko/67ESCC_PCN_integration_07-19-22/immune_cell_T4_M4')



###################################################################

###################################################################

M3T4.celltype_detail <- readRDS('M3toT4_cellchat_celltype_detail.rds')
levels(M3T4.celltype_detail@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(M3T4.celltype_detail@idents)) # number of cells in each cell group
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
M3T4.celltype_detail@DB <- CellChatDB.use

#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
M3T4.celltype_detail <- subsetData(M3T4.celltype_detail) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel
options(future.globals.maxSize= 943718400)

#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
M3T4.celltype_detail <- identifyOverExpressedGenes(M3T4.celltype_detail)
M3T4.celltype_detail <- identifyOverExpressedInteractions(M3T4.celltype_detail)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#M3T4.celltype_detail <- projectData(M3T4.celltype_detail, PPI.human)
#Part II: Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network

M3T4.celltype_detail <- computeCommunProb(M3T4.celltype_detail)

#M3T4.celltype_detail <- computeCommunProb(
#  M3T4.celltype_detail,
#  type = c("triMean", "truncatedMean", "median"),
#  trim = NULL,
#  LR.use = NULL,
#  raw.use = TRUE,
#  population.size = FALSE,
#  do.fast = TRUE,
#  nboot = 100,
#  seed.use = 1L,
#  Kh = 0.5,
#  n = 1
#)


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
M3T4.celltype_detail <- filterCommunication(M3T4.celltype_detail, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
table(M3T4.celltype_detail@meta$celltype_detail)


df.net <- subsetCommunication(M3T4.celltype_detail, sources.use = c(1,13,16,17,19), targets.use = c(2,3,4,5,6,7,8,9,10,11,12,14,15,18,20,21,22)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(M3T4.celltype_detail, signaling = c("CCL")) #gives the inferred cell-cell communications mediated by 
#Infer the cell-cell communication at a signaling pathway level
M3T4.celltype_detail <- computeCommunProbPathway(M3T4.celltype_detail)
#Calculate the aggregated cell-cell communication network
M3T4.celltype_detail <- aggregateNet(M3T4.celltype_detail)

groupSize <- as.numeric(table(M3T4.celltype_detail@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(M3T4.celltype_detail@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(M3T4.celltype_detail@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- M3T4.celltype_detail@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
M3T4.celltype_detail@netP$pathways
#   [1] "MHC-I"    "MIF"      "COLLAGEN" "MHC-II"   "GALECTIN" "CD99"     "FN1"      "MK"       "SPP1"     "LCK"      "ITGB2"    "LAMININ"  "CLEC"     "CD45"     "ICAM"     "VISFATIN" "CD86"    
#[18] "THBS"     "CXCL"     "NECTIN"   "TIGIT"    "CD22"     "SELPLG"   "ANNEXIN"  "ALCAM"    "CD6"      "TNF"      "SEMA7"    "CD80"     "VCAM"     "IL1"      "FASLG"    "CD70"     "TGFb"    
#[35] "THY1"     "TENASCIN" "IFN-II"   "BAFF"     "CD39"     "CSF"      "CD46"     "JAM"      "IL16"     "GAS"      "APRIL"    "CD226"    "ANGPTL"   "CCL"      "FLT3"     "CADM"     "MPZ"     
#[52] "NOTCH"    "PECAM1"          

pathways.show <- c("TIGIT") #"CXCL","PTN","MHC-I","NOTCH"
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "chord")

pathways.show <- c("NECTIN") #"CXCL","PTN","MHC-I","NOTCH"
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "chord")

pathways.show <- c("PD-L1") #"CXCL","PTN","MHC-I","NOTCH"
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "chord")

pathways.show <- c("CXCL") #"CXCL","PTN","MHC-I","NOTCH"
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "chord")

pathways.show <- c("CCL") #"CXCL","PTN","MHC-I","NOTCH"
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "chord")

pathways.show <- c("CD226") 
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(M3T4.celltype_detail, signaling = pathways.show, layout = "chord")



levels(M3T4.celltype_detail@idents) # show factor levels of the cell labels

#  [1] "B cell"           "CD4 T cell"       "CD8 T cell"       "M1 Macrophage"    "M2 Macrophage"    "Mast cell"        "Mixed"            "Monocyte"         "Plasma cell"     
# [10] "T cell"           "effector T cell"  "epithelial"       "exhausted T cell"

# show all the interactions between epithelial cell and T cell 
#netVisual_chord_gene(M3T4.celltype_detail, sources.use = c(12), targets.use = c(2,8,9,11), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cell and Macrophage 
#netVisual_chord_gene(M3T4.celltype_detail, sources.use = c(12), targets.use = c(5,6), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions between epithelial cells and [T cell and all immune]
#netVisual_chord_gene(M3T4.celltype_detail, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,13), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(M3T4.celltype_detail, sources.use = c(6), targets.use = c(11,13),  signaling = c("CDH1",'CDH'), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(M3T4.celltype_detail, sources.use = c(6), targets.use = c(11,13), signaling = c("CDH1","CDH"), remove.isolate = FALSE)

netVisual_chord_gene(M3T4.celltype_detail, signaling = c("TIGIT"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(M3T4.celltype_detail, signaling = c("NECTIN"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(M3T4.celltype_detail, signaling = c("CXCL"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(M3T4.celltype_detail, signaling = c("CCL"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_chord_gene(M3T4.celltype_detail, signaling = c("PD-L1"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


#netVisual_chord_gene(M3T4.celltype_detail, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11),  signaling = c("MIF"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
#netVisual_bubble(M3T4.celltype_detail, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11), signaling = c("MIF"), remove.isolate = FALSE)
#netVisual_chord_cell(M3T4.celltype_detail, sources.use = c(12), targets.use = c(1,2,3,4,5,6,7,8,9,11),  signaling = c("MIF"), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)


netVisual_chord_gene(M3T4.celltype_detail, sources.use = c(11), targets.use = c(2,3,5,6,7,9,10,11,12), lab.cex = 0.5,legend.pos.x = 15,legend.pos.y = 10)
netVisual_bubble(M3T4.celltype_detail, sources.use = c(11), targets.use = c(2,3,5,6,7,9,10,11,13), remove.isolate = FALSE)
netVisual_bubble(M3T4.celltype_detail, sources.use = c(5,6,7,9), targets.use = c(2,3,10,11,12), remove.isolate = FALSE)
netVisual_bubble(M3T4.celltype_detail, sources.use = c(2,3,10,12), targets.use = c(5,6,7,9,11), remove.isolate = FALSE)


#> Comparing communications on a single object

saveRDS(M3T4.celltype_detail, file = 'M3T4.celltype_detail_after_cellchat.rds')
#M3T4.celltype_detail <- readRDS('M3T4.celltype_detail_after_cellchat.rds')
