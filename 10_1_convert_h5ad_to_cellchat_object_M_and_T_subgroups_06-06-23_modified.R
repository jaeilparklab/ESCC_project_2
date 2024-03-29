install.packages("tidyverse")
install.packages("ggplot2")

library(reticulate)
library(Seurat)
library(dplyr)
library(CellChat)
conda_list()
use_condaenv('py38')

#reticulate::py_discover_config()
#install.packages('anndata')
#anndata::install_anndata()
library(anndata)




setwd('F:/Ko/67ESCC_PCN_integration_07-19-22/immune_cell_T4_M4')

##############################################################################################################
##################  This is successful only if follow in order ###############################################
##############################################################################################################



### 1. read adata
ad <- import("anndata", convert = FALSE)


########################################################### 2nd try with overlaid dataset 09-08-2022
 
ad_object_1 <- ad$read_h5ad("M1toT1_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_1$X))
rownames(data.input) <- rownames(py_to_r(ad_object_1$var))
colnames(data.input) <- rownames(py_to_r(ad_object_1$obs))
# access meta data
meta.data <- py_to_r(ad_object_1$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")



saveRDS(cellchat_2, 'M1toT1_cellchat_celltype_detail.rds')



##############################################################################################
ad_object_2 <- ad$read_h5ad("M1toT2_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_2$X))
rownames(data.input) <- rownames(py_to_r(ad_object_2$var))
colnames(data.input) <- rownames(py_to_r(ad_object_2$obs))
# access meta data
meta.data <- py_to_r(ad_object_2$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M1toT2_cellchat_celltype_detail.rds')
##############################################################################################
ad_object_3 <- ad$read_h5ad("M1toT4_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_3$X))
rownames(data.input) <- rownames(py_to_r(ad_object_3$var))
colnames(data.input) <- rownames(py_to_r(ad_object_3$obs))
# access meta data
meta.data <- py_to_r(ad_object_3$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")



saveRDS(cellchat_2, 'M1toT4_cellchat_celltype_detail.rds')
##############################################################################################
ad_object_4 <- ad$read_h5ad("M2toT1_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_4$X))
rownames(data.input) <- rownames(py_to_r(ad_object_4$var))
colnames(data.input) <- rownames(py_to_r(ad_object_4$obs))
# access meta data
meta.data <- py_to_r(ad_object_4$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M2toT1_cellchat_celltype_detail.rds')
##############################################################################################
ad_object_5 <- ad$read_h5ad("M2toT2_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_5$X))
rownames(data.input) <- rownames(py_to_r(ad_object_5$var))
colnames(data.input) <- rownames(py_to_r(ad_object_5$obs))
# access meta data
meta.data <- py_to_r(ad_object_5$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M2toT2_cellchat_celltype_detail.rds')

##############################################################################################
ad_object_5 <- ad$read_h5ad("M2toT3_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_5$X))
rownames(data.input) <- rownames(py_to_r(ad_object_5$var))
colnames(data.input) <- rownames(py_to_r(ad_object_5$obs))
# access meta data
meta.data <- py_to_r(ad_object_5$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M2toT3_cellchat_celltype_detail.rds')
##############################################################################################
ad_object_7 <- ad$read_h5ad("M3toT1_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_7$X))
rownames(data.input) <- rownames(py_to_r(ad_object_7$var))
colnames(data.input) <- rownames(py_to_r(ad_object_7$obs))
# access meta data
meta.data <- py_to_r(ad_object_7$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")

#saveRDS(cellchat, 'P063_cellchat_crude.rds')
#saveRDS(cellchat_1, 'P063_cellchat_celltype.rds')
saveRDS(cellchat_2, 'M3toT1_cellchat_celltype_detail.rds')
##############################################################################################
ad_object_8 <- ad$read_h5ad("M3toT2_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_8$X))
rownames(data.input) <- rownames(py_to_r(ad_object_8$var))
colnames(data.input) <- rownames(py_to_r(ad_object_8$obs))
# access meta data
meta.data <- py_to_r(ad_object_8$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M3toT2_cellchat_celltype_detail.rds')

##############################################################################################
ad_object_9 <- ad$read_h5ad("M3toT4_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_9$X))
rownames(data.input) <- rownames(py_to_r(ad_object_9$var))
colnames(data.input) <- rownames(py_to_r(ad_object_9$obs))
# access meta data
meta.data <- py_to_r(ad_object_9$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")

#saveRDS(cellchat, 'P866_cellchat_crude.rds')
#saveRDS(cellchat_1, 'P866_cellchat_celltype.rds')
saveRDS(cellchat_2, 'M3toT4_cellchat_celltype_detail.rds')


##############################################################################################
ad_object_10 <- ad$read_h5ad("M4toT1_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_7$X))
rownames(data.input) <- rownames(py_to_r(ad_object_7$var))
colnames(data.input) <- rownames(py_to_r(ad_object_7$obs))
# access meta data
meta.data <- py_to_r(ad_object_7$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M4toT1_cellchat_celltype_detail.rds')
##############################################################################################
ad_object_11 <- ad$read_h5ad("M4toT2_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_8$X))
rownames(data.input) <- rownames(py_to_r(ad_object_8$var))
colnames(data.input) <- rownames(py_to_r(ad_object_8$obs))
# access meta data
meta.data <- py_to_r(ad_object_8$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M4toT2_cellchat_celltype_detail.rds')

##############################################################################################
ad_object_12 <- ad$read_h5ad("M4toT4_overlaid_for_cellchat.h5ad")

# access normalized data matrix
data.input <- t(py_to_r(ad_object_9$X))
rownames(data.input) <- rownames(py_to_r(ad_object_9$var))
colnames(data.input) <- rownames(py_to_r(ad_object_9$obs))
# access meta data
meta.data <- py_to_r(ad_object_9$obs)
meta <- meta.data

library(Matrix)
str(data.input)
dgc.Matrix <- as(data.input, "CsparseMatrix")
str(dgc.Matrix)


#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "crude")
#cellchat_1 <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat_2 <- createCellChat(object = dgc.Matrix, meta = meta, group.by = "celltype_detail")


saveRDS(cellchat_2, 'M4toT4_cellchat_celltype_detail.rds')



