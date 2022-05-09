#########----------20211109 细胞通讯分析------------------------------------########

#############  一、WT组 WT placenta 构建cellchat object    #######
rm(list = ls())
library(Seurat)
library(CellChat)
setwd("/Users/zhengyiyi/Desktop/projects/Ob_Qpct_knockout_placenta/Analysis/20211021_CutoffFinalAnalysis/6.QpctPlacentaDEGsAnalysis/CellChat/")
seurat_object <- readRDS("/Users/zhengyiyi/Desktop/projects/Ob_Qpct_knockout_placenta/Analysis/FinalRds/20211022_Final_qpct_placenta.rds")
DefaultAssay(seurat_object)

## 提取出基因表达谱和构建细胞信息 
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
head(data.input)
dim(data.input) # 32285 gene * 9852 cell
class(seurat_object@meta.data)
names(seurat_object@meta.data)

meta <- data.frame(group = seurat_object$orig.ident,
                   nCount_RNA = seurat_object$nCount_RNA,
                   nFeature_RNA = seurat_object$nFeature_RNA,
                   percent.mt = seurat_object$percent.mt,
                   percent.rp = seurat_object$percent.rp,
                   celltype_assign = seurat_object$final_celltype,
                   big_celltype = seurat_object$big_celltype,
                   row.names = rownames(seurat_object@meta.data)) # create a dataframe of the cell labels

## 1.数据预处理部分， 构建cellchat object
cell.use = rownames(meta)[meta$group == "WT placenta"] # extract the cell names from disease data
head(cell.use)
length(cell.use) # 6567 细胞

# 1) Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
dim(data.input)
meta = meta[cell.use, ]
head(meta)
unique(meta$celltype_assign) # check the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, 
                           group.by = "celltype_assign")
rm(seurat_object)
head(cellchat@meta)
# 如果没把meta在构建的时候加进去的话可以用这个代码后续加进去
# cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels") 

# set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "celltype_assign") 
levels(cellchat@idents) # show factor levels of the cell labels

## 2)Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

cellchat@DB <- CellChatDB.use

# 3) 数据预处理
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)

### 第二大部分: 推断细胞互作网络
## Inference of cell-cell communication network
## 1) Compute the communication probability and infer cellular communication network

table(cellchat@meta$celltype_assign)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  Fibroblast Macrophage NFO
# 配体受体
df.net <- subsetCommunication(cellchat, slot.name = "net")
write.csv(df.net, 'df.net_lr.csv')


## 2) Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# at the level of signaling pathways
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, 'df.netP.csv')
## 3) Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# We can also visualize the aggregated cell-cell communication network.
# For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@meta$celltype_assign))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
### 3.save rds
saveRDS(cellchat, file = "cellchat_qpct_placenta_WT.rds")


#############  二、KO组 qct KO placenta 构建cellchat object  ########

rm(list = ls())
library(Seurat)
library(CellChat)
setwd("/Users/zhengyiyi/Desktop/projects/Ob_Qpct_knockout_placenta/Analysis/20211021_CutoffFinalAnalysis/6.QpctPlacentaDEGsAnalysis/CellChat/")
seurat_object <- readRDS("/Users/zhengyiyi/Desktop/projects/Ob_Qpct_knockout_placenta/Analysis/FinalRds/20211022_Final_qpct_placenta.rds")
DefaultAssay(seurat_object)

seurat_object@active.ident <- seurat_object$final_celltype
seurat_object <- subset(seurat_object, idents = c("17_SynT", "18_SpT", "20_GlyT", "9_Neutrophil_Ngp_high"),
                        invert = T)

#### 1.数据预处理部分， 构建cellchat object
# 1) 提取出基因表达谱和构建细胞信息 
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
head(data.input)
dim(data.input) # 32285 gene * 9852 cell
class(seurat_object@meta.data)
names(seurat_object@meta.data)

meta <- data.frame(group = seurat_object$orig.ident,
                   nCount_RNA = seurat_object$nCount_RNA,
                   nFeature_RNA = seurat_object$nFeature_RNA,
                   percent.mt = seurat_object$percent.mt,
                   percent.rp = seurat_object$percent.rp,
                   celltype_assign = seurat_object$final_celltype,
                   big_celltype = seurat_object$big_celltype,
                   row.names = rownames(seurat_object@meta.data)) # create a dataframe of the cell labels

table(meta$celltype_assign, meta$group)

cell.use = rownames(meta)[meta$group == "Qpct ko placenta"] # extract the cell names from disease data
head(cell.use)
length(cell.use) # 3282 细胞

# 2) Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
dim(data.input)
meta = meta[cell.use, ]
head(meta)
unique(meta$celltype_assign) # check the cell labels

levels_define <- c("1_Decidual_Acta2_high","2_Decidual_Htra3_high", "3_Endothelial", "4_T cell", "5_B cell", "6_NK cell",
                   "7_Neutrophil_Csf3r_high", "8_Basophil", "10_Monocyte_Cybb_high", "11_Monocyte_Chil3_high",
                   "12_Monocyte_Treml4_high", "13_Macrophage", "14_Dentritic", "15_Trophoblast_Ccn2_high", 
                   "16_Trophoblast_Fbln1_high", "19_TGC")
meta$celltype_assign <- factor(as.character(meta$celltype_assign),
                               levels = levels_define)

meta$big_celltype <- factor(as.character(meta$big_celltype),
                            levels = c("1_Decidual_Acta2_high", "2_Decidual_Htra3_high", "3_Endothelial", "4_T cell", "5_B cell", "6_NK cell",
                                       "8_Basophil", "Immune cell", "Neutrophil", "19_TGC", "Trophoblast"))

cellchat <- createCellChat(object = data.input, meta = meta, 
                           group.by = "celltype_assign")
rm(seurat_object)
head(cellchat@meta)
# 如果没把meta在构建的时候加进去的话可以用这个代码后续加进去
# cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels") 

# set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "celltype_assign") 
levels(cellchat@idents) # show factor levels of the cell labels

## 3)Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

cellchat@DB <- CellChatDB.use

# 4) 数据预处理
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)

### 第二大部分: 推断细胞互作网络
## Inference of cell-cell communication network
## 1) Compute the communication probability and infer cellular communication network

table(cellchat@meta$celltype_assign)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  Fibroblast Macrophage NFO
# 配体受体
df.net <- subsetCommunication(cellchat, slot.name = "net")
write.csv(df.net, 'Qpct_ko_df.net_lr.csv')


## 2) Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# at the level of signaling pathways
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, 'Qpct_ko_df.netP.csv')
## 3) Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# We can also visualize the aggregated cell-cell communication network.
# For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@meta$celltype_assign))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
### 3.save rds
saveRDS(cellchat, file = "cellchat_qpct_placenta_KO.rds")



#########  三、两组比较下   #######

#########   1.总体上比较了下    ########
## PART I merge
## 1) Load the required libraries
rm(list = ls())

library(CellChat)
library(patchwork)

setwd("/Users/zhengyiyi/Desktop/projects/Ob_Qpct_knockout_placenta/Analysis/20211021_CutoffFinalAnalysis/6.QpctPlacentaDEGsAnalysis/CellChat/")

## 2) Load CellChat object of each dataset and then merge together

cellchat.WT <- readRDS("cellchat_qpct_placenta_WT.rds")
cellchat.KO <- readRDS("cellchat_qpct_placenta_KO.rds")

# Define the cell labels to lift up
group.new = levels(cellchat.WT@meta$celltype_assign)

cellchat.KO <- liftCellChat(cellchat.KO, group.new)
#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots [email protected], [email protected], [email protected] in a single dataset...
object.list <- list(WT = cellchat.WT, KO = cellchat.KO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
#> Warning in mergeCellChat(object.list, add.names = names(object.list),
#> cell.prefix = TRUE): Prefix cell names!
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

cellchat

#An object of class CellChat created from a merged object with multiple datasets 
#1038 signaling genes.
#9841 cells.
data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)

## Part I: Predict general principles of cell-cell communication
# CellChat starts with the big picture to predict general principles of cell-cell communication. When comparing cell-cell communication among multiple biological conditions, it can answer the following biological questions:

# Whether the cell-cell communication is enhanced or not

# The interaction between which cell types is significantly changed

# How the major sources and targets change from one condition to another

## 1) Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1 <- gg1 + gg2

ggsave(p1,filename="F1_compare_total_number_interactions_and_interaction_strength.pdf",
       width = 10, height = 6)
## 2) Compare the number of interactions and interaction strength among different cell populations
#  where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
levels(cellchat@meta$datasets)
# [1] "WT placenta"      "Qpct ko placenta"

# net
par(mfrow = c(1,2), xpd=TRUE)
gg_net1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
gg_net2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
par(no.readonly = T) # 关掉par函数

pdf("F2_net_compare_interaction.pdf", width = 15, height = 10)
gg_net1
dev.off()

pdf("F2_net_compare_interaction_strength.pdf", width = 15, height = 10)
print(gg_net2)
dev.off()

# heatmap
# We can also show differential number of interactions or interaction strength in a greater details using a heatmap.
# The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). 
# The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, 
# red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

pdf("F3_heatmap_compare_interaction&interaction_strength_in_wt_and_ko.pdf", width = 20, height = 10)
gg1 + gg2
dev.off()


# To better control the node size and edge weights of the inferred networks across different datasets,
# we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights)
# across all datasets.
## 这样可以看这个细胞组跟哪个细胞作用最强

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("F4_maximum_interactions.pdf", 20, 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 1, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
par(no.readonly = T) # 关掉par函数
saveRDS(cellchat, "unaggreated_cellchat.rds")

## 3) Differential number of interactions or interaction strength among different cell types
# To simplify the complicated network and gain insights into the cell-cell communication at the cell type level, 
# we can aggregate the cell-cell communication based on the defined cell groups. Here we categorize the cell populations into three cell types,
# and then re-merge the list of CellChat object.
group.cellType <- c(rep("Decidual", 2), "Endothelial", 
                    rep("T_B_NK_Bas_Neu", 6), rep("Immune", 5), rep("Trophoblast", 2), rep("SynT_SpT_TGC_GlyT", 4))

group.cellType <- factor(group.cellType, levels = c("Decidual", "Endothelial", "T_B_NK_Bas_Neu", "Immune",
                                                    "Trophoblast", "SynT_SpT_TGC_GlyT"))

object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

# We then can show the number of interactions or interaction strength between any two cell types in each dataset.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
pdf("F5_maximum_interactions_aggregated.pdf", 20, 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# Simialrly, we can also show the differential number of interactions or interaction strength between any two cell types
# using circle plot. Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset 
# compared to the first one.

pdf("F6_maximum_interactions_strength_aggregated.pdf", 20, 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()
par(no.readonly = T) # 关掉par函数

## 4) Compare the major sources and targets in 2D space
# Comparing the outgoing and incoming interaction strength in 2D space
# allows ready identification of the cell populations with significant changes in sending
# or receiving signals between different datasets.

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})

weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(gg[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("F7_2d_space.pdf", 15, 8)
patchwork::wrap_plots(plots = gg)
dev.off()

object.list <- lapply(object.list, function(x) {netAnalysis_computeCentrality(x)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "3_Endothelial", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "13_Macrophage", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
pdf("F8_differential_outgoing_incoming_chanages.pdf", 15, 8)
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


## Part II: Identify the conserved and context-specific signaling pathways

# 1) Identify signaling networks with larger (or less) difference as well 
# as signaling groups based on their functional/structure similarity

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional", slot.name = "netP")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

pdf("F9_2D_visualization_functional_similarity_signaling_network.pdf", 15, 10)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

## 2) Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

pdf("F10_2D_visualization_structure_similarity_signaling_network.pdf", 15, 10)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2

## 3) Compute and visualize the pathway distance in the learned joint manifold
#> Compute the distance of signaling networks between datasets 1 2
pdf("F11_pathway_distance.pdf", 12, 8)
rankSimilarity(cellchat, type = "functional")
dev.off()


## 4) Identify and visualize the conserved and context-specific signaling pathways
# By comparing the information flow/interaction strengh of each signaling pathway,
# we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, 
# by change their information flow at one condition as compared to another condition.


######## 2. 这个图可以看到每个组富集到的信号通路图，感觉这个图也挺重要的 靶向了具体信号通路 ###########
# 4) Compare the overall information flow of each signaling pathway

# This bar graph can be plotted in a stacked mode or not.
# Significant signaling pathways were ranked based on differences in the overall information flow
# within the inferred networks between NL and LS skin. 
# The top signaling pathways colored red are enriched in NL skin,
# and these colored green were enriched in the LS skin.

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

pdf("F12_overrall_information_flow_of_signaling_pathway.pdf", 12, 8)
gg1 + gg2
dev.off()

################# 3.比较重要的部分，比较了下outgoing, incoming, overall信号 可以看到哪些细胞之间的信号通讯发生了显著变化##############
# 5) Compare outgoing (or incoming) signaling associated with each cell population
# The above analysis summarize the information from the outgoing and 
# incoming signaling together. We can also compare the outgoing (or incoming) 
# signaling pattern between two datasets, allowing to identify signaling pathways/ligand-receptors
# that exhibit different signaling patterns.
# We can combine all the identified signaling pathways from different datasets and thus
# compare them side by side, including outgoing signaling, incoming signaling and overall 
# signaling by aggregating outgoing and incoming signaling together.
# NB: rankNet also shows the comparison of overall signaling, 
# but it does not show the signaling strength in specific cell populations.

library(ComplexHeatmap)
#> Loading required package: grid

#> ComplexHeatmap version 2.7.1.1010
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pdf("F13_outgoing_signal_cell_populations.pdf", 15, 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pdf("F14_incoming_signal_cell_populations.pdf", 15, 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pdf("F15_outgoing_incoming_signal_cell_populations.pdf", 15, 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

############## 4， 靶向了差异的配体受体对 ##############
### Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs

# 1)  Identify dysfunctional signaling by comparing the communication probabities
# We can compare the communication probabilities mediated by ligand-receptor pairs
# from some cell groups to other cell groups. This can be done by setting comparison in the function netVisual_bubble.

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(10:20),  comparison = c(1, 2), angle.x = 45)


# Moreover, we can identify the upgulated (increased) and down-regulated (decreased) signaling ligand-receptor
# pairs in one dataset compared to the other dataset. This can be done by specifying max.dataset and min.dataset
# in the function netVisual_bubble. The increased signaling means these signaling have higher 
# communication probability (strength) in one dataset compared to the other dataset.


gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(9:13),  comparison = c(1, 2), signaling = "TGFb",
                        max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(9:13),  comparison = c(1, 2), signaling = "TGFb",
                        max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("F16_dysfunctional_signaling.pdf", 15, 10)
gg1 + gg2
dev.off()



#NB: The ligand-receptor pairs shown in the bubble plot can be accessed via signaling.LSIncreased = gg1$data.

# 2) Identify dysfunctional signaling by using differential expression analysis
# The above method for identifying the upgulated and down-regulated signaling
# is perfomed by comparing the communication probability between two datasets 
# for each L-R pair and each pair of cell groups. Alternative, we can identify 
# the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene 
# expression analysis. Specifically, we perform differential expression analysis between
# two biological conditions (i.e., NL and LS) for each cell group, and then obtain 
# the upgulated and down-regulated signaling based on the fold change of ligands 
# in the sender cells and receptors in the receiver cells. Such analysis can be done as follows.

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "KO"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "KO",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

write.csv(net.up, "胎盘中_ko_vs_wt_差异上调的配体受体对.csv")
write.csv(net.down, "胎盘中_ko_vs_wt_差异下调的配体受体对.csv")
write.csv(gene.up, "胎盘中_ko_vs_wt_差异上调的基因.csv")
write.csv(gene.down, "胎盘中_ko_vs_wt_差异下调的基因.csv")
#  We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

name <- levels(cellchat@meta$celltype_assign)
for (i in  1:20){
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = i, targets.use = c(1:20), comparison = c(1, 2), 
                          angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling", names(object.list)[2]))
  #> Comparing communications on a merged object
  
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = i, targets.use = c(1:20), comparison = c(1, 2), 
                          angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling", names(object.list)[2]))
  #> Comparing communications on a merged object
  
  pdf(paste0(name[i], "_","Interest_signaling_pathway_upregulated_LRs_.pdf"), 20, 10)
  print(gg1)
  dev.off()
  
  pdf(paste0(name[i], "_","Interest_signaling_pathway_downregulated_LRs.pdf"), 20, 10)
  print(gg2)
  dev.off()
}

# Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
pdf("F18_dysfunction_signaling_sig_using_Chord_diagram.pdf", 15, 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(11:20), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(10:20), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

## Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram

## 看看数据中有的通路名字
# CellChatDB <- CellChatDB.mouse
# colnames(CellChatDB$interaction)
object.list[[1]]@netP$pathways
object.list[[2]]@netP$pathways
# > object.list[[1]]@netP$pathways
# [1] "CCL"        "SPP1"       "IGF"        "MIF"        "GALECTIN"   "IL1"        "MK"         "COMPLEMENT" "CXCL"       "VEGF"       "PARs"       "CSF"       
# [13] "OSM"        "PTN"        "ANNEXIN"    "PDGF"       "TGFb"       "VISFATIN"   "KIT"        "GRN"        "PERIOSTIN"  "ANGPTL"     "BMP"        "IL6"       
# [25] "CHEMERIN"   "PROS"       "TNF"        "IL4"        "GAS"        "IFN-II"     "FASLG"      "WNT"        "TWEAK"      "APRIL"      "EDN"        "CALCR"     
# [37] "CD137"      "SEMA3"      "BRADYKININ" "EGF"        "IL2"        "ncWNT"     
# > object.list[[2]]@netP$pathways
# [1] "CCL"        "GALECTIN"   "MIF"        "MK"         "IGF"        "SPP1"       "IL1"        "CXCL"       "CSF"        "COMPLEMENT" "TGFb"       "VISFATIN"  
# [13] "PTN"        "OSM"        "VEGF"       "PARs"       "ANGPTL"     "ANGPT"      "GRN"        "TNF"        "KIT"        "ANNEXIN"    "PERIOSTIN"  "CHEMERIN"  
# [25] "PDGF"       "IL4"        "IL6"        "APRIL"      "BMP"        "PROS"       "IFN-II"     "CALCR"      "SEMA3"      "WNT"        "IL2"        "FGF"       
# [37] "IL16"       "ncWNT"      "GAS"        "FLT3"       "EDN"        "TWEAK"      "CX3C"       "CD137"      "GDF"        "NRG"        "FASLG"      "HH"        
# [49] "ACTIVIN"   


pathways.show <- c("TGFb", "BMP", 
                   "WNT", "VEGF", "CCL", "MIF",
                   "ANGPTL", "SEMA3") 
for (j in pathways.show){
  print("pathway")
  print(j)
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = j) # control the edge weights across different datasets
  pdf(paste0(j, "_", "signaling_pathway_related_to_vascular.pdf"), 15, 10)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = j, layout = "circle", 
                        edge.weight.max = weight.max[1], edge.width.max = 1,
                        signaling.name = paste(j, names(object.list)[i]))
  }
  dev.off()
}


pathways.show <- c("TGFb", "BMP", 
                   "WNT", "VEGF", "CCL", "MIF",
                   "ANGPTL", "SEMA3") 
for(j in pathways.show){
  par(mfrow = c(1,2), xpd=TRUE)
  ht <- list()
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = j, color.heatmap = "Reds",
                                 title.name = paste(j, "signaling ",names(object.list)[i]))
  }
  
  #> Do heatmap based on a single object 
  #> 
  #> Do heatmap based on a single object
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  pdf(paste0(j, "_", "heatmap_signaling_pathway_related_to_vascular.pdf"), 15, 10)
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()
  
}


# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}




############ 20211115 挑出我们感兴趣的四条通路中差异的配体受体对画图
rm(list = ls())
setwd("/Users/zhengyiyi/Desktop/projects/Ob_Qpct_knockout_placenta/Analysis/20211021_CutoffFinalAnalysis/")
setwd("6.QpctPlacentaDEGsAnalysis/CellChat/comparison/")

library(CellChat)
library(openxlsx)

cellchat <- readRDS("agg_cellchat.rds")
object.list <- readRDS("object.list.rds")
setwd("通路中差异的配体受体/")

list.files()

net.up <- read.xlsx("Summary_胎盘中_ko_vs_wt_差异上调的配体受体对.xlsx", 
                    sheet = 1, rowNames = T, colNames = T)
net.down <- read.xlsx("Summary_胎盘中_ko_vs_wt_差异上调的配体受体对.xlsx", 
                      sheet = 2, rowNames = T, colNames = T)

pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1,2,3,13), targets.use = c(1:20), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1,2,3,13), targets.use = c(1:20), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling", names(object.list)[2]))


pdf("Decidual_Endothelial_Macrophage_TGFb_MIF_SEMA3_upregulated_LRs_.pdf", 20, 10)
print(gg1)
dev.off()

pdf("Decidual_Endothelial_Macrophage_TGFb_MIF_SEMA3_downregulated_LRs.pdf", 20, 10)
print(gg2)
dev.off()

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = setdiff(1:20, c(1,2,3,13)), targets.use = c(1:20), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = setdiff(1:20, c(1,2,3,13)), targets.use = c(1:20), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling", names(object.list)[2]))


pdf("Other_celltypes_TGFb_MIF_SEMA3_upregulated_LRs.pdf", 30, 10)
print(gg1)
dev.off()

pdf("Other_celltypes_TGFb_MIF_SEMA3_downregulated_LRs.pdf", 30, 10)
print(gg2)
dev.off()
