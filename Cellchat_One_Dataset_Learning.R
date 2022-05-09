#####------Part I: Data input & processing and initialization of CellChat object-----######
## 1. 提取出基因表达谱和构建细胞信息 
# 备注： cellchat需要两个信息:
# 1)Normalized data((e.g., library-size normalization and then log-transformed with a pseudocount of 1)),
# 2) cell group information, a dataframe with rownames is required as input for CellChat.

data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
head(data.input)
dim(data.input) # 32285 gene * 6628 cell
class(seurat_object@meta.data)
names(seurat_object@meta.data)

meta <- data.frame(group = seurat_object$orig.ident,
                   nCount_RNA = seurat_object$nCount_RNA,
                   nFeature_RNA = seurat_object$nFeature_RNA,
                   percent.mt = seurat_object$percent.mt,
                   percent.rp = seurat_object$percent.rp,
                   labels = seurat_object$celltype_assign,
                   row.names = rownames(seurat_object@meta.data)) # create a dataframe of the cell labels
unique(meta$labels) # check the cell labels

# [1] Astrocyte      Endothelial    Mature Oligo   Ependymal      Tanycyte       Immature Oligo Microglia     
# [8] Mural          Macrophage     Neuron         NFO            Fibroblast    
# 12 Levels: Astrocyte Endothelial Ependymal Tanycyte Fibroblast Immature Oligo Mature Oligo ... NFO

## 2. 构建cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, 
                           group.by = "labels")
rm(seurat_object)
head(cellchat@meta)
# 如果没把meta在构建的时候加进去的话可以用这个代码后续加进去
# cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels") 

# set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

str(cellchat)

## 3.设置配体受体相互作用库-Set the ligand-receptor interaction database

# CellChatDB in mouse contains 2,021 validated molecular interactions, 
# including 60% of secrete autocrine/paracrine signaling interactions, 
# 21% of extracellular matrix (ECM)-receptor interactions 
# and 19% of cell-cell contact interactions. 

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
str(CellChatDB)

write.xlsx(CellChatDB$interaction, "Mouse_CellChatDB_interaction.xlsx")
write.xlsx(table( CellChatDB$interaction$pathway_name, CellChatDB$interaction$annotation), "Mouse_CellchatDB_patwaygroup.xlsx")
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

cellchat@DB <- CellChatDB.use

##  4.数据预处理
# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
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
cellchat <- projectData(cellchat, PPI.mouse) # 返回结果在cellchat@data.project

dim(cellchat@data.signaling) # 1038 gene * 6628 cell
dim(cellchat@var.features$features.info) # 过表达的基因名字信息

head(cellchat@LR)
dim(cellchat@LR$LRsig) # 1307 * 11

colnames(cellchat@LR$LRsig)
# [1] "interaction_name"   "pathway_name"       "ligand"             "receptor"           "agonist"           
# [6] "antagonist"         "co_A_receptor"      "co_I_receptor"      "evidence"           "annotation"        
# [11] "interaction_name_2"
#####------Part II: Inference of cell-cell communication network 推断细胞互作网络----####

## 1. Compute the communication probability and infer cellular communication network

table(cellchat@meta$labels)
cellchat <- computeCommunProb(cellchat)

# computeCommunProb函数中 type	methods for computing the average gene expression per cell group.
# By default = "triMean", producing fewer but stronger interactions;
# 返回结果在@net中

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# 配体受体
# subsetCommunication函数可以提取出推断的细胞通讯的结果,返回的结果是data.frame
# thresh默认是0.05 认为是显著的互作
df.net <- subsetCommunication(cellchat, slot.name = "net")
write.csv(df.net, 'female.hypo.wt.df.net_lr.csv')


## 2 Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# at the level of signaling pathways
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, 'female.hypo.wt.df.netP.csv')

## 3 Calculate the aggregated cell-cell communication network 计算聚合的cell-cell通信网络
cellchat <- aggregateNet(cellchat) 
# 返回两个结果，count:两组细胞之间的相互作用数目； weight:两组细胞之间相互作用的加权数
head(cellchat@net$count)
head(cellchat@net$weight)

cellchat@netP$pathways # 显著的信号通路
head(cellchat@LR$LRsig)
dim(cellchat@LR$LRsig) # 显著的相互作用对

### 3.save rds
saveRDS(cellchat, file = "cellchat_female_hypo_WT.rds")


#####------Part III: Visualization of cell-cell communication network-----######
## 1.showing the number of interactions or the total interaction strength (weights) 看整体互作关系
# between any two cell groups using circle plot.

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


## 2.每个group发出的信号
# 作者通过参数参数 edge.weight.max控制了最大的边的权重，是mat中最大的数目，这样就可以比较不同网络之间的边权重。
mat <- cellchat@net$weight
group_length <- length(unique(cellchat@meta$labels))
row_number <- ceiling(group_length/4)
par(mfrow = c(row_number,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

par(no.readonly = T) # 关掉par函数
## 3.Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram 
## a.层次图
# 使用层次图、圆图或弦图可视化每个信号通路
cellchat@netP$pathways

pathways.show <- c("LAMININ") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout= "hierarchy")

## b. 圆图
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

## c.和弦图 Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
## d.热图
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

## 4.Compute the contribution of each ligand-receptor pair to 
# the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair

# 计算每个配体 - 受体对对整体信号通路的贡献，并可视化由单个配体 - 受体对介导的细胞 - 细胞通讯
netAnalysis_contribution(cellchat, signaling = pathways.show)

# 我们还可以可视化由单个配体 - 受体对介导的细胞 - 细胞通讯。 
# 我们提供了一个函数 extractEnrichedLR 来提取给定信号通路的所有重要相互作用（L-R 对）和相关信号基因。

pairLR.LAMININ <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.LAMININ[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver, layout = "hierarchy")

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Automatically save the plots of the all inferred network for quick exploratio
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways[1:5]
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

## 5.Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# 可视化由多个配体受体或信号通路介导的细胞间通讯
# a.Bubble plot

# a1. show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# a2. show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("PSAP","PTN"), 
                 remove.isolate = FALSE)


# a3. show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("PSAP","PTN","NCAM")) 
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)

# b.Chord diagram

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11),
                     lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by Macrophage
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

# c.Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "CXCL")

#By default, plotGeneExpression only shows the expression of signaling genes related to the inferred significant communications.
# USERS can show the expression of all signaling genes related to one signaling pathway by
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

#####------Part IV: Systems analysis of cell-cell communication network---#####

## 1.Identify signaling roles (e.g., dominant senders, receivers) of cell groups
## as well as the major contributing signaling

## a.Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap,
#allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 8, height = 2.5, font.size = 10)

## b.Visualize the dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

## c.Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

## 2.Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
library(NMF)
library(ggalluvial)

# a1.Identify and visualize outgoing communication pattern of target cells
selectK(cellchat, pattern = "outgoing")

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

# a2.Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "cellchat_female_hypo_WT.rds")
