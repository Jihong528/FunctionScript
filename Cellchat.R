### 20210923 Cellchat 尝试

library(Seurat)
library(CellChat)

setwd("/Users/zhengyiyi/Desktop/projects/Ob_Hypo/Ob_13wk_female/")
dir.create("20210923_CellChatTry")
setwd("20210923_CellChatTry")
seurat_object <- readRDS("/Users/zhengyiyi/Desktop/projects/Ob_Hypo/Ob_13wk_female/Final_Rds/20210301_female_Aggregated_seurat.sct.rds")

DefaultAssay(seurat_object)

## 1.提取出基因表达谱和构建细胞信息
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
head(data.input)
dim(data.input) # 19420 gene * 15449 cell
class(seurat_object@meta.data)
names(seurat_object@meta.data)

meta <- data.frame(group = seurat_object$orig.ident,
                   nCount_RNA = seurat_object$nCount_RNA,
                   nFeature_RNA = seurat_object$nFeature_RNA,
                   percent.mt = seurat_object$percent.mt,
                   percent.rp = seurat_object$percent.rp,
                   celltype_assign = seurat_object$celltype_assign,
                   row.names = rownames(seurat_object@meta.data)) # create a dataframe of the cell labels


## 2.先试下用ob组的数据 构建CellChat 对象
cell.use = rownames(meta)[meta$group == "ob/ob"] # extract the cell names from disease data
head(cell.use)
length(cell.use) # 6550 细胞

# Prepare input data for CelChat analysis
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

## 3.Set the ligand-receptor interaction database
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

## 4.数据预处理

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


## Inference of cell-cell communication network

## 5.Compute the communication probability and infer cellular communication network

table(cellchat@meta$celltype_assign)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# The cell-cell communication related with the following cell groups are excluded due to the few number of cells:  Fibroblast Macrophage NFO

# 配体受体
df.net <- subsetCommunication(cellchat, slot.name = "net")
write.csv(df.net, 'df.net_lr.csv')



# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# at the level of signaling pathways
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, 'df.netP.csv')
# Calculate the aggregated cell-cell communication network
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

# Due to the complicated cell-cell communication network, 
# we can examine the signaling sent from each cell group.
# Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.

mat <- cellchat@net$weight[c(1,2,3,4,6,7,9,11), c(1,2,3,4,6,7,9,11)]
groupSize <- groupSize[c(1,2,3,4,6,7,9,11)]
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize, 
                   weight.scale = T,
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

mat <- cellchat@net$weight[c(1,2,3,4,6,7,9,11), c(1,2,3,4,6,7,9,11)]
groupSize <- groupSize[c(1,2,3,4,6,7,9,11)]
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize, 
                   weight.scale = T,
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}


## Visualization of cell-cell communication network

# 查看都有哪些信号通路
opar <- par(no.readonly = TRUE)
par(opar)
cellchat@netP$pathways
levels(cellchat@idents)
pathways.show <- c("CXCL")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = c(1, 4, 6, 7, 9, 12) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    vertex.receiver = vertex.receiver, 
                    layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    vertex.receiver = vertex.receiver, 
                    layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    vertex.receiver = vertex.receiver, 
                    layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show,
                  color.heatmap = "Reds")

# Chord diagram
group.cellType <- c(rep("Astrocyte", 4), rep("Microglia", 4),
                    rep("Neuron", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,12) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show, 
                     pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver)

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver


vertex.receiver = seq(1, 3)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], 
            vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat,
                                 signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i],
                         "_L-R_contribution.pdf"), 
               plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# Bubble plot
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, 
                 targets.use = c(5:11),
                 signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object

pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF"))
netVisual_bubble(cellchat, sources.use = c(3,4), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)


# Part IV: Systems analysis of cell-cell communication network

# Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,
                                  width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat,
                                         signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))
