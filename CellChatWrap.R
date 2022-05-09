## Author: zhengyiyi
## Time: 202111228
## Functions: run cellchat



########------------用法----------------####

usage <- function( input = NULL){
  #############  一、WT组 构建cellchat rds    #######
  rm(list = ls())
  library(Seurat)
  library(CellChat)
  library(patchwork)
  options(stringsAsFactors = FALSE)
  setwd("/Users/zhengyiyi/Desktop/projects/Db/RerunResult/20211228_Cellchat/")
  
  path <- "/Users/zhengyiyi/Desktop/projects/Db/RerunResult/FinalRds/"
  
  seurat_object <- readRDS(paste0(path, "f_db_from_combined.rds"))
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- NormalizeData(seurat_object)
  
  # 提取出 WT 数据
  seurat_object <- subset(seurat_object, orig.ident == "Female WT")
  
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
  
  prefix <- "female_hypo_wt"
  f_wt_cellchat <- generateCellChatRds(data.input, meta, prefix)
  
  #############  二、KO组 构建cellchat rds    #######
  
  rm(list = ls())
  library(Seurat)
  library(CellChat)
  library(patchwork)
  options(stringsAsFactors = FALSE)
  setwd("/Users/zhengyiyi/Desktop/projects/Db/RerunResult/20211228_Cellchat/")
  
  path <- "/Users/zhengyiyi/Desktop/projects/Db/RerunResult/FinalRds/"
  
  seurat_object <- readRDS(paste0(path, "f_db_from_combined.rds"))
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object <- NormalizeData(seurat_object)
  
  # 提取出 KO 数据
  seurat_object <- subset(seurat_object, orig.ident == "Female db/db")
  
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
  
  prefix <- "female_hypo_db"
  f_db_cellchat <- generateCellChatRds(data.input, meta, prefix)
  
  #############  三、两组之间细胞通讯的比较   #########
  rm(list = ls())
  library(Seurat)
  library(CellChat)
  library(patchwork)
  options(stringsAsFactors = FALSE)
  setwd("/Users/zhengyiyi/Desktop/projects/Db/RerunResult/20211228_Cellchat/")
  source("/Users/zhengyiyi/Desktop/code/FunctionScripts/CellChatWrap.R")
  
  ## 1.load objects and merge togethor 
  
  cellchat.Ctrl <- readRDS("female_hypo_wt_cellchat.rds")
  cellchat.Treat <- readRDS("female_hypo_db_cellchat.rds")
  
  table(cellchat.Ctrl@meta$labels)
  table(cellchat.Treat@meta$labels)
  
  # Define the cell labels to lift up
  group.new = levels(cellchat.Ctrl@meta$labels)
  cellchat.Treat <- liftCellChat(cellchat.Treat, group.new)
  
  object.list <- list(WT = cellchat.Ctrl, KO = cellchat.Treat)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
  
  
  data.dir <- './female_db_wt_comparison'
  dir.create(data.dir)
  setwd(data.dir)
  
  ## 2. 分析
  prefix <- "female_hypo"
  ## Part I: Predict general principles of cell-cell communication
  object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]])
  object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]])
  cellchat <- generalCellchatCompare(prefix, object.list, cellchat)
  
  ## Part II: Identify the conserved and context-specific signaling pathways
  
  cellchat <- identifyContextSpecificPathway(prefix,object.list, cellchat)
  
  ## Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
  pos.dataset <- "KO"
  cellchat <- identifyDysfunctionalLRs(prefix, cellchat, object.list, pos.dataset)
  
  ## Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
  cellchat <- generatePathwayFig(prefix, cellchat, object.list)
  
}
########------------generateCellChatRds---------#####

generateCellChatRds <- function(data.input, meta, prefix, CellChatDB){
  
  library(CellChat)
  library(openxlsx)
  print("meta must be dataframe and have a column named with labels that it's content is celltype information")
  cat("Gene number: ", nrow(data.input))
  cat("Cell number: ", ncol(data.input))
  
  print(dim(meta))
  print(head(meta))
  ## 1.create cellchat object
  cellchat <- createCellChat(object = data.input, meta = meta, 
                             group.by = "labels")
  ## 2.set "labels" as default cell identity
  cellchat <- setIdent(cellchat, ident.use = "labels") 
  print(" # show factor levels of the cell labels")
  print(levels(cellchat@idents))
  
  print("# number of cells in each cell group")
  groupSize <- as.numeric(table(cellchat@idents)) 
  
  ## 3.设置配体受体相互作用库-Set the ligand-receptor interaction database
  
  # CellChatDB in mouse contains 2,021 validated molecular interactions, 
  # including 60% of secrete autocrine/paracrine signaling interactions, 
  # 21% of extracellular matrix (ECM)-receptor interactions 
  # and 19% of cell-cell contact interactions. 
  print("1. Set the ligand-receptor interaction database")
  #CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data

  # use a subset of CellChatDB for cell-cell communication analysis
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  
  # use all CellChatDB for cell-cell communication analysis
  # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  
  cellchat@DB <- CellChatDB.use
  
  ##  4.数据预处理
  # Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 4) # do parallel
  
  print("2. run identifyOverExpressedGenes and identifyOverExpressedInteractions function")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.mouse) # 返回结果在cellchat@data.project
  
  print("LRs")
  print(head(cellchat@LR))
  print(dim(cellchat@LR$LRsig)) # 1307 * 11
  
  ## 5. Compute the communication probability and infer cellular communication network
  
  print(table(cellchat@meta$labels))
  cellchat <- computeCommunProb(cellchat)
  
  print("Filter out the cell-cell communication if there are only few number of
  cells in certain cell groups")
  
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  print(table(cellchat@idents))
  # subsetCommunication函数可以提取出推断的细胞通讯的结果,返回的结果是data.frame
  # thresh默认是0.05 认为是显著的互作
  print("3. run computeCommunProbPathway function")
  df.net <- subsetCommunication(cellchat, slot.name = "net")

  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)

  # at the level of signaling pathways
  df.netP <- subsetCommunication(cellchat, slot.name = "netP")

  ## Calculate the aggregated cell-cell communication network 计算聚合的cell-cell通信网络
  print("4. run aggregateNet function")
  cellchat <- aggregateNet(cellchat) 
  
  res <- list(net = df.net, netP = df.netP, LRsig = cellchat@LR$LRsig)
  write.xlsx(res, paste0(prefix, "_cellchat_net_netP_LRsis.xlsx"), rowNames = T, 
                         colNames = T, overwrite = T)
  ### 3.save rds
  print("5. save cellchat rds")
  saveRDS(cellchat, file = paste0(prefix, "_cellchat.rds"))
  return(cellchat)
}

########------------generalCellchatCompare-----------#######
generalCellchatCompare <- function(prefix, object.list, cellchat){
  
  ### Folder generate ####
  print("Part I Analysis")
  folder_name="1.PartI_GeneralCompareFig"
  if (file.exists(folder_name)){
    print("PartI_GeneralCompareFig existed.")
  }else{
    dir.create(folder_name)
  } 
  
  ### 1) Compare the total number of interactions and interaction strength ####
  print("1. Compare the total number of interactions and interaction strength")
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  p1 <- gg1 + gg2
  filename1 <- paste(prefix, "F1_total_num_interaction&strengh.pdf", sep = "_")
  ggsave(p1, filename = filename1, width = 10, height = 6)
  
  # Copy files to 2.Cluster
  file.copy(filename1, folder_name ,overwrite = T)
  file.remove(filename1) 
  
  ### 2) Compare the number of interactions and interaction strength among different cell populations ####
  #  where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
  
  print("2. Compare the number of interactions and interaction strength among different cell populations")
  print(levels(cellchat@meta$datasets))
  cat("First dataset: ", levels(cellchat@meta$datasets)[1])
  cat("Second dataset: ", levels(cellchat@meta$datasets)[2])
  
  par(mfrow = c(1,2), xpd=TRUE)
  gg_net1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
  gg_net2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  par(no.readonly = T) # 关掉par函数
  
  filename2 <- paste(prefix, "F2_net_compare_interaction&strength.pdf", sep = "_")
  pdf(filename2, width = 15, height = 10)
  print(gg_net2)
  dev.off()
  
  # Copy files to 1.PartI_GeneralCompareFig
  file.copy(filename2, folder_name ,overwrite = T)
  file.remove(filename2) 
  
  ##  heatmap
  print("We can also show differential number of interactions or interaction strength in a greater details using a heatmap.
The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). 
The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, 
red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.")
  
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  plots <- gg1 + gg2
  
  filename3 <- paste(prefix, "F3_heatmap_compare_interaction&strength.pdf", sep = "_")
  pdf(filename3, width = 20, height = 10)
  print(plots)
  dev.off()
  
  # Copy files to 1.PartI_GeneralCompareFig
  file.copy(filename3, folder_name ,overwrite = T)
  file.remove(filename3)
  
  
  # To better control the node size and edge weights of the inferred networks across different datasets,
  # we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights)
  # across all datasets.
  ## 这样可以看这个细胞组跟哪个细胞作用最强
  filename4 <- paste(prefix, "F4_maximum_interactions.pdf", sep = "_")
  weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
  pdf(filename4, 20, 10)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                     edge.weight.max = weight.max[2], edge.width.max = 10, 
                     title.name = paste0("Number of interactions - ", names(object.list)[i]))
  }
  dev.off()
  par(no.readonly = T) # 关掉par函数
  
  # Copy files to 1.PartI_GeneralCompareFig
  file.copy(filename4, folder_name ,overwrite = T)
  file.remove(filename4)
  
  ### 3) Compare the major sources and targets in 2D space ####
  
  print("3. Compare the major sources and targets in 2D space ####")
  ## 行值加列值然后减去对角线上自己的值
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
    gg[[i]] <- netAnalysis_signalingRole_scatter(gg[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  
  filename5 <- paste(prefix, "F5_2d_space.pdf", sep = "_")
  pdf(filename5, 15, 8)
  print(patchwork::wrap_plots(plots = gg))
  dev.off()
  
  # Copy files to 1.PartI_GeneralCompareFig
  file.copy(filename5, folder_name ,overwrite = T)
  file.remove(filename5)
  
  object.list <- lapply(object.list, function(x) {netAnalysis_computeCentrality(x)})
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  return(cellchat)
}

#######------------identifyContextSpecificPathway---------------#######

identifyContextSpecificPathway <- function(prefix,object.list, cellchat){
print("Part II Analysis")
print("Identify the conserved and context-specific signaling pathways")

folder_name="2.PartII_ContextSpecificPathwaysFig"
if (file.exists(folder_name)){
  print("2.PartII_ContextSpecificPathwaysFig existed.")
}else{
  dir.create(folder_name)
} 



### 3, Compare the overall information flow of each signaling pathway
print("3) Identify and visualize the conserved and context-specific signaling pathways")
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

filename3 <- paste(prefix, "F3_overrall_information_flow_of_signaling_pathway.pdf", sep = "_")

pdf(filename3, 12, 20)
print(gg1 + gg2)
dev.off()

file.copy(filename3, folder_name ,overwrite = T)
file.remove(filename3)

### 4.Compare outgoing (or incoming) signaling associated with each cell population

print("4) Compare outgoing (or incoming) signaling associated with each cell population")

library(ComplexHeatmap)

object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]])
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]])

## outgoing 
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15)
filename4 <- paste(prefix, "F4_outgoing_signal_cell_populations.pdf", sep = "_")
pdf(filename4, 15, 15)
print(draw(ht1 + ht2, ht_gap = unit(0.5, "cm")))
dev.off()
file.copy(filename4, folder_name ,overwrite = T)
file.remove(filename4)

## incoming 
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15, color.heatmap = "GnBu")
filename5 <- paste(prefix, "F5_incoming_signal_cell_populations.pdf", sep = "_")
pdf(filename5, 15, 15)
print(draw(ht3 + ht4, ht_gap = unit(0.5, "cm")))
dev.off()
file.copy(filename5, folder_name ,overwrite = T)
file.remove(filename5)

## incoming+outgoing
ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 15, color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 15, color.heatmap = "OrRd")
filename6 <- paste(prefix, "F6_outgoing_incoming_signal_cell_populations.pdf", sep = "_")
pdf(filename6, 15, 15)
print(draw(ht5 + ht6, ht_gap = unit(0.5, "cm")))
dev.off()
file.copy(filename6, folder_name ,overwrite = T)
file.remove(filename6)

return(cellchat)
}

######-----------identifyDysfunctionalLRs-----------------------########

identifyDysfunctionalLRs <- function(prefix, cellchat, object.list, pos.dataset){
  print("Part III Analysis")
  
  folder_name="3.PartIII_UpDownLRsFigByCommunicationProbabities"
  if (file.exists(folder_name)){
    print("3.PartIII_UpDownLRsFigByCommunicationProbabities existed.")
  }else{
    dir.create(folder_name)
  } 
  
  print("Identify dysfunctional signaling by comparing the communication probabities")

  gg1 <- netVisual_bubble(cellchat,  comparison = c(1, 2), 
                          angle.x = 90, remove.isolate = T, max.dataset = 2, 
                          title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  gg2 <- netVisual_bubble(cellchat,
                          comparison = c(1, 2),  angle.x = 90, remove.isolate = T,  max.dataset = 1, 
                          title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
  #> Comparing communications on a merged object
  res <- list(signaling.TreatIncreased = gg1$data, signaling.TreatDecreased = gg2$data)
  filename1 <- paste(prefix, "dysfunctional_signaling_by_comparing_communication_probabities_sig_by_communication_probabities.xlsx", 
                     sep = "_")
  write.xlsx(res, filename1, rowNames = T, colNames = T, overwrite = T)
  file.copy(filename1, folder_name ,overwrite = T)
  file.remove(filename1)
  
  if ( min(table(cellchat@idents$joint)) < 30){
    len <- sum(table(cellchat@idents$joint)!=min(table(cellchat@idents$joint)))
  }else{
    len <- length(levels(cellchat@idents$joint))
  }
  
  for (i in 1:len){
    print(levels(cellchat@idents$joint)[i])
    gg1_tm <- netVisual_bubble(cellchat,  comparison = c(1, 2), max.dataset = 2,
                               sources.use = i, targets.use = 1:len,
                               angle.x = 90, remove.isolate = T,
                               title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    gg2_tm <- netVisual_bubble(cellchat, 
                               sources.use = i, targets.use = 1:len,
                               comparison = c(1, 2), max.dataset = 1, angle.x = 90, remove.isolate = T,
                               title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    filename2 <- paste(prefix, levels(cellchat@idents$joint)[i], "to_othercelltype_up_down-regulated signaling in Treat.pdf", sep = "_")
    ggsave(gg1_tm + gg2_tm, filename = filename2, width = 20, height = 10)
    file.copy(filename2, folder_name ,overwrite = T)
    file.remove(filename2)
  }
  
  # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
  # pos.dataset = "KO" # 外部输入
  # neg.dataset = "WT"
  
  folder_name="4.PartIII_UpDownLRsFigByDEAnalysis"
  if (file.exists(folder_name)){
    print("4.PartIII_UpDownLRsFigByDEAnalysis existed.")
  }else{
    dir.create(folder_name)
  } 
  
  print("First, we use identifyOverExpressedGenes to identify DEGs")
  print("please input pos.dataset name")
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  # thresh.pc: Threshold of the percent of cells expressed in one cluster
  # thresh.fc: Threshold of Log Fold Change
  # thresh.p: Threshold of p-values
  # res: object@var.features
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                         pos.dataset = pos.dataset, features.name = features.name, 
                                         only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged CellChat object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  print("Second, we use netMappingDEG to map the results of differential expression analysis onto
      the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest")
  
  net <- netMappingDEG(cellchat, features.name = features.name, thresh = 0.05)
  
  print("extract the ligand-receptor pairs with upregulated ligands in Treat")
  net.up.Treat <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = NULL)
  print("extract the ligand-receptor pairs with downregulated in Treat")
  net.down.Treat <- subsetCommunication(cellchat, net = net, datasets = pos.dataset, ligand.logFC = -0.1, receptor.logFC = NULL)
  
  gene.up <- extractGeneSubsetFromPair(net.up.Treat, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down.Treat, cellchat)
  
  print("Third, we use netVisual_bubble to obtain dysfunction ligand-receptor pairs")
  
  pairLR.use.up = net.up.Treat[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,  comparison = c(1, 2), 
                          angle.x = 90, remove.isolate = T, max.dataset = 2, 
                          title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  #> Comparing communications on a merged object
  pairLR.use.down = net.down.Treat[, "interaction_name", drop = F]
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                          comparison = c(1, 2),  angle.x = 90, remove.isolate = T,  max.dataset = 1, 
                          title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
  #> Comparing communications on a merged object
  res <- list(net.up.Treat = net.up.Treat, net.down.Treat = net.down.Treat, 
              signaling.TreatIncreased = gg1$data, signaling.TreatDecreased = gg2$data)
  filename1 <- paste(prefix, "dysfunctional_signaling_by_comparing_communication_probabities_sig.xlsx", sep = "_")
  write.xlsx(res, filename1, rowNames = T, colNames = T, overwrite = T)
  file.copy(filename1, folder_name ,overwrite = T)
  file.remove(filename1)
  
  if ( min(table(cellchat@idents$joint)) < 30){
    len <- sum(table(cellchat@idents$joint)!=min(table(cellchat@idents$joint)))
  }else{
    len <- length(levels(cellchat@idents$joint))
  }
  
  for (i in 1:len){
    print(levels(cellchat@idents$joint)[i])
    gg1_tm <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,  comparison = c(1, 2), 
                               sources.use = i, targets.use = 1:len,
                               angle.x = 90, remove.isolate = T, max.dataset = 2, 
                               title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    gg2_tm <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                               sources.use = i, targets.use = 1:len,
                               comparison = c(1, 2),  angle.x = 90, remove.isolate = T,  max.dataset = 1, 
                               title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    filename2 <- paste(prefix, levels(cellchat@idents$joint)[i], "to_othercelltype_up_down-regulated signaling in Treat.pdf", sep = "_")
    ggsave(gg1_tm + gg2_tm, filename = filename2, width = 20, height = 10)
    file.copy(filename2, folder_name ,overwrite = T)
    file.remove(filename2)
  }
  return(cellchat)
}

#####-------------generatePathwayFig-----------########

generatePathwayFig <- function(prefix, cellchat, object.list ){
  print("Part IV Analysis")
  folder_name="5.PartIV_PathwayShowFig"
  if (file.exists(folder_name)){
    print("5.PartIV_PathwayShowFig existed.")
  }else{
    dir.create(folder_name)
  } 
  
  object.list[[1]]@netP$pathways
  object.list[[2]]@netP$pathways
  pathways.show <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
  
  if (file.exists(paste0(folder_name, "/CircleFig"))){
    print(paste0(folder_name, "/CircleFig"))
  }else{
    dir.create(paste0(folder_name, "/CircleFig"))
  } 
  
  for (j in pathways.show){
    
    print("pathway")
    print(j)
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = j) # control the edge weights across different datasets
    tmp_name <- paste(prefix, j, "signaling_pathway.pdf", sep = "_")
    pdf(tmp_name, 15, 10)
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = j, layout = "circle", 
                          edge.weight.max = weight.max[1], edge.width.max = 8,
                          signaling.name = paste(j, names(object.list)[i]))
    }
    dev.off()
    
    file.copy(tmp_name, paste(folder_name, "CircleFig", sep = "/"), overwrite = T)
    file.remove(tmp_name)
  }
  
  
  if (file.exists(paste0(folder_name, "/HeatmapFig"))){
    print(paste0(folder_name, "/HeatmapFig"))
  }else{
    dir.create(paste0(folder_name, "/HeatmapFig"))
  } 
  
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
    
    tmp_name <- paste(prefix, j, "heatmap_signaling_pathway_related_to_vascular.pdf", sep = "_")
    pdf(tmp_name , 15, 10)
    print(ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm")))
    dev.off()
    
    file.copy(tmp_name, paste(folder_name, "HeatmapFig", sep = "/"), overwrite = T)
    file.remove(tmp_name)
  }
  return(cellchat)
}


