## 20201122 分析ob 6wk female的数据
## 导入rds 
rm(list=ls())
setwd("/Users/zhengyiyi/Desktop/projects/Ob/Ob_female_6wk")

source("/Users/zhengyiyi/Desktop/code/R/Seurat_function.R")
source("/Users/zhengyiyi/Desktop/code/R/diffprop_functions.R")
library(Seurat)
library(ggplot2)

## 读入rds
Aggregated_seurat <- readRDS("female_new_Aggregated_seurat.rds")

## 细胞类型鉴定
Aggregated_seurat$celltype_try <- Aggregated_seurat$celltype_assign
markers <- c("Agt","Flt1","Ccdc153","Rax","Col1a2","Pdgfra","Mrc1","Mobp","Cx3cr1","Vtn","Snap25","Syt1","Fyn")
FeaturePlot(Aggregated_seurat, features = c("Col1a2", "Vtn", "Fyn") , reduction = "tsne")
FeaturePlot(Aggregated_seurat, features = c("Cx3cr1", "Mrc1") , reduction = "tsne")
Pdgfra_label <- Aggregated_seurat@active.ident==8
Aggregated_seurat$celltype_try[Pdgfra_label] <- "Immature Oligo" 
fyn_label <- (Aggregated_seurat@active.ident==8) & (Aggregated_seurat@reductions$tsne@cell.embeddings[, 2] < -19)
Aggregated_seurat$celltype_try[fyn_label] <- "NFO"

## 画图
# 聚类图 未鉴定细胞类型
p1 <- DimPlot(Aggregated_seurat,reduction="tsne",
              label = TRUE,pt.size=0.8,
              label.size=4, repel=TRUE)
p1 <- p1 + labs(title="female resolution 0.5") +  theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(Aggregated_seurat,reduction="tsne",label = TRUE,
              pt.size=0.8,label.size=4, repel=TRUE, split.by="orig.ident")
pdf("female_Clustering_tsne_resolution_0.5_plot.pdf",7,7)
print(p1)
dev.off()
pdf("female_Clustering_tsne_resolution_0.5_plot_split_by_orig_ident.pdf",12,6)
print(p2)
dev.off()

# 鉴定完细胞类型的图片
p1 <- DimPlot(Aggregated_seurat,reduction="tsne",label = TRUE,
              pt.size=0.8,label.size=4, repel=TRUE,group.by="ident") + NoLegend()
p1 <- p1 + labs(title="Clustering resolution 0.5")
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5)) 

p2 <- DimPlot(Aggregated_seurat,reduction="tsne",label = TRUE,
              pt.size=1.2,label.size=4.5, repel=TRUE,group.by="celltype_try") + NoLegend()
p2  <- p2  + labs(title="Identification_celltype_by_Seurat_lable_transfer") +theme(plot.title = element_text(hjust = 0.5)) 

pdf("Figure1_Clustering_tsne_celltype_by_Seurat_resolution_0.5.pdf",15,6)
plot_grid(p1,p2,ncol=2,rel_widths = c(1.8,2),labels=c("A","B"),label_size = 20)
dev.off()

p2 <- DimPlot(Aggregated_seurat,reduction="tsne",label = TRUE,group.by="celltype_try",
              pt.size=0.8,label.size=4, repel=TRUE, split.by="orig.ident") + NoLegend()
pdf("Figure2_female_Clustering_tsne_celltype_split_by_orig_ident.pdf",12,6)
print(p2)
dev.off()


# marker 基因的表达
markers <- c("Agt","Flt1","Mia","Rax","Col1a2","Pdgfra","Mrc1","Mobp","Cx3cr1","Vtn","Snap25","Syt1","Fyn")

names <- c(levels(Aggregated_seurat$celltype_try)[1:10],"Neurons","Neurons", "NFO")
for (i in 1:13){
  p_plot <- paste0("p_",i)
  temp <- VlnPlot(object = Aggregated_seurat, features = markers[i], ncol = 1, group.by = "seurat_clusters", pt.size = 0) 
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # 
  temp <- temp +  NoLegend()
  temp <- temp + labs(title=paste0(markers[i],"_", names[i])) + theme(axis.title.y=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
  assign(p_plot,temp)
}
library(cowplot)
library(ggplot2)
plots_cluster <- plot_grid(p_1, p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11,p_12,p_13,ncol=2)

pdf("Figure3_violin_plot_marker_each_cluster.pdf",13,13)
plots_cluster
dev.off()


Aggregated_seurat$celltype <- Aggregated_seurat$celltype_try
markers <- c("Agt","Flt1","Ccdc153","Mia","Rax","Col1a2","Pdgfra","Mrc1","Mobp","Cx3cr1","Vtn","Snap25","Syt1","Fyn")
names <- c(levels(Aggregated_seurat$celltype_try)[1:3], "Ependymal",
           levels(Aggregated_seurat$celltype_try)[4:10],"Neurons","Neurons", "NFO")
for (i in 1:14){
  p_plot <- paste0("p_",i)
  temp <- VlnPlot(object = Aggregated_seurat, features = markers[i], ncol = 1, group.by = "celltype", pt.size = 0) 
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # 
  temp <- temp +  NoLegend() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
  temp <- temp + labs(title=paste0(markers[i],"_", names[i])) + theme(axis.title.y=element_blank()) + theme(plot.title = element_text(hjust = 0.5))
  assign(p_plot,temp)
}

library(cowplot)
library(ggplot2)
plots_celltype <- plot_grid(p_1, p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11,p_12,p_13,p_14,ncol=2)

pdf("Figure4_Violin_plot_marker_celltype.pdf",13,14)
plots_celltype
dev.off()

# marker 基因的featureplot图
markers <- c("Agt","Flt1","Ccdc153","Mia","Rax","Col1a2","Pdgfra","Mrc1","Mobp","Cx3cr1","Vtn","Snap25","Syt1","Fyn")
p <- FeaturePlot(Aggregated_seurat, features = markers,
                 reduction = "tsne", ncol = 4, col = c("lightgrey", "red"))
ggsave("Figure5_featureplot_plot_marker.pdf", p, width = 22, height = 16)

# 细胞数目检验
source("diffprop_functions.R")
# 确定每一个类别的Obesity以及Control之间的细胞数目以及比例并进统计学检验
sample_info <- as.factor(Aggregated_seurat$orig.ident)

cell_type <- Aggregated_seurat$celltype_try
cellnumber <- Sample_celltype_hist(sample_info,cell_type,"Condition_cell_type")

## Read in file of counts of cells in each population across conditions
obs.counts = as.matrix(read.csv("Condition_cell_type.csv", row.names = 1))
print(obs.counts)

res.table.ControlVsObesity = c()
## Go through a series of error probabilities
for (err_prob in c(0.05)) {
  tip.exp <- generateNull(obs.counts, n=100000, p=err_prob);
  ## Control vs Obesity
  res.1 = two.class.test(obs.counts, tip.exp, cond.control="Control", cond.treatment="Obesity",to.plot=F)
  res.table.ControlVsObesity = rbind(res.table.ControlVsObesity, res.1)
}
rownames(res.table.ControlVsObesity) = as.character(c(0.05))
write.csv(res.table.ControlVsObesity,"cell_type_test.csv")

#把检验的值标上去
cellnumber <- t(cellnumber)
cellnumber_pct <- data.frame(control_percent = round(cellnumber[,1]/sum(cellnumber[,1]),4) * 100,
                             obesity_percent = round(cellnumber[,2]/sum(cellnumber[,2]),4) * 100,
                             control_cellnumber = cellnumber[,1],
                             obesity_cellnumber = cellnumber[,2])
write.csv(cellnumber_pct,"female_cellnumber_pct.csv")

p_value_all <- as.matrix(read.csv("cell_type_test.csv",row.names=1,header=T))
p_value_errorp_0.05 <- p_value_all["0.05", ]
p_value <- as.numeric(round(p_value_errorp_0.05, 3)) ##四舍五入小数点后3位
p_value_mark <- c("ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns")

labels <- as.factor(c(rep("Control",nrow(cellnumber_pct)),rep("Obesity",nrow(cellnumber_pct))))
pct_number <- c(cellnumber_pct[,1],cellnumber_pct[,2])
celltypes <- c(rep(rownames(cellnumber_pct),2))
plot_data <- data.frame(pct = pct_number, label = labels, celltype = celltypes)

plot_data$label <- factor(plot_data$label,levels=c("Control","Obesity"))
plot_data$celltype <- factor(plot_data$celltype,levels=rownames(cellnumber))

# cellnumber_pct <- cellnumber_pct[sort(rownames(cellnumber_pct)),]
positions <- c()
for (i in 1:nrow(cellnumber_pct)){
  if(cellnumber_pct[i,1] > cellnumber_pct[i,2]){
    temp <- cellnumber_pct[i,1] + 2}
  else{
    temp <- cellnumber_pct[i,2] + 2}
  positions <- c(positions,temp)
}
x_min_value <- seq(0.75, 12.75, 1)[1:12] #要画的线的x轴的左边的值
x_max_value <- seq(1.25,13.25,1)[1:12] #要画的线的x轴的右边的值
p <- ggplot(plot_data, aes(x=celltype,y=pct_number,fill=labels)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() 
p <- p + geom_signif(annotations = p_value_mark, y_position = positions, xmin = x_min_value, xmax = x_max_value)
p1 <- p + theme(axis.text.x = element_text(size = 10, color = "black", face="bold",
                                           vjust = 0.7, hjust = 0.6, angle = 30))
p1 <- p1 + ggtitle("Cell population percentages across different conditions") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size = 20))
p1 <- p1 + xlab("celltype") + ylab("Percentage")
p1 <- p1 + theme(axis.title.x = element_text(size = 18)) # x轴字体大小
p1 <- p1 + theme(axis.title.y = element_text(size = 15)) # y轴字体大小
p1 <- p1 + scale_fill_manual(values=c("#377EB8","#E41A1C"))
pdf("Figure5_cellnumber_difference_between_control_obesity.pdf",10,6)
p1
dev.off()


## 找每一类细胞的ob 和control的差异基因以及差异的印记基因

##------秩和检验找差异基因

## 导入印迹基因
imprinted_genes <- read.csv("/Users/zhengyiyi/Desktop/ImprintedGene/imprinted_gene_list.csv",
                            header = T)
imprinted_genes$Gene_name <- as.character(imprinted_genes$Gene_name)

## 找所有细胞中Obesity 和 Control之间的差异基因

## 不卡logFC找差异基因 是为了得到全部的差异基因为了画热图
DEGs_Obesity_VS_Control <- FindMarkers(Aggregated_seurat,slot="data",
                                       logfc.threshold = 0, 
                                       ident.1 = 'Obesity', 
                                       ident.2 = 'Control',
                                       min.pct = 0.1,group.by = "orig.ident",test.use="wilcox")
write.csv(DEGs_Obesity_VS_Control, "Female_6wk_ob_vs_control_unfilter_DEGs.csv")


DEGs_Obesity_VS_Control <- FindMarkers(Aggregated_seurat,slot="data",
                                       logfc.threshold = 0.1,ident.1 = 'Obesity', 
                                       ident.2 = 'Control',min.pct = 0.25,group.by = "orig.ident",
                                       test.use="wilcox")
DEGs_Obesity_VS_Control <- DEGs_Obesity_VS_Control[DEGs_Obesity_VS_Control$p_val_adj <= 0.05, ]


#DEGs_Obesity_VS_Control <- DEGs_Obesity_VS_Control[DEGs_Obesity_VS_Control$p_val_adj <= 0.05, ]

## 找这些差异基因中是否有印迹基因
DEGs <- rownames(DEGs_Obesity_VS_Control)
imprinted_DEGs <- intersect(imprinted_genes$Gene_name,DEGs)
imprinted_DEGs_informations <- DEGs_Obesity_VS_Control[imprinted_DEGs,]

## 保存结果
results_DEGs <- list(all_cells = DEGs_Obesity_VS_Control)
results_imprinting_DEGs <- list(imprinting_DEGs_all_cells = imprinted_DEGs_informations)

## 找每一类Obesity 和 Control之间的差异基因
library(openxlsx)

#cell_type_old <- unique(Aggregated_seurat$celltype)

cell_type <- sort(unique(Aggregated_seurat$celltype_try))

save_name_1 <- c("all_cells",as.character(cell_type))
save_name_2 <- c("all_cells",as.character(cell_type))
for (i in 1:length(cell_type)){
  j  <- i + 1
  ## 找每一类Obesity 和 Control的差异基因
  subset_seurat <- subset(x = Aggregated_seurat, subset = celltype_assign == cell_type[i])
  DEGs_Obesity_VS_Control_cell_type <- FindMarkers(subset_seurat,slot="data",
                                                   logfc.threshold = 0.1,ident.1 = 'Obesity', 
                                                   ident.2 = 'Control',min.pct = 0.25,
                                                   group.by = "orig.ident",test.use="wilcox")
  DEGs_Obesity_VS_Control_cell_type <- DEGs_Obesity_VS_Control_cell_type[DEGs_Obesity_VS_Control_cell_type$p_val_adj <= 0.2, ]
  results_DEGs[[j]] <- DEGs_Obesity_VS_Control_cell_type
  
  ## 找这些差异基因中是否有印迹基因
  DEGs <- rownames(DEGs_Obesity_VS_Control_cell_type)
  imprinted_DEGs <- intersect(imprinted_genes$Gene_name,DEGs)
  imprinted_DEGs_informations <- DEGs_Obesity_VS_Control_cell_type[imprinted_DEGs,]
  results_imprinting_DEGs[[j]] <- imprinted_DEGs_informations
  
}

## 把list的列名重命名
names(results_DEGs) <-  save_name_1#重命名
names(results_imprinting_DEGs) <-  save_name_2#重命名

## 将结果写到excel中
write.xlsx(results_DEGs, file = "female_6wk_DEGs_Obesity_VS_Normal.xlsx", col.names = T, row.names = T)
write.xlsx(results_imprinting_DEGs, file = "female_6wk_imprinting_DEGs_Obesity.xlsx", col.names = T, row.names = T)

## 找出脑中印迹的基因,并把这些基因按细胞类型合并到excel中来。
brain_imprinted_genes <- read.csv("Imprinted_genes_of_brain.csv",header = T)

imprinted_DEGs_genes <- data.frame()
for (i in 1:length(results_imprinting_DEGs)){
  row_names <- rownames(results_imprinting_DEGs[[i]])
  if(length(row_names)== 0){
    next
  }
  else{
    temp = results_imprinting_DEGs[[i]]
    celltype_names <- names(results_imprinting_DEGs[i])
    temp$celltype <- rep(celltype_names,length(row_names))
    temp$genename <- row_names
    imprinted_DEGs_genes <- rbind(imprinted_DEGs_genes,temp)
  }
}
write.csv(imprinted_DEGs_genes,"female_6wk_imprinted_DEGs_genes.csv")

same_genenames <- intersect(as.character(brain_imprinted_genes$Gene_name), as.character(imprinted_DEGs_genes$genename))

brain_imprinted_DEGs_genes <- data.frame()
for (i in 1:length(same_genenames)){
  index <- imprinted_DEGs_genes$genename==same_genenames[i]
  temp <- imprinted_DEGs_genes[index, ]
  brain_imprinted_DEGs_genes <- rbind(brain_imprinted_DEGs_genes,temp)
}

o <- order(brain_imprinted_DEGs_genes[,"celltype"])
brain_imprinted_DEGs_genes_order <- brain_imprinted_DEGs_genes[o, ]
write.csv(brain_imprinted_DEGs_genes_order,"female_6wk_brain_imprinted_DEGs_genes_order.csv")

saveRDS(Aggregated_seurat, file = "female_new_Aggregated_seurat.rds")

## 所有细胞类型的差异印迹基因的热图
library(Seurat)
library(pheatmap)
library(RColorBrewer)

# features_1 <- c("Meg3", "Ndn", "Nap1l5", "Snrpn", "Gnas", "Peg3", "Impact", "Usp29", "Ube3a", "Nnat")
brain_imprinted_DEGs_genes_order <- read.csv("female_6wk_brain_imprinted_DEGs_genes_order.csv",
                                             header = T, row.names = 1)
features_1 <- brain_imprinted_DEGs_genes_order[1:8,7]

features_1 <- c("Nnat", "Nap1l5", "Ndn", "Peg3", "Snrpn", "Ube3a", "Gnas")
annotation_col <- data.frame(Aggregated_seurat@meta.data$orig.ident,
                             row.names=rownames(Aggregated_seurat@meta.data))
colnames(annotation_col) <- 'group'
heatmap_byGroup <- FetchData(object=Aggregated_seurat,vars=features_1,slot='scale.data')
heatmap_byGroup<-cbind(heatmap_byGroup,annotation_col)
heatmap_byGroup<-heatmap_byGroup[, 1:7]

pdf("female_heatmap_imprinted_genes.pdf",12,5)
pheatmap::pheatmap(t(heatmap_byGroup), #fontsize_row=3, 
                   color=colorRampPalette(c("#CD00CD","black","yellow"))(100),
                   breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=10, treeheight_col=2, 
                   border_color='grey',cluster_cols = F,cluster_rows = F,
                   fontsize_row = 18,
                   annotation_col = annotation_col,
                   show_colnames = F,scale='none',
                   border=TRUE, gaps_col = 5629,clustering_method = "average")
dev.off()

## 所有细胞中的差异印记基因的featureplot 图
features_1 <- c("Nnat", "Nap1l5", "Ndn", "Peg3", "Snrpn", "Ube3a", "Gnas")
FeaturePlot(Aggregated_seurat, features = features_1, 
            ncol = 4, cols = c("lightgrey", "red"),
            reduction = "tsne")
ggsave("Female_imprinted_DEGs_FeatuePlot.pdf",height = 8, width = 16)

## 各细胞类型的印迹基因热图
Aggregated_seurat <- readRDS("female_6wk_Aggregated_seurat.rds")
features_1 <- c("Nap1l5", "Ndn", "Peg3", "Snrpn", "Ube3a", "Gnas")
annotation_col <- data.frame(subtype = Aggregated_seurat$celltype,
                             row.names = rownames(Aggregated_seurat@meta.data))
library(dplyr)
library(pheatmap)
library(RColorBrewer)
heatmap_byGroup <- FetchData(object=Aggregated_seurat,
                             vars=features_1,
                             slot='scale.data')
heatmap_byGroup <- cbind(heatmap_byGroup,annotation_col)
heatmap_byGroup <- arrange(heatmap_byGroup, subtype)

heatmap_byGroup_order <- heatmap_byGroup[, 1:6]
colorCount=18
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
col1=getPalette(colorCount)[1:12]

names(col1) <- levels(annotation_col$subtype)

ann_colors = list(group = col1)

library(pheatmap)
pdf("female_6wk_Imprinted_DEGs_celltypes.pdf",11,5)
pheatmap(t(heatmap_byGroup_order), #fontsize_row=3, 
         colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
         breaks=seq(-1, 1, length.out = 10),
         treeheight_row=0.5, treeheight_col=1, 
         border_color='grey', cluster_cols = F,cluster_rows = F,
         fontsize_row = 18,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F,show_rownames = T, scale='none',
         border=TRUE, angle_col = "0",
         fontsize = 10)
dev.off()

## neurons 中的这几个印记基因的表达情况
Imprinted_genes <- features_1

annotation_col <- data.frame(subtype = HFD_neuron$celltype_assign_small,
                             row.names = rownames(HFD_neuron@meta.data))
library(dplyr)
heatmap_byGroup <- FetchData(object=HFD_neuron,vars=Imprinted_genes,slot='scale.data')
heatmap_byGroup <- cbind(heatmap_byGroup,annotation_col)
heatmap_byGroup <- arrange(heatmap_byGroup, subtype)
a <- levels(heatmap_byGroup$subtype)
levs <- c(a[1], a[3:10], a[2], a[11:15])
levels(heatmap_byGroup$subtype) <- levs
heatmap_byGroup_order <- heatmap_byGroup[, 1:(length(Imprinted_genes)-1)]

colorCount=18
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
col1=getPalette(colorCount)[1:15]

names(col1) <- levels(annotation_col$subtype)

ann_colors = list(group = col1)

library(pheatmap)
pdf("20201106_HFD_neurons_subtypes_Imprinted_DEGs.pdf",11,4)
pheatmap(t(heatmap_byGroup_order), #fontsize_row=3, 
         colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10),
         breaks=seq(-1, 1, length.out = 10),
         treeheight_row=0.5, treeheight_col=1, 
         border_color='grey', cluster_cols = F,cluster_rows = F,
         fontsize_row = 18,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F,show_rownames = T, scale='none',
         border=TRUE, angle_col = "0",
         fontsize = 10)
dev.off()


###########---------20201125 画差异基因数目的barplot图------------######
setwd("/Users/zhengyiyi/Desktop/projects/Ob_female_6wk/20201125_Compare_DEGs_with_13wk/")
name <- getSheetNames('female_6wk_DEGs_Obesity_VS_Normal.xlsx')
file = "female_6wk_DEGs_Obesity_VS_Normal.xlsx"
data <- gene_number_stats(file,length(name),name)  

## plot
library(ggplot2)
plot_data <- data.frame(number=c(data$up,-1*data$down),
                        celltype=rep(data$celltype,2),
                        ud=c(rep("up",dim(data)[1]),rep("down",dim(data)[1])))

plot_data$celltype <- factor(plot_data$celltype,levels=name)

range(plot_data$number)

min_number <- min(plot_data$number)-30
max_number <- max(plot_data$number)+30
p1 <- ggplot(plot_data, aes(x=celltype, y=number, fill = ud)) +
  geom_col(position = position_dodge(width = 0), width = 0.6, size = 0.3, colour = "black")+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = "black", 
                                        fill = "transparent"),
        legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "celltype", y = "number") +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_fill_manual(values=c("#377EB8","#E41A1C"))+
  scale_y_continuous(breaks = seq(min_number, max_number, 50),
                     labels = as.character(abs(seq(min_number, max_number, 50))),
                     limits = c(min_number, max_number))
ggsave("Ob_female_6wk_DEG_number_barplot.pdf", p1, width = 8, height = 4)


#########------20201125 比较6 周和13周的差异基因------------------######
sheetname_6wk <- getSheetNames("female_6wk_DEGs_Obesity_VS_Normal.xlsx")
sheetname_13wk <- getSheetNames("female_13wk_DEGs_Obesity_VS_Normal.xlsx")
DEGs_6wk <- read.xlsx("female_6wk_DEGs_Obesity_VS_Normal.xlsx",
                      colNames = T)
DEGs_13wk <- read.xlsx("female_13wk_DEGs_Obesity_VS_Normal.xlsx",
                       )

female_13wk <- readRDS("female_13wk_Aggregated_seurat.rds")
female_6wk <- readRDS("/Users/zhengyiyi/Desktop/projects/Ob_female_6wk/female_new_Aggregated_seurat.rds")

cellnumber_13wk <- table(female_13wk$orig.ident, female_13wk$celltype_assign)
cellnumber_6wk <- table(female_6wk$orig.ident, female_6wk$celltype_try)

cellnumber <- rbind(cellnumber_13wk, cellnumber_6wk)

library(openxlsx)
write.xlsx(cellnumber, "summary_cellnumber.xlsx", colNames = T, rowNames = T)

## 在6wk 13wk中 都上调/下调的差异基因
upgenelist <- list()
downgenelist <- list()
for (i in 1:length(sheetname_6wk)){
group1 <- read.xlsx("female_6wk_DEGs_Obesity_VS_Normal.xlsx", 
                      sheet = i, 
                      rowNames = T)
group2 <- read.xlsx("female_13wk_DEGs_Obesity_VS_Normal.xlsx", 
                      sheet = i, 
                      rowNames = T)
up_gene <- intersect(rownames(group1)[group1$avg_logFC >= 0],
                     rownames(group2)[group2$avg_logFC >= 0])
down_gene <- intersect(rownames(group1)[group1$avg_logFC < 0],
                     rownames(group2)[group2$avg_logFC < 0])

table_names <- c(paste0("6wk_", colnames(group1)),
                 paste0("13wk_", colnames(group2)))
up_table <- cbind(group1[up_gene, ], group2[up_gene, ])
down_table <- cbind(group1[down_gene, ], group2[down_gene, ])
colnames(up_table) <- table_names
colnames(down_table)  <- table_names
upgenelist[[i]] <- up_table
downgenelist[[i]] <- down_table
}
names(upgenelist) <- sheetname_6wk
names(downgenelist) <- sheetname_6wk
write.xlsx(upgenelist, "female_6wk_up_13wk_up_DEGs_across_cellypes.xlsx", rowNames = T)
write.xlsx(downgenelist, "female_6wk_down_13wk_down_DEGs_across_cellypes.xlsx", rowNames = T)


## 在6wk 13wk中 上下调方向相反的差异基因
upgenelist <- list()
downgenelist <- list()
for (i in 1:length(sheetname_6wk)){
  group1 <- read.xlsx("female_6wk_DEGs_Obesity_VS_Normal.xlsx", 
                      sheet = i, 
                      rowNames = T)
  group2 <- read.xlsx("female_13wk_DEGs_Obesity_VS_Normal.xlsx", 
                      sheet = i, 
                      rowNames = T)
  up_gene <- intersect(rownames(group1)[group1$avg_logFC >= 0],
                       rownames(group2)[group2$avg_logFC < 0])
  down_gene <- intersect(rownames(group1)[group1$avg_logFC < 0],
                         rownames(group2)[group2$avg_logFC >= 0])
  
  table_names <- c(paste0("6wk_", colnames(group1)),
                   paste0("13wk_", colnames(group2)))
  up_table <- cbind(group1[up_gene, ], group2[up_gene, ])
  down_table <- cbind(group1[down_gene, ], group2[down_gene, ])
  colnames(up_table) <- table_names
  colnames(down_table)  <- table_names
  upgenelist[[i]] <- up_table
  downgenelist[[i]] <- down_table
}
names(upgenelist) <- sheetname_6wk
names(downgenelist) <- sheetname_6wk
write.xlsx(upgenelist, "female_6wk_up_13wk_down_DEGs_across_cellypes.xlsx", rowNames = T)
write.xlsx(downgenelist, "female_6wk_down_13wk_up_DEGs_across_cellypes.xlsx", rowNames = T)


########------ 20201127 提取出六种细胞亚型做聚类分析 并按control 和 ob 分开----------------#####
## 画up-down 和down-up的feature plot 图 和 up-up 和 down-down 的feature plot图

# 工作路径
setwd("/Users/zhengyiyi/Desktop/projects/Ob_female_6wk/20201125_Compare_DEGs_with_13wk/")

library("openxlsx")
library("Seurat")
library("dplyr")
subtypes <- read.xlsx("table/20201126_Summary_Compare_DEGs_between_6wk_13wk_marked.xlsx",
                      colNames = T, rowNames = F, sheet = 8)
DEGs <- read.xlsx("table/20201126_Summary_Compare_DEGs_between_6wk_13wk_marked.xlsx",
                  colNames = T, rowNames = F, sheet = 7)

same_direc_DEGs <- DEGs[DEGs$direction=="up-up" | DEGs$direction=="down-down", ]
same_direc_DEGs <- arrange(same_direc_DEGs, desc(direction), Celltypes)

reverse_direc_DEGs <- DEGs[DEGs$direction=="up-down" | DEGs$direction=="down-up", ]
reverse_direc_DEGs <- arrange(reverse_direc_DEGs, desc(direction), Celltypes)

#-----6week 数据
female_6wk <- readRDS("female_6wk_Aggregated_seurat.rds")
metadata_names <- colnames(female_6wk@meta.data)
remove_index <- grep(pattern="prediction", metadata_names)
remove_name <- metadata_names[remove_index]
female_6wk@meta.data <- female_6wk@meta.data[, -which(metadata_names %in% remove_name)]
female_6wk$predicted.id <- NULL
female_6wk$celltype <- female_6wk$celltype_try
female_6wk$celltype_assign <- NULL
saveRDS(female_6wk, file = "female_6wk_Aggregated_seurat.rds")

female_6wk@active.ident <- female_6wk$celltype
subtype_6wk <- subset(x = female_6wk, idents =  subtypes$celltypes)

# clustering 
source("/Users/zhengyiyi/Desktop/code/Seurat_function.R")
source("/Users/zhengyiyi/Desktop/code/diffprop_functions.R")
library(Seurat)
library(ggplot2)

condition <- "6wk_subtype_clustering"
npc_used <- 20
resolution_number <- 0.5
subtype_6wk <- Plot_cluster(subtype_6wk,condition,npc_used,resolution_number)

## cluster 标号重新设置
for (i in 14:1){
  print(levels(subtype_6wk@active.ident)[i])
  levels(subtype_6wk@active.ident)[i] <- as.character(i)
}

library(ggplot2)
p1 <- DimPlot(subtype_6wk, reduction="tsne",
              label = TRUE, pt.size=1.2,
              label.size=4.5, repel=TRUE, group.by = "celltype_try")
p1 <- p1 + labs(title="resolution 0.5") +  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(subtype_6wk, reduction="tsne", label = TRUE,
              pt.size=1.2, label.size=4.5, repel=TRUE, 
              group.by = "celltype_try", split.by="orig.ident")

pdf("female_6wk_subtype_Clustering_tsne_resolution_0.5_plot.pdf",7,7)
print(p1)
dev.off()

pdf("female_6wk_subtype_Clustering_tsne_resolution_0.5_plot_split_by_orig_ident.pdf",12,6)
print(p2)
dev.off()

# up-up, down-down
intertest_DEGs <- same_direc_DEGs  
gene_num <- nrow(same_direc_DEGs)
pic_name <- "female_6wk_only_Subtypes_same_direction_DEGs_Featureplot_plot_marker.pdf"
subtype_6wk <- FeaturePlot_SpecificGene(subtype_6wk, intertest_DEGs, gene_num, pic_name)

# up-down, down-up
intertest_DEGs <- reverse_direc_DEGs  
gene_num <- nrow(reverse_direc_DEGs)
pic_name <- "female_6wk_only_Subtypes_reverse_direction_DEGs_Featureplot_plot_marker.pdf"
subtype_6wk <- FeaturePlot_SpecificGene(subtype_6wk, intertest_DEGs, gene_num, pic_name)
saveRDS(subtype_6wk, file = "female_subtype_6wk.rds")

#----13week 数据
female_13wk <- readRDS("female_13wk_Aggregated_seurat.rds")
metadata_names <- colnames(female_13wk@meta.data)
remove_index <- grep(pattern="prediction", metadata_names)
remove_name <- metadata_names[remove_index]
female_13wk@meta.data <- female_13wk@meta.data[, -which(metadata_names %in% remove_name)]
female_13wk$predicted.id <- NULL
female_13wk$celltype_GSE113576 <- NULL
female_13wk$celltype <- female_13wk$celltype_assign
saveRDS(female_13wk, file = "female_13wk_Aggregated_seurat.rds")
female_13wk@active.ident <- female_13wk$celltype_assign
subtype_13wk <- subset(x = female_13wk, idents =  subtypes$celltypes)

# clustering 
condition <- "13wk_subtype_clustering"
npc_used <- 20
resolution_number <- 0.5
subtype_13wk <- Plot_cluster(subtype_13wk,condition,npc_used,resolution_number)

## cluster 标号重新设置
for (i in 15:1){
  print(levels(subtype_13wk@active.ident)[i])
  levels(subtype_13wk@active.ident)[i] <- as.character(i)
}

library(ggplot2)
p1 <- DimPlot(subtype_13wk, reduction="tsne",
              label = TRUE, pt.size=1.2,
              label.size=4.5, repel=TRUE, group.by = "celltype_assign")
p1 <- p1 + labs(title="resolution 0.5") +  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(subtype_13wk, reduction="tsne", label = TRUE,
              pt.size=1.2, label.size=4.5, repel=TRUE, 
              group.by = "celltype_assign", split.by="orig.ident")
pdf("female_13wk_subtype_Clustering_tsne_resolution_0.5_plot.pdf",7,7)
print(p1)
dev.off()
pdf("female_13wk_subtype_Clustering_tsne_resolution_0.5_plot_split_by_orig_ident.pdf",12,6)
print(p2)
dev.off()

# up-up, down-down
intertest_DEGs <- same_direc_DEGs  
gene_num <- nrow(same_direc_DEGs)
pic_name <- "female_13wk_only_Subtypes_same_direction_DEGs_Featureplot_plot_marker.pdf"
subtype_13wk <- FeaturePlot_SpecificGene(subtype_13wk, intertest_DEGs, gene_num, pic_name)

# up-down, down-up
intertest_DEGs <- reverse_direc_DEGs  
gene_num <- nrow(reverse_direc_DEGs)
pic_name <- "female_13wk_only_Subtypes_reverse_direction_DEGs_Featureplot_plot_marker.pdf"
subtype_13wk <- FeaturePlot_SpecificGene(subtype_13wk, intertest_DEGs, gene_num, pic_name)
saveRDS(subtype_13wk, file = "female_subtype_13wk.rds")

## 合并图片
same_genename <- same_direc_DEGs$GeneName
p_6wk_same <- FeaturePlot(subtype_6wk, features = same_genename, 
                        col = c("lightgrey", "red"),
                        split.by = "orig.ident",
                        reduction = "tsne", min.cutoff = 'q10')
p_13wk_same <- FeaturePlot(subtype_13wk, features = same_genename, 
                          col = c("lightgrey", "red"),
                          split.by = "orig.ident",
                          reduction = "tsne", min.cutoff = 'q10')

reverse_genename <- reverse_direc_DEGs$GeneName
p_6wk_reverse <- FeaturePlot(subtype_6wk, features = reverse_genename, 
                          col = c("lightgrey", "red"),
                          split.by = "orig.ident",
                          reduction = "tsne", min.cutoff = 'q10')
p_13wk_reverse <- FeaturePlot(subtype_13wk, features = reverse_genename, 
                           col = c("lightgrey", "red"),
                           split.by = "orig.ident",
                           reduction = "tsne", min.cutoff = 'q10')

library(cowplot) 
p_same <- cowplot::plot_grid(p_6wk_same, p_13wk_same, ncol = 2)
p_reverse <- cowplot::plot_grid(p_6wk_reverse, p_13wk_reverse, ncol = 2)
pic_name_same <- "female_6wk_13wk_same_direction_DEGs_in_subtypes.pdf"
pic_name_reverse<- "female_6wk_13wk_reverse_direction_DEGs_in_subtypes.pdf"

ggsave(pic_name_same, p_same, width = 20, height = 40)
ggsave(pic_name_reverse, p_reverse, width = 20, height = 40)


#########---------20201129 ob female 6wk 神经元亚型的分析-------------------#########
rm(list=ls())
setwd("/Users/zhengyiyi/Desktop/projects/Ob_female_6wk/")
Aggregated_seurat <- readRDS("female_6wk_neuron.rds")
source("/Users/zhengyiyi/Desktop/code/Seurat_function.R")
source("/Users/zhengyiyi/Desktop/code/diffprop_functions.R")
library(Seurat)
library(ggplot2)
setwd("20201129_Neuron_analysis/")
condition <- "female_6wk_neurons_clustering"
npc_used <- 20
resolution_number <- 0.4
Aggregated_seurat <- Plot_cluster(Aggregated_seurat,condition,npc_used,resolution_number)
## cluster 标号重新设置
for (i in 5:1){
  print(levels(Aggregated_seurat@active.ident)[i])
  levels(Aggregated_seurat@active.ident)[i] <- as.character(i)
}
library(ggplot2)
p1 <- DimPlot(Aggregated_seurat,reduction="tsne",
              label = TRUE,pt.size=1.2,
              label.size=4.5, repel=TRUE)
p1 <- p1 + labs(title="resolution 0.4") +  theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(Aggregated_seurat,reduction="tsne",label = TRUE,
              pt.size=1.2,label.size=4.5, repel=TRUE, split.by="orig.ident")
pdf("female_6wk_neurons_Clustering_tsne_resolution_0.4_plot.pdf",7,7)
print(p1)
dev.off()
pdf("female_6wk_neurons_Clustering_tsne_resolution_0.4_plot_split_by_orig_ident.pdf",12,6)
print(p2)
dev.off()
condition <- "female_6wk_neurons"
Aggregated_seurat_markers <- FindmarkerForCluster(Aggregated_seurat,condition)
save_name <- paste0(condition,"_DEGsbetweenClusters.csv")
write.csv(Aggregated_seurat_markers,save_name)

setwd("/Users/zhengyiyi/Desktop/projects/Ob/Ob_female_6wk/20201129_Neuron_analysis/")
features_1 <- c("Ly6h","Snap25", "Syt1", "Slc17a6", "Slc32a1", "Negr1", 
                "Gad1", "Gad2", "Hcrt", "Pomc", "Bace2", "Apoe", "Slc2a13", "Scg2")
p1 <- VlnPlot(Aggregated_seurat, features_1, ncol = 5, group.by = "test")
ggsave("VlnPlot_of_neuron_makrer_in_neuron_cells.pdf", p1, width = 10, height = 20)
p1 <- FeaturePlot(Aggregated_seurat, features = features_1, ncol = 2, 
                  col = c("lightgrey", "red"), reduction = "tsne")
ggsave("FeaturePlot_of_neuron_makrer_in_neuron_cells.pdf", p1, width = 10, height = 20)
condition <- "female_6wk_neurons"
top_num <- 20
TopMarkersInClusters <- TopMarkersInCluster(Aggregated_seurat_markers,condition,top_num)
p1 <- DoHeatmap(Aggregated_seurat, features = TopMarkersInClusters$gene)
ggsave("Top20Gene_of_eachcluster.pdf", p1, width = 15, height = 15)
a <- read.xlsx("neurons_marker.xlsx", colNames = F, rowNames = F)
p1 <- DoHeatmap(Aggregated_seurat, features = a$X1)

#########---------20201203  & 20201205 神经元亚型细胞类型鉴定---------------------------------######

rm(list=ls())
setwd("/Users/zhengyiyi/Desktop/projects/Ob/Ob_female_6wk/20201129_Neuron_analysis/")
Aggregated_seurat <- readRDS("female_6wk_neuron.rds")
source("/Users/zhengyiyi/Desktop/code/Seurat_function.R")
source("/Users/zhengyiyi/Desktop/code/diffprop_functions.R")
library(Seurat)
library(ggplot2)
## remove cluster 5
selected <- c(1,2,3,4)
subtype_6wk <- subset(x = Aggregated_seurat, idents =  selected)
Aggregated_seurat <- subtype_6wk
rm(subtype_6wk)
## recluster
condition <- "reclustering"
npc_used <- 20
resolution_number <- 0.4
k_param <- 20
Aggregated_seurat <- Plot_cluster(Aggregated_seurat,condition,npc_used,resolution_number,k_para)

genes <- read.csv("HFD_control_DEGs_GABA_VS_Glu.csv", header = T)

# RunTSNE and RunUMAP
Aggregated_seurat <- RunTSNE(Aggregated_seurat,dims = 1:npc_used,verbose = FALSE,
                             check_duplicates = FALSE, seed.use = 3,
                              dim.embed = 3)
Aggregated_seurat <- RunUMAP(Aggregated_seurat, dims = 1:npc_used,verbose = FALSE)

Aggregated_seurat <- FindNeighbors(Aggregated_seurat,
                                  dims = 1:npc_used, 
                                  verbose = FALSE,
                                  k.param = k_param,
                                  nn.method = "annoy",
                                  annoy.metric = "euclidean")
Aggregated_seurat <- FindClusters(Aggregated_seurat,
                                 verbose = FALSE,
                                 resolution = resolution_number)

# BiocManager::install("tidyverse")
library(tidyverse)
tsne = Aggregated_seurat@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>% cbind(tx = Aggregated_seurat@meta.data$orig.ident,
                            clusters = Aggregated_seurat@active.ident)
library(dplyr)
tsne <- arrange(tsne, clusters)
##做成3D图
exp_count=tsne #只放3列数据
library(rgl)
attach(mtcars)
cols=c(rep("#a6cee3",250),rep("#1f78b4",189),rep("#b2df8a",179),
       rep("#33a02c",96))
x=exp_count[,1]
y=exp_count[,2]
z=exp_count[,3]
plot3d(x,y,z,col=cols,size=5)


library(ggplot2)
p1 <- DimPlot(Aggregated_seurat,reduction="tsne",
              label = TRUE,pt.size=1.2,
              label.size=4.5, repel=TRUE)
p1 <- p1 + labs(title="resolution 0.4") +  theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(Aggregated_seurat,reduction="tsne",label = TRUE,
              pt.size=1.2,label.size=4.5, repel=TRUE, split.by="orig.ident")
pdf("female_6wk_neurons_Clustering_tsne_resolution_0.4_plot.pdf",7,7)
print(p1)
dev.off()
pdf("female_6wk_neurons_Clustering_tsne_resolution_0.4_plot_split_by_orig_ident.pdf",12,6)
print(p2)
dev.off()
condition <- "female_6wk_neurons"
Aggregated_seurat_markers <- FindmarkerForCluster(Aggregated_seurat,condition)
save_name <- paste0(condition,"_DEGsbetweenClusters.csv")
write.csv(Aggregated_seurat_markers,save_name)
setwd("/Users/zhengyiyi/Desktop/projects/Ob/Ob_female_6wk/20201129_Neuron_analysis/")
features_1 <- c("Snap25", "Syt1", "Slc17a6", "Slc32a1",  
                "Gad1", "Gad2", "Hcrt", "Pomc",  "Apoe")
features_1 <- c("Snap25", "Syt1", "Slc17a6", "Slc32a1", "Apoe")
p1 <- VlnPlot(Aggregated_seurat, features_1, ncol = 5)
ggsave("VlnPlot_of_neuron_makrer_in_neuron_cells.pdf", p1, width = 10, height = 5)

p1 <- FeaturePlot(Aggregated_seurat, features = features_1, ncol = 2, 
                  col = c("lightgrey", "red"), reduction = "tsne")

ggsave("FeaturePlot_of_neuron_makrer_in_neuron_cells.pdf", p1, width = 10, height = 20)


condition <- "female_6wk_neurons"
top_num <- 20
TopMarkersInClusters <- TopMarkersInCluster(Aggregated_seurat_markers,condition,top_num)

p1 <- DoHeatmap(Aggregated_seurat, features = TopMarkersInClusters$gene)
ggsave("Top20Gene_of_eachcluster.pdf", p1, width = 15, height = 15)

a <- read.xlsx("neurons_marker.xlsx", colNames = F, rowNames = F)
p1 <- DoHeatmap(Aggregated_seurat, features = a$X1)

## assign neurons types
new.cluster.ids <- c("GABA", "GABA", "Glu", "Apoe")
names(new.cluster.ids) <- levels(Aggregated_seurat)
Aggregated_seurat <- RenameIdents(Aggregated_seurat, new.cluster.ids)
DimPlot(Aggregated_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

plot <- DimPlot(Aggregated_seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)
head(select.cells)
Aggregated_seurat <- CellSelector(plot = plot, 
                                  object = Aggregated_seurat,
                                  ident = "Glu")
Aggregated_seurat <- CellSelector(plot = plot, 
                                  object = Aggregated_seurat,
                                  ident = "GABA")
Aggregated_seurat <- CellSelector(plot = plot, 
                                  object = Aggregated_seurat,
                                  ident = "Glu")
p1 <- DimPlot(Aggregated_seurat, reduction = "tsne")
ggsave("b.pdf", p1, width = 6, height = 6)
features_1 <- c("Snap25", "Syt1", "Slc17a6", "Slc32a1",  
                "Gad1", "Gad2", "Apoe", "Hcrt", "Pomc")
p2 <- FeaturePlot(Aggregated_seurat, features_1, reduction = "tsne")
p1 <- VlnPlot(Aggregated_seurat, features_1, ncol = 3, pt.size = 0.5)
ggsave("a.pdf", p2, width = 12, height = 12)

Aggregated_seurat$celltype <- NULL
Aggregated_seurat$celltype_Seurat <- NULL
metadata_names <- colnames(Aggregated_seurat@meta.data)
remove_index <- grep(pattern="expression", metadata_names)
remove_name <- metadata_names[remove_index]
Aggregated_seurat@meta.data <- Aggregated_seurat@meta.data[, -which(metadata_names %in% remove_name)]
Aggregated_seurat$celltype_assign <- Aggregated_seurat@active.ident

saveRDS(Aggregated_seurat, "female_6wk_neuron_remove_cluster5_and_reassign_cell.rds")

## 找神经元亚型之间的ob对obesity之间的差异基因
imprinted_genes <- read.csv("/Users/zhengyiyi/Desktop/ImprintedGene/imprinted_gene_list.csv",
                            header = T)
imprinted_genes$Gene_name <- as.character(imprinted_genes$Gene_name)
cell_type <- sort(unique(Aggregated_seurat$celltype_assign))
treat <- "Obesity"
control <- "Control"
file_prefix <- "female_6wk_neurons"
FindDEGsFromCelltypes(Aggregated_seurat, cell_type, 
                      treat, control, imprinted_genes, file_prefix)

#######-------20201205 ob female 6wk 神经元分析最终的图稿------------###### 
library(pheatmap)
reduction_type <- "umap"
group_by_type <- "celltype_assign"
pdf_prefix <- "Ob_female_6wk"
celltype_markers <- c("Snap25", "Syt1", "Slc17a6", 
                      "Slc32a1", "Gad1", "Gad2", 
                      "Apoe", "Hcrt", "Pomc",
                      "Agrp", "Cartpt", "Npy",
                      "Th", "Slc18a2", "Nr4a2", "Fgf2")
plot_data <- table(Aggregated_seurat$orig.ident, Aggregated_seurat$celltype_assign)
plot_data <- data.frame(Percentage = c(plot_data[1, ]/sum(plot_data[1, ]), 
                                       plot_data[2, ]/sum(plot_data[2, ])),
                        Treat = c(rep("Control", 3), rep("Ob", 3)),
                        Type = rep(colnames(plot_data), 2))
plot_data$Type <- factor(plot_data$Type,
                         levels = c("Apoe", "GABA", "Glu"))
imprinted_genes <- c("Nap1l5", "Ndn", "Peg3", "Snrpn", "Ube3a", "Gnas")

## 画图
NeuronSubtypePlots(Aggregated_seurat, reduction_type, group_by_type, 
                   celltype_markers, plot_data, imprinted_genes)


#######-------20201207 看下多巴胺神经元的表达情况
setwd("/Users/zhengyiyi/Desktop/projects/Ob/Ob_female_6wk/20201129_Neuron_analysis/")
Aggregated_seurat <- readRDS("female_6wk_neuron_remove_cluster5_and_reassign_cell.rds")

