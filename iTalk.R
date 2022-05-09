
###################################  iTALK  #####################

# BiocManager::install("iTALK")
install.packages("/Users/zhengyiyi/Desktop/iTALK-master/",repos=NULL,type="source")

install.packages("scde")
BiocManager::install("scde")
BiocManager::install("DEsingle")
BiocManager::install("MAST")
BiocManager::install("scater")
# github上下载
install.packages("/Users/zhengyiyi/Desktop/iTALK-master/",repos=NULL,type="source")

# if(!require(devtools)) install.packages("devtools");
# devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)

library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
library(network)
library(scde)
library(DEsingle)
library(MAST)
library(scater)


setwd("~/Desktop/Ob")
load("0_ALL_TSNE.Rdata")
OB_free=subset(readsCountSM_TSNE,cells=as.vector(read.table("OB_sample.txt")[,1]))

iTalk_data <- as.data.frame(t(OB_free@assays$RNA@counts))

# iTALK Ҫ??????cell_type?У??ҵ?ϸ????Ⱥ?洢??seurat_cluster
iTalk_data$cell_type <- OB_free@meta.data$celltype

# iTALK 要对比gender之间
iTalk_data$compare_group <- OB_free@meta.data$age_class

unique(iTalk_data$cell_type)

unique(iTalk_data$compare_group)

my10colors <- my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282')
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")

#comm_list<-c('growth factor','other','cytokine','checkpoint')
comm_list<-c('growth factor','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='DEG', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}

iTalk_res1 <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:100,]        

iTalk_res1=read.csv("iTalk_res_rank_nocancerfree.csv")
pdf("iTalk_res_NetView_20_100_nocancerfree.pdf",6,6)
print(NetView(iTalk_res1[1:100,],col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5))
dev.off()
######## Groups#####

deg_nk<-DEG(iTalk_data %>% filter(cell_type=='AP'),method='Wilcox',contrast=c('young', 'old'))
deg_AP=deg_nk
deg_B<-DEG(iTalk_data %>% filter(cell_type=='B'),method='Wilcox',contrast=c('young', 'old'))
deg_T<-DEG(iTalk_data %>% filter(cell_type=='T'),method='Wilcox',contrast=c('young', 'old'))
deg_NK<-DEG(iTalk_data %>% filter(cell_type=='NK'),method='Wilcox',contrast=c('young', 'old'))
deg_Granulocytes<-DEG(iTalk_data %>% filter(cell_type=='Granulocytes'),method='Wilcox',contrast=c('young', 'old'))
deg_Macrophage<-DEG(iTalk_data %>% filter(cell_type=='Macrophage'),method='Wilcox',contrast=c('young', 'old'))
deg_DC<-DEG(iTalk_data %>% filter(cell_type=='DC'),method='Wilcox',contrast=c('young', 'old'))
deg_Endothelial<-DEG(iTalk_data %>% filter(cell_type=='Endothelial'),method='Wilcox',contrast=c('young', 'old'))
deg_Pericytes<-DEG(iTalk_data %>% filter(cell_type=='Pericytes/SMC'),method='Wilcox',contrast=c('young', 'old'))
deg_HSC<-DEG(iTalk_data %>% filter(cell_type=='HSC'),method='Wilcox',contrast=c('young', 'old'))
#
all_list=rbind(deg_AP,deg_B,deg_T, deg_NK,deg_Granulocytes,deg_Macrophage,deg_DC,
               deg_Endothelial,deg_Pericytes,deg_HSC)
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(all_list,datatype='DEG',comm_type=comm_type)
  
  res<-rbind(res,res_cat)
}
write.csv(res,"AP_allcelltype_ageclass_DEG_italk.csv")
######做一个图
# FindLR DEG类型的数据，可以输入一个基因集合，结果为相应基因内的配体-受体列表
# 如果有超过20组配体-受体结果，取前20进行展示
pdf("AP_allcelltype_ageclass_DEG_italk.pdf")
res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]

print(LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC))
dev.off()

#########

##各个配体-受体分类分别作图

res<-NULL
for(comm_type in comm_list){
  #comm_type='other' 
  #多个细胞类型之间显著表达的配体-受体，结果会过滤同一细胞类型内的配体-受体关系
  res_cat <- FindLR(deg_3, datatype='DEG', comm_type=comm_type)
  #单个细胞类型显著表达的配体-受体
  #res_cat <- FindLR(data_1=deg_cd4T, datatype='DEG', comm_type=comm_type)
  res_cat <- res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC, decreasing=T),]
  write.csv(res_cat, paste0('LRpairs_DEG_',comm_type,'.xls'))
  #plot by ligand category
  png(paste0('LRpairs_DEG_',comm_type,'_circ.png'), width=600, height=650)
  if(nrow(res_cat)==0){
    next
  }else if(nrow(res_cat>=PN)){
    LRPlot(res_cat[1:PN,], datatype='DEG', link.arr.lwd=res_cat$cell_from_logFC[1:PN],
           cell_col=cell_col, link.arr.width=res_cat$cell_to_logFC[1:PN])
  }else{
    LRPlot(res_cat, datatype='DEG', link.arr.lwd=res_cat$cell_from_logFC,
           cell_col=cell_col, link.arr.width=res_cat$cell_to_logFC)
  }
  dev.off()
  png(paste0('LRpairs_DEG_',comm_type,'_net.png'), width=600, height=650)
  NetView(res_cat, col=cell_col, vertex.label.cex=1.2, edge.label.cex=0.9, 
          vertex.size=30, arrow.width=3, edge.max.width=10, margin = 0.2)
  dev.off()
  res<-rbind(res,res_cat)
}
