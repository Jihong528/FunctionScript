library(Seurat)
library(celltalker)
library(stringr)
##########################1. 读入seurat对象##########################
seurat_obj <- readRDS("D:/scRNAseq/hy/GPCR/result2/OB_seurat/seurat_TSNE.rds")

#调整meta.data，方便满足celltalker的数据要求
#对照和处理需要多个组
cell_meta = seurat_obj@meta.data
names(cell_meta)[which(names(cell_meta)=='orig.ident')] <- "sample.type"
names(cell_meta)[which(names(cell_meta)=='celltype')] <- "cell.type"
cell_meta$sample.id <- c(paste0(seurat_obj$orig.ident[1:1000],"_1"),paste0(seurat_obj$orig.ident[1001:2160],"_2"),paste0(seurat_obj$orig.ident[2161:3660],"_1"),paste0(seurat_obj$orig.ident[3661:5435],"_2"))
seurat_obj@meta.data <- cell_meta

######################2. 鉴定差异表达配体和受体######################
##ramilowski_pairs是配体受体数据集
#如果自己有受体配体对，也可以整理成这样一个三列的格式，导入进来
a <-ramilowski_pairs
a$ligand <- str_to_title(a$ligand)
a$receptor <- str_to_sentence(a$receptor)
a$ligand <- paste0(a$ligand,"_")
write.csv(a,"lr.csv",quote = F)

ramilowski_pairs <- read.csv("lr.csv",row.names = 1)
ligs <- as.character(unique(ramilowski_pairs$ligand))
recs <- as.character(unique(ramilowski_pairs$receptor))

ligs.present <- rownames(seurat_obj)[rownames(seurat_obj) %in% ligs]
recs.present <- rownames(seurat_obj)[rownames(seurat_obj) %in% recs]
genes.to.use <- union(ligs.present,recs.present) 

##差异表达配体受体
Idents(seurat_obj) <- "sample.type"
markers <- FindAllMarkers(seurat_obj, assay="RNA",features=genes.to.use,only.pos=TRUE)
ligs.recs.use <- unique(markers$gene)
length(ligs.recs.use)

interactions.forward1 <- ramilowski_pairs[as.character(ramilowski_pairs$ligand) %in% ligs.recs.use,]
interactions.forward2 <- ramilowski_pairs[as.character(ramilowski_pairs$receptor) %in% ligs.recs.use,]
interact.for <- rbind(interactions.forward1, interactions.forward2)
dim(interact.for)

######################3. 创建一致表达的配体受体tibble格式数据######################
#Create data for celltalker
expr.mat <- GetAssayData(seurat_obj,slot="counts")
defined.clusters <- seurat_obj@meta.data$cell.type
defined.groups <- seurat_obj@meta.data$sample.type
defined.replicates <- seurat_obj@meta.data$sample.id

reshaped.matrices <- reshape_matrices(count.matrix=expr.mat,clusters=defined.clusters,groups=defined.groups,replicates=defined.replicates,ligands.and.receptors=interact.for)
unnest(reshaped.matrices,cols="samples")
names(pull(unnest(reshaped.matrices,cols="samples"))[[1]])
#为每个组创建一致表达的配体和受体 
# cells.reqd=10：每个cluster中至少有10个细胞表达了配体/受体
# freq.pos.reqd=0.5：至少有50％重复个体中表达的配体/受体
consistent.lig.recs <- create_lig_rec_tib(exp.tib=reshaped.matrices,
                                          clusters=defined.clusters,groups=defined.groups,
                                          replicates=defined.replicates,
                                          cells.reqd=0,freq.pos.reqd=0,
                                          ligands.and.receptors=interact.for)
unnest(consistent.lig.recs[1,2],cols="lig.rec.exp")

######################4. 生成celltalk对象并可视化######################
# freq.group.in.cluster: 只对包含细胞数大于总细胞数5%的簇进行互作分析
put.int <- putative_interactions(ligand.receptor.tibble=consistent.lig.recs,
                                 clusters=defined.clusters,groups=defined.groups,
                                 freq.group.in.cluster=0.2,
                                 ligands.and.receptors=interact.for)

##鉴定分组间特异出现的配体/受体并作图
unique.ints <- unique_interactions(put.int,group1="Control",group2="Obesity",interact.for)

#plot
control.to.plot <- pull(unique.ints[1,2])[[1]]
for.circos.control <- pull(put.int[1,2])[[1]][control.to.plot]
circos_plot(interactions=for.circos.control,clusters=defined.clusters)

obesity.to.plot <- pull(unique.ints[1,2])[[1]]
for.circos.obesity <- pull(put.int[1,2])[[1]][obesity.to.plot]
circos_plot(interactions=for.circos.obesity,clusters=defined.clusters)


pdf("ligand_receptor_interaction.pdf",height = 6, width = 6)
print(circos_plot(interactions=for.circos.control,clusters=defined.clusters))
print(circos_plot(interactions=for.circos.obesity,clusters=defined.clusters))
dev.off()

