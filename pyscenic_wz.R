
cd /home/wangzhe/miniconda3/bin
source activate
conda activate R

#

library(Seurat)
setwd("/home/zhengjh/other/pan")
list.files()
mouse <- readRDS("mouse3pan_int_all_new.rds")
counts <- mouse@assays$RNA@counts
rowsums <- rowSums(counts)
counts <- counts[]
write.csv(t(as.matrix(data@assays$RNA@counts)),file = "data.csv")

############# python  ###########
# pyscenic 我安装在 env bioinfo下面了
conda activate bioinfo

python
#
import os, sys
os.getcwd()
os.listdir(os.getcwd()) 
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs);


##### 三步走 ####
###############
nohup pyscenic grn --num_workers 10 --output adj.sample.tsv --method grnboost2 sample.loom hs_hgnc_tfs.txt &   #转录因子文件

################
nohup pyscenic ctx adj.sample.tsv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname sample.loom --mode "dask_multiprocessing" --output reg.csv --num_workers 3 --mask_dropouts &

################
pyscenic aucell sample.loom reg.csv --output sample_SCENIC.loom --num_workers 3


###### 可视化  ######
rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")  
library(SCopeLoomR)
scenicLoomPath='sample_SCENIC.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
sce=readRDS("./fibo_1000.rds")
sce
library(pheatmap) 
n=t(scale(t( getAUC(regulonAUC[,] )))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
dim(n) 
ac=data.frame(group= as.character( Idents(sce)))
n[1:4,1:4]
n=n[,colnames(n) %in% colnames(sce)]
rownames(ac)=colnames(n) 
cg=read.table('choose_tf.txt')[,1]
cg
cg_n=n[rownames(n) %in% cg,]
pheatmap(cg_n,show_colnames =F,show_rownames = T,
         annotation_col=ac)
table(ac$group)
# 尊重作者，进行二分类！
ac$group=ifelse(ac$group %in% c(2:5,7,9),'mCAF','iCAF')
pheatmap(cg_n,show_colnames =F,show_rownames = T,
         annotation_col=ac)
pheatmap(cg_n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         filename = 'heatmap_choose_regulon.png')
dev.off()

