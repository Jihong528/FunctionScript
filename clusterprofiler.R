library(clusterProfiler)
BiocManager::install("org.Mm.eg.db") # 安装人类的注释数据
library(org.Mm.eg.db) # 测试，载入人类的注释数据

data <- read.table("DEGs.txt",header=T,row.names = 1) # 从文件中载入数据，并赋值给data
# 我们需要的是第二列的基因列表，因此，只筛选第二列即可
gene <- rownames(data) # 抽取第二列，赋值给gene
gene <- bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db) # groupGO默认使用"ENTREZID"作为输入数据，因此首先将“SYMBOL”转化为“ENTREZID”，以防出错
gene <- gene[,2] # 抽取第二列作为input

ego=enrichGO(OrgDb="org.Hs.eg.db", gene = gene,pvalueCutoff = 0.01,readable=TRUE)#一行代码完成GO富集，gene1就是上图的entrezID的list
write.csv(as.data.frame(ego),"GO-enrich.csv",row.names =F) #写入文件

ekk <- enrichKEGG(gene = gene,organism = 'hsa',pvalueCutoff = 0.05)
write.csv(as.data.frame(ekk),"KEGG-enrich.csv",row.names =F)
dotplot(ego,showCategory=10,title="Enrichment GO Top10") #泡泡图
barplot(ego, showCategory=20,title="EnrichmentGO")  #柱状图




#GSEA
library(dplyr)
#对数据进行预处理-排序
data<-markers_gsea[markers_gsea$cluster=='subtype2',]
gene <- data$gene
Fold.Change<-as.data.frame(data$avg_logFC)
gene<-as.data.frame(gene)
data_list<-cbind(gene,Fold.Change)
data.sort <- arrange(data_list, desc(data$avg_logFC)); head(data.sort)
gene <- data_list$gene
gene <- bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db) # groupGO默认使用"ENTREZID"作为输入数据，因此首先将“SYMBOL”转化为“ENTREZID”，以防出错
genelist <- gene[,2]# 抽取第二列作为input

#使用broad GSEA提供的gene sets 来提供TERM2GENE
gmtfile <- system.file("extdata", "c5.all.v7.2.entrez.gmt", package="clusterProfiler")
c2 <- read.gmt(gmtfile)
head(c2)
### 先使用基于超几何分布的富集分析
enrich <- enricher(genelist, TERM2GENE=c2); head(enrich)
#再做富集分析
## assume 1st column is ID
## 2nd column is FC
head(data_list)
## feature 1: numeric vector
glist <- data_list[,2];head(glist)
## feature 2: named vector
names(glist) <- as.character(genelist);head(glist)
## feature 3: decreasing order
glist <- sort(glist,decreasing = T); head(glist)
gsea <- GSEA(glist, TERM2GENE=c2, verbose=FALSE, pvalueCutoff = 1); head(gsea)
write.csv(as.data.frame(gsea),"GSEA_c5_subtype2.csv",row.names =F)
gseaplot2(gsea,geneSetID=2)
#
library(ReactomePA)
Go_gseresult <- gseGO(glist, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
KEGG_gseresult <- gseKEGG(glist, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
Go_Reactomeresult <- gsePathway(glist, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gseaplot2(Go_gseresult, 1)



#################
library(Pi)
data<-data_list
rownames(data)<- data_list$gene
eGSEA <-
  Pi::xPierGSEA(
    data,
    ontology = "MsigdbH",
    size.range = c(20, 5000),
    nperm = 20000,
    fast = F,
    RData.location = "data/ontology_Rdata")
