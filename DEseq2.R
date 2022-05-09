library(Seurat)
library(openxlsx)
setwd("/home/zhengjh/projects/Ob/Results/20210110_ob_female_13wk_2")
female <- readRDS("20210110_Aggregated_seurat.rds")
DefaultAssay(female) <- "RNA"
female[["RNA"]]@counts<-as.matrix(female[["RNA"]]@counts)+1
female_rep2_DEGS <- FindMarkers(female, ident.1 = "ob/ob", ident.2 = "WT",
                                test.use = "DESeq2", slot = "counts", group.by = "orig.ident")
log2FC <- log2(exp(female_rep2_DEGS$avg_logFC))
female_rep2_DEGS$avg_log2FC <- log2FC
write.xlsx(female_rep2_DEGS, "try2_singlecell_female_rep2_DEGS_ob2ctrl_DEseq2.xlsx", 
           rowNames = T)

