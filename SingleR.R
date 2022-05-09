###singleR ANNO
ref <- readRDS("GSE113576_ref_seurat.rds")

refdata <- GetAssayData(ref, slot="data")
testdata <- GetAssayData(mob, slot="data")
clusters <- mob@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = ref@meta.data$celltype, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

for(i in 1:nrow(celltype)){
  mob@meta.data[which(mob@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(mob@meta.data$celltype)