## 20201201  scLearn
## 两个依赖的R包
BiocManager::install("M3Drop")
BiocManager::install("SingleCellExperiment")

library(devtools)
install_github("bm2-lab/scLearn")
library("scLearn")

##-----1.Data preprocessing

# loading the reference dataset
data<-readRDS('/Users/zhengyiyi/Desktop/projects/GSE113576_ref_seurat.rds')

data<-readRDS("projects/Ob/Data/Cell_type_identification_data/ExoExp_RefSeuratObject/GSE113576_ref_seurat.rds")
rawcounts <- data@assays$RNA@counts
refe_ann<-as.character(data$celltype)
names(refe_ann)<-colnames(data)
# cell quality control and rare cell type filtered and feature selection
data_qc<-Cell_qc(rawcounts,refe_ann,species="Mm")
saveRDS(data_qc, "data_qc.rds")
data_type_filtered<-Cell_type_filter(data_qc$expression_profile,
                                     data_qc$sample_information_cellType,
                                     min_cell_number = 10)

high_varGene_names <- Feature_selection_M3Drop(data_type_filtered$expression_profile)

# training the model
scLearn_model_learning_result<-scLearn_model_learning(high_varGene_names,data_type_filtered$expression_profile,data_type_filtered$sample_information_cellType,bootstrap_times=10)
saveRDS(scLearn_model_learning_result, "scLearn_model_learning_result.rds")
# Cell assignment
# loading the quary cell and performing cell quality control.
data2<-readRDS('female_Aggregated_seurat.rds')
rawcounts2 <- data2@assays$RNA@counts
#query_ann<-as.character(data2$cell_type1)
#names(query_ann)<-colnames(data2)
#query_ann<-query_ann[query_ann %in% c("alpha","beta","delta","gamma")]
#rawcounts2<-rawcounts2[,names(query_ann)]
#data_qc_query<-Cell_qc(rawcounts2,query_ann,species="Hs")
### 
data_qc_query<-Cell_qc(rawcounts2,species="Mm",gene_low=0,umi_low=0,mito_high = 40)
# Assignment with trained model above
scLearn_predict_result<-scLearn_cell_assignment(scLearn_model_learning_result,data_qc_query$expression_profile)
saveRDS(scLearn_predict_result, "scLearn_predict_result.rds")

sum(names(data2$celltype) == scLearn_predict_result$Query_cell_id)
data2$celltype_scLearn <- scLearn_predict_result$Predict_cell_type

sum(data2$celltype_scLearn == data2$celltype)
saveRDS(data2, "female_Aggregated_seurat.rds")
write.csv(data2@meta.data, "female_6wk_meatadata.csv")

