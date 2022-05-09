# This is a pipeline to convert the ensemble ID to genename
# Version: 2020-12-03; Jihong Zheng
rm(list=ls())
if (!require('biomaRt')  ) 
{
  print("Install the package: stringi")
  BiocManager::install('biomaRt')
}
library(biomaRt)
listMarts()
CovertEnsembl2Genename <- function(filename,EnsemblDataBase){
ensembl_gene_id <- read.csv(filename, header = T)
ensembl_gene_id <- ensembl_gene_id[,1]
print("Number of ensemble gene id")
print(length(ensembl_gene_id))
# 1. use useMart to select used database
ensembl<-useMart("ensembl")
# 2. 用listDatasets()函数显示当前数据库所含的基因组注释
Datasets <- listDatasets(ensembl)
# 3. 用useDataset()函数选定数据库中的基因组
# 选定ensembl数据库中的Anole lizard genes (AnoCar2.0)基因组
mart <- useDataset(EnsemblDataBase, useMart("ensembl"))
# 4. 选定我们需要获得的注释类型
# 用lsitFilters()函数查看可选择的类型，选定要获取的注释类型，以及已知注释的类型
Annotate_type <- listFilters(mart)
# 5.用getBM()函数获取注释
gene_symbols<- getBM(attributes=c('ensembl_gene_id',
                                    'external_gene_name'),
                       filters= 'ensembl_gene_id',
                       values = ensembl_gene_id, mart = mart)
write.csv(gene_symbols, "EnsemblID2Genename.csv")
}

Args <- commandArgs(TRUE)
filename <- Args[1]
EnsemblDataBase <- Args[2]
CovertEnsembl2Genename(filename, EnsemblDataBase)

## test
## setwd("/Users/zhengyiyi/Desktop/projects/Ob/SmartSeq2_Ob_P7/")
## ensembl_gene_id <- read.csv("P7_ob_ensemble_id.csv",header = T)$ensemble_id
## EnsemblDataBase <- "mmusculus_gene_ensembl"
## filename='P7_ob_ensemble_id.csv'
## EnsemblDataBase <- "mmusculus_gene_ensembl"