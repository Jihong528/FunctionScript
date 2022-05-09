BiocManager::install("sva")
BiocManager::install("bladderbatch")

library(bladderbatch)
library(FactoMineR)
library(factoextra)


library("sva")
### combat
head(int@assays$RNA@data)
nordata <- as.data.frame(int@assays$RNA@data)
intdata <- as.data.frame(int@assays$integrated@data)
head(nordata[1:10])
dim(nordata)
head(nordata[1:10])
saminfo <- read.csv("sample_info.csv", row.names = 1)
head(saminfo)

pca.plot = function(dat,col){
  
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}
p1=pca.plot(nordata,factor(saminfo$geo))
pdf("pca_batchremove_bf.pdf", 18, 12) 
p1
dev.off()

p2=pca.plot(intdata,factor(saminfo$geo))
pdf("pca_batchremove_af.pdf", 18, 12) 
p2
dev.off()

#设置model（可选）
mod = model.matrix(~as.factor(species), data=saminfo)
#校正
combat_tpm <- ComBat(dat = as.matrix(log2(nordata+1)), batch = saminfo$batch_GSE, par.prior = F)

p2=pca.plot(combat_tpm,factor(sample_infor$GSE_name))

