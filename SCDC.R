## 20210113 学习 SCDC R 包
## https://github.com/meichendong/SCDC
## https://meichendong.github.io/SCDC/articles/SCDC.html
## SCDC: Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA Sequencing References

## STEP 1: install 
## github 上这样安不上 后来我选择自己下载源代码自己安装
## 再安装几个依赖包就可以了

if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("meichendong/SCDC")


install.packages("/Users/zhengyiyi/Desktop/SCDC-master/",
                 repos=NULL,type="source")

library(SCDC)
library(dplyr)

## STEP 2:SCDC Pre-process of scRNA-seq Data
## Run example data fot test


ct1 <- c("mediumorchid1","mediumpurple1","lightskyblue","seagreen1","yellow","tan1","azure3")
seger <- readRDS("/Users/zhengyiyi/Desktop/code/R/SCDC/SCDC-master/vignettes/data/segerstolpe.rds")
DemoPlot(seger, cluster = "cluster", sample = "sample", select.ct = c("alpha","beta","delta","gamma","ductal","acinar"), Palette = ct1)

seger.qc <- SCDC_qc(seger, ct.varname = "cluster", sample = "sample", scsetname = "Segerstolpe",
                    ct.sub = c("alpha","beta","delta","gamma","ductal","acinar"), qcthreshold = 0.7) 
seger.qc$heatfig


qc.seger <- readRDS("/Users/zhengyiyi/Desktop/code/R/SCDC/SCDC-master/vignettes/data/qc_segerstolpe.rds")
qc.baron <- readRDS("/Users/zhengyiyi/Desktop/code/R/SCDC/SCDC-master/vignettes/data/qc_baron.rds")
qc.xin <- readRDS("/Users/zhengyiyi/Desktop/code/R/SCDC/SCDC-master/vignettes/data/qc_xin.rds")

fadista_77 <- readRDS("/Users/zhengyiyi/Desktop/code/R/SCDC/SCDC-master/vignettes/data/fadista_77.rds")
# setting the search.length as 0.01 
# might take several minutes to finish the ENSEMBLE procedure.
fadista.healthy.ens <- SCDC_ENSEMBLE(bulk.eset = fadista_77[,fadista_77$hba1c_class2 == "Normal"], 
                                     sc.eset.list = list(baronh = qc.baron$sc.eset.qc, 
                                                         segerh = qc.seger$sc.eset.qc), ct.varname = "cluster",
                                     sample = "sample", truep = NULL, ct.sub =  c("alpha","beta","delta","gamma","acinar","ductal"), search.length = 0.01, grid.search = T)  

fadista.t2d.ens <- SCDC_ENSEMBLE(bulk.eset = fadista_77[,fadista_77$hba1c_class2 == "T2D"], sc.eset.list = list(baronh = qc.baron$sc.eset.qc, segerh = qc.seger$sc.eset.qc), ct.varname = "cluster",
                                 sample = "sample", truep = NULL, ct.sub =  c("alpha","beta","delta","gamma","acinar","ductal"), search.length = 0.01, grid.search = T)  




library(reshape2)
library(ggplot2)
getPropBox <- function(ens_h, ens_d, metric = NULL, ref, input_wt = NULL){
  if (!is.null(input_wt)){
    prop.h = as.data.frame(wt_prop(input_wt, ens_h$prop.only))
    prop.d = as.data.frame(wt_prop(input_wt, ens_d$prop.only))
  } else {
    prop.h = as.data.frame(wt_prop(ens_h$w_table[metric,1:2], ens_h$prop.only))
    prop.d = as.data.frame(wt_prop(ens_d$w_table[metric,1:2], ens_d$prop.only))
  }
  prop2 <- rbind(prop.h, prop.d)
  prop2$condition <- c(rep("Normal", nrow(prop.h)), rep("T2D", nrow(prop.d)))
  
  dtmelt <- melt(prop2, id.vars = "condition")
  dtmelt$ref <- as.factor(ref)
  return(dtmelt)
}

fdt.ens.spearman <- getPropBox(ens_h = fadista.healthy.ens,
                               ens_d = fadista.t2d.ens, metric = 1,
                               ref = "ENSEMBLE+SpearmanR")
fdt.seger <- getPropBox(ens_h = fadista.healthy.ens, 
                        ens_d = fadista.t2d.ens,  ref = "Segerstolpe", 
                        input_wt = c(0,1))
fdt.baron <- getPropBox(ens_h = fadista.healthy.ens, ens_d = fadista.t2d.ens,  
                        ref = "Baron", input_wt = c(1,0))

dtall_fadista <- rbind(fdt.ens.spearman, fdt.seger, fdt.baron)
dtall_fadista$refcond <- paste(dtall_fadista$ref, dtall_fadista$condition)
colfunc <- colorRampPalette(c("red", "white"))
colfunc2 <- colorRampPalette(c("blue", "white"))
pfa2 <- ggplot(dtall_fadista[dtall_fadista$refcond %in% c("ENSEMBLE+SpearmanR Normal", "Segerstolpe Normal","Baron Normal",
                                                          "ENSEMBLE+SpearmanR T2D","Segerstolpe T2D","Baron T2D"),], 
               aes(x=variable, y=value, color = factor(refcond, levels = c("ENSEMBLE+SpearmanR Normal", "Segerstolpe Normal","Baron Normal",
                                                                           "ENSEMBLE+SpearmanR T2D","Segerstolpe T2D","Baron T2D")))) +
  geom_boxplot(outlier.size=-1)+
  geom_jitter(aes(x=variable,
                  color = factor(refcond,
                                 levels = c("ENSEMBLE+SpearmanR Normal", "Segerstolpe Normal","Baron Normal",
                                            "ENSEMBLE+SpearmanR T2D","Segerstolpe T2D","Baron T2D"))),
              position = position_dodge(0.75), alpha = 0.5,cex = 0.25)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size=10),
        axis.text.y = element_text(size = 10),
        text = element_text(size = 10),
        plot.title = element_text(size=10, face = "bold"),
        plot.margin=unit(c(1,1,-5,0), "mm"),
        legend.position="top",legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.box.spacing = unit(0, "mm"))+
  scale_color_manual(values=c(colfunc(10)[seq(1,9,3)],colfunc2(10)[seq(1,9,3)])) +
  ylab("") + xlab("")
pfa2


qc.3cl <- readRDS("/Users/zhengyiyi/Desktop/code/R/SCDC/SCDC-master/vignettes/data/MIX3cl_scESET.rds")

sc3cl.basis <- SCDC_basis_ONE(qc.3cl$sc.eset.qc, ct.varname = "md_cluster", sample ="orig.ident")
df.d <- dist(t(log(sc3cl.basis$basis.mvw + 1e-10)), method = "euclidean")
hc1 <- hclust(df.d, method = "complete")
plot(hc1, cex = 0.8, hang = -1, main = "log(basis matrix)")