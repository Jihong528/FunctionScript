#if (!requireNamespace("BiocManager", quietly = TRUE)) 
#install.packages("BiocManager")

#BiocManager::install(c("GENIE3", "AUCell", "RcisTarget","NMF","rbokeh","mixtools",
#"pheatmap", "Rtsne", "R2HTML","SingleCellExperiment"),force = TRUE)

#install.packages(c("zoo","R2HTML","doMC","doRNG")
#install.packages("devtools")
#devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#devtools::install_github("aertslab/SCENIC", ref="v1.1.0",force = TRUE)

#
library(leidenbase)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(zoo)
library(mixtools)
library(rbokeh)
library(NMF)
library(pheatmap)
library(Rtsne)
library(SCopeLoomR)      
library(SingleCellExperiment)
library(R2HTML)
library(doMC)  
library(doRNG)
library(SCENIC) 

setwd("/home/wangzhe/Dengxuan")
#a=load("0_COAD_all.Rdata")
#Epi=subset(scRNA_TSNE,celltype11=="Epi")
load("Epi.Rdata")
Epi=subset(Epi, nFeature_RNA>1000 & percent.mt < 15)
singleCellMatrix <- as.matrix(Epi@assays$RNA@counts)
singleCellMatrix= singleCellMatrix[rowSums(singleCellMatrix)>100,]
cellInfo <- data.frame(seuratCluster=Epi@meta.data$location) #CellType
colnames(cellInfo) <- "CellType"
rownames(cellInfo) <- rownames(Epi@meta.data)
cbind(table(cellInfo$CellType))
########


colVars <- list(CellType=c("L"="forestgreen",
                           "N"="darkorange","D45493C"="magenta4",

                   "C"="magenta4"))
                  #  "D10375N"="hotpink",
                #  "D16966N"="red3",
                #  "D16966L"="yellow","D10375L"="black", "D4549L" ="black"
                #  ))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")

#
library(SCENIC)
org="hgnc" # hgnc, mgi, or dmel
#dbDir="C:/Users/WZ/Desktop/CTC/SCIENIC" # RcisTarget databases location
#dbDir="/share/home/wz20/Project/CTC/SCIENIC"
#dbDir="/Users/wangzhe/Desktop/CTC1/SCIENIC"
dbDir="/home/wangzhe/SCENIC"
myDatasetTitle="SCENIC on BCE" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,nCores=4)
scenicOptions@settings[["dbs"]][["500bp"]]="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
scenicOptions@settings[["dbs"]][["10kb"]]="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"


#
# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
####
# (Adjust minimum values according to your dataset)
exprMat=singleCellMatrix
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))
#??ע?????Ƿ???ȥ??
#interestingGenes <- c("SOX9", "SOX10", "DLX5")
#interestingGenes[which(!interestingGenes %in% genesKept)]
#
# Run GENIE3
exprMat_filtered <- exprMat[genesKept, ]
#memory.limit(15000) #��??һ???ڴ?
corrMat <- cor(t(exprMat_filtered), method="spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))

exprMat_filtered_log <- log2(exprMat_filtered+1) 
data(motifAnnotations_hgnc) 
motifAnnotations_hgnc_v8=motifAnnotations_hgnc
runGenie3(exprMat_filtered_log, scenicOptions) #????GENIE3?õ?Ǳ??ת¼????TF


library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 4
scenicOptions@settings$seed <- 123
# For a very quick run:
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...
runSCENIC_1_coexNetwork2modules(scenicOptions)
#???￪ʼ
runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- initializeScenic(org="hgnc", dbDir=dbDir,nCores=1)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log) 

################

.openDev <- function(fileName, devType, ...)
{
  if(devType=="pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType=="png")
    png(paste0(fileName, ".png", type="cairo"), ...)
  
  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
}

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}

runSCENIC_4_aucell_binarize <- function(scenicOptions, 
                                        skipBoxplot=FALSE, skipHeatmaps=FALSE, skipTsne=FALSE, exprMat=NULL)
{
  nCores <- getSettings(scenicOptions, "nCores")
  regulonAUC <- tryCatch(loadInt(scenicOptions, "aucell_regulonAUC"),
                         error = function(e) {
                           if(getStatus(scenicOptions, asID=TRUE) < 3) 
                             e$message <- paste0("It seems the regulons have not been scored on the cells yet. Please, run runSCENIC_3_scoreCells() first.\n", 
                                                 e$message)
                           stop(e)
                         })
  thresholds <- loadInt(scenicOptions, "aucell_thresholds")
  thresholds <- getThresholdSelected(thresholds)
  
  # Assign cells
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(regulonAUC)[x,]>trh))
                                   }),names(thresholds))
  ### Convert to matrix (regulons with zero assigned cells are lost)
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"
  saveRDS(binaryRegulonActivity, file=getIntName(scenicOptions, "aucell_binary_full"))
  
  # Keep only non-duplicated thresholds
  # (e.g. only "extended" regulons if there is not a regulon based on direct annotation)
  binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDuplicatedExtended(rownames(binaryRegulonActivity))),]
  saveRDS(binaryRegulonActivity_nonDupl, file=getIntName(scenicOptions, "aucell_binary_nonDupl"))
  
  minCells <- ncol(binaryRegulonActivity) * .01
  msg <- paste0("Binary regulon activity: ",
                nrow(binaryRegulonActivity_nonDupl), " TF regulons x ",
                ncol(binaryRegulonActivity), " cells.\n(",
                nrow(binaryRegulonActivity), " regulons including 'extended' versions)\n",
                sum(rowSums(binaryRegulonActivity_nonDupl)>minCells),
                " regulons are active in more than 1% (", minCells, ") cells.")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  
  if(!skipBoxplot)
  {
    .openDev(fileName=getOutName(scenicOptions, "s4_boxplotBinaryActivity"),
             devType=getSettings(scenicOptions, "devType"))
    par(mfrow=c(1,2))
    boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon",
            sub='number of cells \nthat have the regulon active',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell",
            sub='number of regulons \nactive per cell',
            col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
    dev.off()
  }
  
  ################################################################################
  # Binary activity heatmap
  if(!skipHeatmaps)
  {
    regulonSelection <- loadInt(scenicOptions, "aucell_regulonSelection", ifNotExists="null", verbose=FALSE)
    if(is.null(regulonSelection)) 
      regulonSelection <- regulonSelections(binaryRegulonActivity, binaryRegulonActivity_nonDupl, minCells)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
    cellInfo <- data.frame(cellInfo)
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    
    ### Plot heatmap:
    for(selRegs in names(regulonSelection$labels))
    {
      if(length(regulonSelection[[selRegs]])>0)
      {
        regulonSelection[[selRegs]] <- regulonSelection[[selRegs]][which(regulonSelection[[selRegs]] %in% rownames(binaryRegulonActivity))]
        binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
        
        if(nrow(binaryMat)>0) 
        {
          fileName <- paste0(getOutName(scenicOptions, "s4_binaryActivityHeatmap"),selRegs)
          fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
          
          rowv <- ifelse(nrow(binaryMat) >= 2, T, NA)
          colv <- ifelse(ncol(binaryMat) >= 2, T, NA)
          
          NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,   
                        annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                        annColor=colVars,
                        Rowv=rowv,
                        Colv=colv,
                        color = c("white", "black"),
                        filename=fileName)
          if(getSettings(scenicOptions, "devType")!="pdf") dev.off()
        }else{
          if(getSettings(scenicOptions, "verbose")) message(paste0("No regulons to plot for regulon selection '", selRegs, "'. Skipping."))
        }
      }
    }
  }
  
  ################################################################################
  # Tsne - on binary activity
  if(!skipTsne)
  {
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="Binary", filePrefix=getIntName(scenicOptions, "tsne_prefix"), onlyHighConf=FALSE) # default: nPcs, perpl, seed
    if(!is.null(tSNE_fileName))
    {
      tSNE <- readRDS(tSNE_fileName)
      
      # AUCell (activity) as html: 
      fileName <- getOutName(scenicOptions, "s4_binarytSNE_colAct")
      plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally
      
      # Plot cell properties:
      sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
      cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
      colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
      pdf(paste0(getOutName(scenicOptions, "s4_binarytSNE_colProps"),".pdf"))
      plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
      dev.off()
    }
  }
  
  # Finished. Update status.
  scenicOptions@status$current <- 4
  invisible(scenicOptions)
}


################################################################################
# Regulon orders/selection for plots
#' @export
regulonSelections <- function(binaryRegulonActivity, binaryRegulonActivity_nonDupl, minCells)
{
  #binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_full")
  #binaryRegulonActivity_nonDupl <- loadInt(scenicOptions, "aucell_binary_nonDupl")
  
  ### Select regulons:
  regulonSelection <- list(labels=c(all="All regulons \n (including duplicated regulons)",
                                    corr="Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)",
                                    onePercent="Regulons active in more than 1% of cells",
                                    notCorr="Regulons with no other regulons correlated\n abs(cor)>0.30 \n or active in fewer than 1% of cells"))
  
  # All regulons.
  regulonSelection[["all"]] <- rownames(binaryRegulonActivity)
  
  # Active in > 1% cells
  regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
  regulonSelection[["onePercent"]] <- regMinCells
  
  # Correlation across regulons (based on binary cell activity)
  reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
  reguCor[which(is.na(reguCor))] <- 0
  diag(reguCor) <- 0
  
  # Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored.
  corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
  regulonSelection[["corr"]]  <- corrRegs
  missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
  regulonSelection[["notCorr"]]  <- missingRegs
  saveRDS(regulonSelection, file=getIntName(scenicOptions, "aucell_regulonSelection"))
  
  ## Set regulon order (only plotting most correlated regulons)
  reguCor_dist <- as.dist(1-reguCor[corrRegs,corrRegs])
  if(length(reguCor_dist) >= 2) 
  {
    binaryRegulonOrder <- hclust(reguCor_dist)
    binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]
  } else 
  {
    binaryRegulonOrder <- labels(reguCor_dist)
  }
  saveRDS(binaryRegulonOrder, file=getIntName(scenicOptions, "aucell_binaryRegulonOrder"))
  
  return(regulonSelection)
}

runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

####ʶ??ϸ????????????regulons


#########
calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
{
  if(any(is.na(cellAnnotation))) stop("NAs in annotation")
  if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}


plotRSS <- function(rss, labelsToDiscard=NULL, zThreshold=1,
                    cluster_columns=FALSE, order_rows=TRUE, trh=0.01, varName="cellType",
                    col.low="grey90", col.mid="darkolivegreen3", col.high="darkgreen",
                    revCol=FALSE)
{
  varSize="RSS"
  varCol="Z"
  if(revCol) {
    varSize="Z"
    varCol="RSS"
  }
  
  rssNorm <- scale(rss) # scale the full matrix...
  rssNorm[rssNorm < zThreshold] <- 0
  rssNorm <- rssNorm[,which(!colnames(rssNorm) %in% labelsToDiscard)] # remove after calculating...
  
  ## to get topic order (easier...)
  tmp <- .plotRSS_heatmap(rssNorm, trh=trh, cluster_columns=cluster_columns, order_rows=order_rows)
  rowOrder <- rev(tmp@row_names_param$labels)
  
  ## Dotplot
  rss.df <- reshape2::melt(rss)
  head(rss.df)
  colnames(rss.df) <- c("Topic", varName, "RSS")
  rssNorm.df <- reshape2::melt(rssNorm)
  colnames(rssNorm.df) <- c("Topic", varName, "Z")
  rss.df <- merge(rss.df, rssNorm.df)
  
  rss.df <- rss.df[which(rss.df$Z >= 1.5),]
  rss.df <- rss.df[which(!rss.df[,varName] %in% labelsToDiscard),] # remove after calculating...
  # dim(rss.df)
  
  rss.df[,"Topic"] <- factor(rss.df[,"Topic"], levels=rowOrder)
  p <- dotHeatmap(rss.df, 
                  var.x=varName, var.y="Topic", 
                  var.size=varSize, min.size=.5, max.size=5,
                  var.col=varCol, col.low=col.low, col.mid=col.mid, col.high=col.high)
  
  invisible(list(plot=p, df=rss.df, rowOrder=rowOrder))
}

#' @aliases plotRSS
#' @export 
plotRSS_oneSet <- function(rss, setName, n=5)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    geom_point(color = "blue", size = 1) + 
    ggtitle(setName) + 
    geom_label_repel(aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     na.rm=TRUE) +
    theme_classic()
}



## Internal functions:
.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

.plotRSS_heatmap <- plotRSS_heatmap <- function(rss, trh=NULL, row_names_gp=gpar(fontsize=5), order_rows=TRUE, cluster_rows=FALSE, name="RSS", ...)
{
  if(is.null(trh)) trh <- signif(quantile(rss, p=.97),2)
  
  library(ComplexHeatmap)
  rssSubset <- rss[rowSums(rss > trh)>0,]
  rssSubset <- rssSubset[,colSums(rssSubset > trh)>0]
  message("Showing regulons and cell types with any RSS > ", trh, " (dim: ", nrow(rssSubset), "x", ncol(rssSubset),")")
  
  if(order_rows)
  {
    maxVal <- apply(rssSubset, 1, which.max)
    rss_ordered <- rssSubset[0,]
    for(i in 1:ncol(rssSubset))
    {
      tmp <- rssSubset[which(maxVal==i),,drop=F]
      tmp <- tmp[order(tmp[,i], decreasing=FALSE),,drop=F]
      rss_ordered <- rbind(rss_ordered, tmp)
    }
    rssSubset <- rss_ordered
    cluster_rows=FALSE
  }
  
  Heatmap(rssSubset, name=name, row_names_gp=row_names_gp, cluster_rows=cluster_rows, ...)
} 

dotHeatmap <- function (enrichmentDf,
                        var.x="Topic", var.y="ID", 
                        var.col="FC", col.low="dodgerblue", col.mid="floralwhite", col.high="brown1", 
                        var.size="p.adjust", min.size=1, max.size=8,
                        ...)
{
  require(data.table)
  require(ggplot2)
  
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  p <- ggplot(data=enrichmentDf, mapping=aes_string(x=var.x, y=var.y)) + 
    geom_point(mapping=aes_string(size=var.size, color=var.col)) +
    scale_radius(range=c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1),
          ...)
  return(p)
}

# temporary- TODO:delete
#' @export
dotheatmap <- dotHeatmap
#######
#ϸ????????????regulon
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
write.csv(getAUC(regulonAUC),"regulonAUC.csv")
cellAnnotation=cellInfo[colnames(regulonAUC),'CellType']
write.csv(cellAnnotation,"regulonAUC_meta.csv")
rss=calcRSS(AUC=getAUC(regulonAUC),cellAnnotation)
write.csv(rss,"rss.csv")
#
library(ComplexHeatmap)
rssPlot=plotRSS(rss)
pdf("rssPlot.pdf",8,30)
print(rssPlot$plot)
dev.off()
write.csv(rssPlot$df,"rssPlot_df.csv")

##
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
#
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
write.csv(regulonActivity_byCellType_Scaled,"regulonActivity_byCellType_Scaled.csv")
#չʾϸ?????ͺ͵????ӵ?????????ͼ??
pdf("regulonActivity_byCellType_Scaled.pdf",8,30)
print(pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA))
dev.off()

#
pdf("regulonActivity_byCellType.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
write.csv(topRegulators,"RelativeActivity.csv")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

minPerc <- .7  #????????0.7?ȽϺã?˵????70%?õ???????????ϸ????????
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pdf("binaryActPerc_subset1.pdf",10,10)   #????ϸ??????ȡ????intersct
print(pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA))
dev.off()

#????Ȥ????ͼrssPlot_df.csv
pheatmap

###########

setwd("C:/Users/WZ/Documents")
data=read.csv("auc_HEATMAP.csv",row.names = 1)
aov1 <- aov(data[,2]~x, data)  
#a=summary(aov1)
a=TukeyHSD(aov1)[["x"]]


f<-function(x) sum(x>0)
c=aggregate(data[,2],list(data[,1]),f)


for (i in 3:18) {
  aov1 <- aov(data[,i]~x, data)  
  
  b=TukeyHSD(aov1)[["x"]]
  a=cbind(a,b)
  
}
write.csv(a,"ANOVA.csv")

f<-function(x) sum(x>0)
c=aggregate(data[,2],list(data[,1]),f)


for (i in 3:328) {
  d=aggregate(data[,i],list(data[,1]),f)
  c=cbind(c,d)
  
}
write.csv(c,"c.csv")
#
dbLoadingAttempt <- function(dbFilePath){
  ret <- FALSE
  ret <- tryCatch({
    md <- feather::feather_metadata(dbFilePath)
    md$path
    md$dim[2] == length(md$types)
    randomCol <- sample(names(md$types),1)
    rnk <- importRankings(dbFilePath, randomCol)
    TRUE
  }
  , error=function(e){
    print(e$message)
    return(FALSE)
  }
  )
  
  return(ret)
}