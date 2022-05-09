#!/home/zhengjh/miniconda3/bin/Rscript

# Function: This is a pipeline to find DEGs using DESeq2
# Input: 1.rds odject which includes counts variable and label variable; 2. outputname to save the DEGs
# Output: csv of DESeq2 DEGs
# Version: 2021-03-05, Jihong Zheng

library(optparse)
library(getopt)

option_list <- list(
  
  make_option(c("-o", "--output"), type = "character", default = FALSE,
              action = "store", help = "This is the output directory."
  ),
  make_option(c("--rds_name"), type = "character", default = FALSE,
              action = "store", help = "This is the name of rds to input."
  ),
  make_option(c("--output_name"), type = "character", default = FALSE,
              action = "store", help = "This is the name of ouput file."
  )
)

# -help 
opt = parse_args(OptionParser(option_list = option_list, 
                              usage = "This is a script to find DEGs using DESeq2"))
print(opt)

# set the output path
library(DESeq2)
DESeq2_Object <- readRDS(opt$rds_name)

setwd(opt$output)
timestart<-Sys.time() 
print("The time begining:")
print(timestart)

##  Main function 
print("The level and number of input data:")
# The levels of condition need to assign by yourself. DESeq2 will use the first level as the condtion.
#condition<-factor(c(rep("Control",7104), rep("Obesity",6965)), levels = c("Control","Obesity"))

condition <- factor(DESeq2_Object$label)
print(levels(condition))
print(table(condition))

print("The contrl is")
print(table(condition)[1])
print("The treat is")
print(table(condition)[2])


count <- DESeq2_Object$counts
count <- as.matrix(count) # Need to transform to matrix.

coldata <- data.frame(row.names=colnames(count), condition)                    
count_2 <- count[which(rowSums(count) > 0), ] # Exclude the gene which is 0 in all cells.
count <-  count_2+1
dds <- DESeqDataSetFromMatrix(count,coldata,design=~condition)
dds <- DESeq(dds)
res <- results(dds)
outputname <- opt$output_name

print("Finished to calculate DEGs and now writing the DEGs results.")
write.csv(res, outputname)

timeend<-Sys.time()
print("The time ending:")
print(timeend)
runningtime<-timeend-timestart
print("The time DESeq2 used is ")
print(runningtime) 

