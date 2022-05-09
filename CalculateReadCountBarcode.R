CalculateReadCountBarcode <- function(mol.info.file,  row.prefix){
  
  library(Seurat)
  library(DropletUtils)
  molinfo <- read10xMolInfo(mol.info.file)
  selected_info <- molinfo$data[1:5]
  selected_info$cell <- factor(selected_info$cell)
  cell_read_sum <- tapply(selected_info$reads,selected_info$cell,sum)
  
  cell_read_sum <- as.matrix(cell_read_sum)
  rownames(cell_read_sum) <- paste0(row.prefix,  rownames(cell_read_sum))
  colnames(cell_read_sum) <- "reads_sum"
  print(dim(cell_read_sum)) 
  print(head(cell_read_sum))
  return(cell_read_sum)
}