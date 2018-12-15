#Loading Data: download the raw RNAseq "dge_raw.txt, normalized RNAseq 
#"dge_normalized.txt", binarized insitu hybridization and geometry data first...
#And them the data on to a DistMap object
####Single-cell genomic
load.data <- function(directory){
  library(DistMap)
  raw.data = read.csv("dge_raw.txt",sep = "\t",header = F)
  rownames(raw.data) = raw.data$V1
  raw.data$V1 = NULL
  
  normalized.data = read.csv("dge_normalized.txt", sep = "\t")
  insitu.matrix = read.csv("binarized_bdtnp.csv",check.names=F)
  
  geometry = read.csv("geometry.txt",sep = " ")
  distmap.obj = new("DistMap",
           raw.data=as.matrix(raw.data),
           data=as.matrix(normalized.data),
           insitu.matrix=as.matrix(insitu.matrix),
           geometry=as.matrix(geometry))
  distmap.obj <- binarizeSingleCellData(distmap.obj, seq(0.15, 0.5, 0.01))
  distmap.obj <- mapCells(distmap.obj)
  return(distmap.obj)
}
