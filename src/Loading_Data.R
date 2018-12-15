#Dependencies: DistMap (can be downloaded from https://github.com/rajewsky-lab/distmap)

#load.data: download the raw RNAseq "dge_raw.txt, normalized RNAseq 
#"dge_normalized.txt", binarized insitu hybridization and geometry data first...
#And them get the data on to a DistMap object
###Input: directory
###   directory-directory path to the data files ("dge_raw.txt", 
###               "dge_normalized.txt", "binarized_bdtnp.csv", "geometry.txt")

###Output: distmap.obj
###  distmap.obj-object that stores the distmap object

load.data <- function(directory){
  library(DistMap)
  
  file.data <- paste(directory, "dge_raw.txt", sep = "")
  raw.data = read.csv(file.data, sep = "\t", header = F)
  rownames(raw.data) = raw.data$V1
  raw.data$V1 = NULL
  
  file.normalized.data <- paste(directory, "dge_normalized.txt", sep = "")
  normalized.data = read.csv(file.normalized.data, sep = "\t")
  
  file.insitu.matrix <- paste(directory, "binarized_bdtnp.csv", sep = "")
  insitu.matrix = read.csv(file.insitu.matrix, check.names=F)
  
  file.geometry <- paste(directory, "geometry.txt", sep = "")
  geometry = read.csv(file.geometry, sep = " ")
  distmap.obj = new("DistMap",
           raw.data=as.matrix(raw.data),
           data=as.matrix(normalized.data),
           insitu.matrix=as.matrix(insitu.matrix),
           geometry=as.matrix(geometry))
  distmap.obj <- binarizeSingleCellData(distmap.obj, seq(0.15, 0.5, 0.01))
  distmap.obj <- mapCells(distmap.obj)
  return(distmap.obj)
}
