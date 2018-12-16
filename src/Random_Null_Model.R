#Random null model for performance evaluation
#Dependencies: DistMap (can be downloaded from https://github.com/rajewsky-lab/distmap)
##Random.null.model - randomly sample gene subsets and predict cellular positions 
##                    by computing matthews correlation coefficient for every 
##                    cell, location pair.

##Input: X.RNAseq, distamp.obj, nboot, num.genes
###X.RNAseq- single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.
###distmap.obj - DistMap (R package) object that contains expression values, and
###              MCC values.
###nboot - number of random samples generated.
###num.genes - vector of number of genes for each sub-challenge on the 
###            DREAM Single Cell Trnascriptomic Challenge
###dir.MCC.rand: Directory where to save MCC values
###dir.vISH.rand: Directory where to save virtual InSitu Hybridization values

##Output: saves MCC matrices for each random sample; also saves matrix of vISH 
##        mismatch for all random samples
Random.null.model <- function( X.RNAseq, distmap.obj, nboot, num.genes, 
                               dir.MCC.rand, 
                               dir.vISH.rand){
  library(DistMap)
  if( missing(nboot) ){
    nboot <- 100
  }
  mismatch.vISH.rnd <- array(0, c(length(num.genes), nboot, nrow(dm@geometry)))
  for( j in 1:(length(num.genes)) ){
    for( i in 1:nboot ){
      genes <- sample(colnames(X.RNAseq), num.genes[j], replace = F)
      id.match <- match(colnames(X.RNAseq), as.character(unlist(genes[1:num.genes[j]])))
      dm.tmp = new("DistMap",
                   raw.data=as.matrix(distmap.obj@raw.data),
                   data=as.matrix(distmap.obj@data),
                   insitu.matrix=distmap.obj@insitu.matrix[, !is.na(id.match)],
                   geometry=distmap.obj@geometry)
      dm.tmp <- binarizeSingleCellData(dm.tmp, seq(0.15, 0.5, 0.01))
      dm.tmp <- mapCells(dm.tmp)
      dm.tmp@mcc.scores[is.nan(dm.tmp@mcc.scores)] <- 0
      # dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
      file <- paste(dir.MCC.rand, "/Random_MCC/MCC_Subchallenge_", j, "_boot_", i, 
                    ".RData", sep = "")
      MCC.rand <- dm.tmp@mcc.scores
      save(MCC.rand, file = file)
      mismatch.tmp <- matrix(0, ncol(X.RNAseq), nrow(dm@geometry))
      for( k in 1:ncol(X.RNAseq) ){
        cat("Num genes = ", num.genes[j], " nboot = ", i,
            " Gene = ", k, "\n")
        pha <- computeVISH(dm, colnames(X.RNAseq)[k], threshold = 0.75)
        pha.tmp <- computeVISH(dm.tmp, colnames(X.RNAseq)[k], threshold = 0.75)
        
        VISH.diff <- abs(pha - pha.tmp)/max(pha)
        id.nnz <- which(VISH.diff >= 1e-1)
        mismatch.tmp[k, id.nnz] <- 1
      }
      mismatch.vISH.rnd[j, i, ] <- colMeans(mismatch.tmp)
      # dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
      file <- paste(dir.vISH.rand, "/mismatch_VFISH_random.RData", sep = "")
      save(mismatch.vISH.rnd, file = file)
    }
  }
  
}