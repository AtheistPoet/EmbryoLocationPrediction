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
###dir.random: Directory path where to save MCC and vISH values for the random model
###            This is a subdirectory in the Data directory

##Output: saves MCC matrices for each random sample; also saves matrix of vISH 
##        mismatch for all random samples
Random.null.model <- function( X.RNAseq, distmap.obj, nboot, num.genes, 
                               dir.random){
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
      file <- paste(dir.random, "/MCC_Subchallenge_", j, "_boot_", i, 
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
      file <- paste(dir.random, "/mismatch_VFISH_random.RData", sep = "")
      save(mismatch.vISH.rnd, file = file)
    }
  }
  return(1)
}




##Random.MCC.dist: Calculate average distance of top 10 MCC locations 
##random null model
Random.MCC.dist <- function(X.RNAseq, dir.random, dir.data, num.genes, nboot){
  for( j in 1:length(num.genes) ){
    mcc.dist.rand <- array(0, c(nboot, nrow(X.new)))
    for( b in 1:nboot ){
      cat("Num genes = ", num.genes[j], " boot loop = ", b, "\n")
      dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
      file <- paste(dir.random, "/MCC_Subchallenge_", j, "_boot_", b,
                    ".RData", sep = "")
      load(file)
      id.rand <- apply(MCC.rand, 2, function(x){order(x, decreasing = T)[1:10]})
      mcc.top10 <- apply(MCC.rand, 2, function(x){sort(x, decreasing = T)[1:10]})
      id.maxmcc <- apply(dm@mcc.scores, 2, function(x){order(x, decreasing = T)[1]})
      for( i in 1:nrow(X.new) ){
        # mcc.corr.norm.rand[b, i] <- cor(dm@mcc.scores[, i],
        #                                MCC.mat[, i])
        # id.rand <- order(MCC.mat[, i], decreasing = T)[1:10]
        # id.maxmcc <- order(dm@mcc.scores[, i], decreasing = T)[1]
        if( nnzero(mcc.top10[, i]) != 0 ){
          mcc.dist.rand[b, i] <- log10(1 + sqrt(sum((colSums(dm@geometry[id.rand[, i], ]*
                                                               (cbind(mcc.top10[, i])%*%matrix(1, 1, 3)))/
                                                       sum(mcc.top10[, i]) -
                                                       (dm@geometry[id.maxmcc[i], ]))^2)))
        }else{
          mcc.dist.rand[b, i] <- log10(1 + sqrt(sum((colSums(dm@geometry[id.rand[, i], ]*
                                                               (cbind(rep(1, length(mcc.top10[, i])))%*%matrix(1, 1, 3)))/
                                                       length(mcc.top10[, i]) -
                                                       (dm@geometry[id.maxmcc[i], ]))^2)))
        }
        
      }
      # mcc.corr.norm.rand[b, ] <- corr.maxmcc.rand
      
    }
    file <- paste(dir.data, "/MCC_distance_random_model_Subchallenge_", j, ".RData", 
                  sep = "")
    save(mcc.dist.rand, file = file)
  }
  return(1)
}