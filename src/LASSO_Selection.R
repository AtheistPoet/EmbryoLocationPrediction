#Dependencies: glmnet, foreach, doParallel
#LASSO.Selection: Uses RNAseq data for driver genes across a bunch of cells and returns a 
#ranking of the genes in increasing order of importance for predicting embryonic 
#locations of cells. For each gene i, this builds a LASSO model, where gene i is 
#regressed on all the other genes, and the magnitude of LASSO coefficients 
#determines the significance or strength of each gene in predicting the expression
#pattern of genes i

##Input: X.RNAseq,lambda, nfolds, clsters
###X.RNAseq- single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.
###lambda  - Vector of penalty values for LASSO to select from
###nfolds  - number of folds for cross-validation
###clsters - number of clusters to run the code in parallel

##Output: genes.rank.lasso
###genes.rank.lasso-rank of all the genes as determined by LASSO coefficients and 
###                 averaged over all gene models
LASSO.Selection <- function( X.RNAseq, lambda, nfolds, clsters ){
  
  ##Load required Libraries
  library(glmnet)
  library(foreach)
  library(doParallel)
  
  ##initialize penalty vector, number of cross validation folds
  if( missing(lambda) ){
    lambda <- 10^seq(-5, 5, length = 10)
  }
  if( missing(nfolds) ){
    nfolds <- 5
  }
  if( missing(clsters) ){
    clsters <- 2
  }
  
  ##Register parallel cluster
  cl <- makeCluster(clsters)
  registerDoParallel(cl)
  
  ##Matrix to store LASSO coefficients
  lasso.mat.coef <- matrix(0, ncol(X.RNAseq), ncol(X.RNAseq))
  for( i in 1:ncol(X.RNAseq) ){
    # cat("Working on Gene number = ", i, "\n")
    cv.LASSO.mdl <- cv.glmnet(X.RNAseq[, -i], X.RNAseq[, i], lambda = lambda, 
                        nfolds = nfolds, parallel = F,
                        standardize = F)
    mdl.LASSO.tuned <- glmnet(X.RNAseq[, -i], X.RNAseq[, i], 
                              lambda = cv.LASSO.mdl$lambda.min,
                  standardize = F)
    lasso.mat.coef[-i, i] <- (c(abs(matrix(mdl.LASSO.tuned$beta, ncol = 1))))
    
  }
  stopCluster(cl)
  
  genes.rank.lasso <- rep(0, ncol(X.RNAseq))
  rank.list <- matrix(0, ncol(X.RNAseq), ncol(X.RNAseq))
  for( i in 1:ncol(X.RNAseq) ){
    rank.list[-i, i] <- rank(-lasso.mat.coef[-i, i])
  }
  for( i in 1:ncol(X.new) ){
    genes.rank.lasso[i] <- mean(rank.list[i, -i])
  }
  return(genes.rank.lasso)
}






