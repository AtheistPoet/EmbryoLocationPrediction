##PCA.Selection: Ranks genes by factor loadings along principal components.

###Input: X.RNAseq, num.PC
### X.RNAseq-single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.
### num.PC-Number of top principal components to use for ranking genes according to 
###        factor loadings

###Output: genes.rank.PCA
###  genes.rank.PCA- Average rank of genes 
PCA.Selection <- function( X.RNAseq, num.PC ){
  if( missing(num.PC) ){
    num.PC <- 3
  }
  PCA.obj <- prcomp(X.RNAseq)
  rank.PC <- matrix(0, num.PC, ncol(X.new))
  for( i in 1:num.PC ){
    rank.PC[i, ] <- rank(-abs(mdl.PCA$rotation[, i]))
  }
  genes.rank.PCA <- colMeans(rank.PC)
  return(genes.rank.PCA)
}