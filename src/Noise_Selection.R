#Noise.Selection: Rank genes in increasing order of gene expression noise
##Input: X.RNASeq
### X.RNAseq - single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.

##Output: genes.rank.noise
###genes.rank.noise-rank of all the genes in increasing order of noise


Noise.Selection <- function( X.RNAseq ){
  
  #Compute coefficient of variation (CV) and mean expression levels
  CV <- (apply(X.RNAseq, 2, sd)/colMeans(X.RNAseq))
  mn <- colMeans(X.RNAseq)
  
  ##fit model line to CV vs mean plot
  x <- log10(mn)
  y <- log10(CV)
  mdl <- lm(y ~ x)
  
  ##Rank genes in increasing order of noise
  genes.rank.noise <- rank(y - mdl$fitted.values)
  
  return(genes.rank.noise)
}