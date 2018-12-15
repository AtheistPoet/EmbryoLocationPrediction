#Dependencies: party
#Random.Forest.Selection: Uses RNAseq data for driver genes across a bunch of cells 
#and returns a ranking of the genes in increasing order of importance for predicting
#embryonic locations of cells. For each gene i, this builds a Ramdom Forest 
#regre model, where gene i is regressed on all the other genes, and importance
#is computed based on 'permutation importance' measure. We use the 'cforest' function
#from the party package with unbiased trees, which has been shown to be 
#resilient to selection bias generally exhibited by Random Forests [1]. Further, 
#we add garbage variables (either randomly sampled form the data matrix or 
#gaussian random variables with mean and standard deviation match that of the data
#matrix) to assess the relative performance of the model in selcting relevant features.
#Importance scores for all the genes 
#determine the significance or strength of each gene in predicting the expression
#pattern of gene i
#[1] Strobl, Carolin, et al. "Bias in random forest variable importance measures: 
#Illustrations, sources and a solution." BMC bioinformatics 8.1 (2007): 25.


##Input: X.RNAseq, mtry, ntrees, nboot, num.garbage
### X.RNAseq - single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.
###    mtry  - Number of variables available for splitting at each tree node.
###  ntrees  - Number of trees to grow. 
###   nboot  - number of times to repeat random forest model for each gene. The
###            importance score for each gene model is the average of all nboot models.
### num.garbage - Number of garbage variables to include; half of these would be 
###               randomly sampled from X.RNAseq, the other half are gaussian 
###               random variables with mean = mean(X.RNAseq) and standard deviation = 
###               sd(X.RNAseq)

##Output:  RF.obj
### RF.obj - A list object with two elements: genes.rank.RF and test. if the test variable 
###          is 1, genes.rank.RF is the rank of the genes assigned by Random forest
###         and none of the garbage variables have been selected by the model in 
###         top 60 positions, otherwise garbage variables have been selected
Random.Forest.Selection <- function( X.RNAseq, mtry, ntrees, nboot, 
                                      num.garbage){
  library(party)
  set.seed(1000089)
  library(party)
  if( missing(ntrees) ){
    ntrees <- 1000
  }
  if( missing(mtry) ){
    mtry <- 10
  }
  if( missing(nboot) ){
    nboot <- 100
  }
  if( missing(num.dummy) ){
    num.garbage <- 10
  }
  rf.importance.mat <- matrix(0, ncol(X.RNAseq), ncol(X.RNAseq))
  
  varimp.list <- matrix(0, ncol(X.RNAseq), ncol(X.RNAseq) + num.garbage)
  for( i in 1:ncol(X.RNAseq) ){
    for( b in 1:nboot ){
      cat("Running Gene loop = ", i, " and boot loop = ", b, "\n")
      rf.dat.mat <- cbind(X.RNAseq, matrix(X.RNAseq[sample(1:length(c(X.RNAseq)), 
                                                    (num.garbage/2)*nrow(X.RNAseq), 
                                              replace = T)], nrow(X.RNAseq)))
      rf.dat.mat <- cbind(rf.dat.mat, matrix(rnorm((num.garbage/2)*nrow(X.RNAseq), 
                                            mean = mean(c(X.RNAseq)), 
                                     sd(c(X.RNAseq))), nrow(X.RNAseq)))
      colnames(rf.dat.mat)[(ncol(X.RNAseq) + 1):(ncol(X.RNAseq) + num.garbage)] <- 
        1:num.garbage
      rf.dat.mat <- data.frame(rf.dat.mat)
      eval(substitute(rf.mdl <- cforest(variable ~ ., data = rf.dat.mat,
                                     control = cforest_unbiased(ntree = ntrees, 
                                                                mtry = mtry)), 
                      list(variable = as.name(colnames(rf.dat.mat)[i]))))
      varimp.list[i, -i] <- varimp.list[i, -i] + varimp(rf.mdl)
    }
    varimp.list[i, -i] <- varimp.list[i, -i]/nboot
    save(varimp.list, file = "RF_variable_importance_84_genes.RData")
  }
  
  genes.rank.RF <- rep(0, ncol(X.RNAseq) + num.garbage)
  genes.rank.RF.list <- matrix(0, ncol(X.RNAseq), ncol(X.RNAseq) + num.garbage)
  for( i in 1:ncol(X.RNAseq) ){
    genes.rank.RF.list[i, ] <- rank(-varimp.list[1, ])
  }
  genes.rank.RF <- colMeans(genes.rank.RF.list)
  id.order <- order(genes.rank.RF, decreasing = F)
  if( length(intersect(id[1:60], 
                       (ncol(X.RNAseq) + 1):(ncol(X.RNAseq) + num.garbage))) == 0 ){
    test <- 1
  }else{
    test <- 0
  }
  RF.obj <- list()
  RF.obj$genes.rank.RF <- genes.rank.RF
  RF.obj$test <- test
  return(RF.obj)
}