library(MASS)

### Function to estimate prediction matrix
## Param:
# norm_count_for_estim: named normalized count matrix (genes on rows, cells on columns) used for estimation
# known_genes: gene names of genes to be used in the prediction
# predict_genes: gene names of genes to predict
## Return:
# matrix (of #predicted genes X #cells elements) containg the element for the prediction
estimate_func <- function(norm_count_for_estim, known_genes, predict_genes){
  
  
  data_estimation <- norm_count_for_estim
  
  # list of known genes
  index_X <- known_genes
  N_known <- length(index_X) # number of known genes (60, 40 or 20)
  
  
  # list of predicted genes
  index_Y <- predict_genes
  N_est <- length(index_Y)
  
  
  X_estimation <- data_estimation[index_X,] # matrix of N_known * 1297 elements
  Y_estimation <- data_estimation[index_Y,] # matrix of N_est * 1297 elements
  
  # calculate this term only one time
  m_inv <- ginv(X_estimation %*% t(X_estimation)) %*% X_estimation 
  

  
  A <- matrix(0, nrow = N_est, ncol = N_known)
  for (i in 1:N_est){
    
    y <- Y_estimation[i,]
    
    a_LS <- m_inv %*% y
    
    A[i, ] <- a_LS
    
  }
  
  # assign row anc col names
  rownames(A) <- predict_genes
  colnames(A) <- known_genes
  
  return (A)
  
}


### Function to predict
## Param:
# norm_count_for_predict: named normalized count matrix (genes on rows, cells on columns) used for prediction
# estimated_matrix: matrix (of #predicted genes X #cells elements) to use in the prediction process
# known_genes: gene names of genes to be used in the prediction
# predict_genes: gene names of genes to predict
## Return:
# named matrix (of #predict_genes X #cell element) contaning the prediction of gene expression values for genes in predict_genes
prediction_func  <- function(norm_count_for_predict, estimated_matrix, known_genes, predict_genes ){
  

  X_prediction = norm_count_for_predict # matrix of N_known * 3039 elements
  
  Y_predicted <- matrix(nrow = length(predict_genes), ncol = ncol(norm_count_for_predict))
  rownames(Y_predicted) <- rownames(estimated_matrix)
  colnames(Y_predicted) <- colnames(norm_count_for_predict)
  
  for (j in c(1:ncol(norm_count_for_predict))){ #for i in 1:3039
    
    Y_predicted[,j] <- estimated_matrix %*% X_prediction[,j]
  }
  
  return(Y_predicted)
}


### Function to predict the expression of some genes, given the expression level of other genes
## Param:
# norm_count_for_estim: named normalized count matrix (genes on rows, cells on columns) used for estimation
# norm_count_for_predict: named normalized count matrix (genes on rows, cells on columns) used for prediction
# known_genes: gene names of genes to be used in the prediction
# predict_genes: gene names of genes to predict
## Return:
# named matrix (of #predict_genes X #cell element) contaning the prediction of gene expression values for genes in predict_genes
predict_genes <- function(norm_count_for_estim, norm_count_for_predict, known_genes, predict_genes){
  
  est_mat <- estimate_func(norm_count_for_estim=norm_count_for_estim, known_genes=known_genes, predict_genes=predict_genes)
  
  pred_mat <- prediction_func(norm_count_for_predict=norm_count_for_predict, estimated_matrix=est_mat, predict_genes=predict_genes )
  
  
  return(pred_mat)
  
}


