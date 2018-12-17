#Denoising.Autoencoder.Selection
#Dependencies: Keras for R; Tensorflow for R
##Uncomment below to install keras for R
#library(keras)
#install_keras()   ###Installs tensorflow backend as well

##Denoising.Autoencoder.Selection: Ranks gene by training denoising autoencoder models
##for all the genes. For gene i, we train an autoencoder, while keeping the input for
##gene i either clamped to zero or randomly shuffling the input. The error on
##test set then gives an indication of the importance of gene i

##Inputs: X.RNaseq, train.ratio, val.ratio, number_layers, hidden_neurons, dropout, 
##        learn_rate, corrupt, batch_size, batch_size_val, epochs, verbose, 
##        activation, Type

###X.RNAseq - single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.
##train.ratio - percentage of data to use for training
##val.ratio - percentage of data to use for validation
##number_layers - number of hidden layers to use. If this is greater than one, then 
##                the effective number of layers is twice this number (the encoder and 
##                decoder are symmetrical.)
##hidden_neurons - number of neurons in each hidden layer.
##dropout - dropout rate for regularization.
##learn_rate - learning rate
## corrupt - corruption level for input genes.
## batch_size - training batch size.
##batch_size_val - validation batch size.
##epochs - number of epochs to train.
##verbose - either '0' or '1', Whether to print progress or not.
##activation - activation function for th hidden layers
##Type - Type of gene importance technique. If this is 'zero', then for each gene,
##       its neural network has the expression for gene i clamped to zero, and if this 
##       is set to 'random', gene i's expression is randomly shuffled.

##Output: genes.rank.DA (gene ranks)

Denoising.Autoencoder.Selection <- function( X.RNAseq, train.ratio, 
                                             val.ratio, number_layers, 
                                             hidden_neurons, 
                                             dropout, learn_rate, corrupt, 
                                             batch_size, batch_size_val, 
                                             epochs, verbose, 
                                             activation, 
                                             Type){
  set.seed(1)
  if( missing(train.ratio) ){
    train.ratio <- 0.7
  }
  if( missing(val.ratio) ){
    val.ratio <- 0.2
  }
  if( missing(number_layers) ){
    number_layers <- 1
  }
  if( missing(hidden_neurons) ){
    hidden_neurons <- 256
  }
  if( missing(dropout) ){
    dropout <- 0.1
  }
  if( missing(learn_rate) ){
    learn_rate <- 0.001
  }
  if( missing(corrupt) ){
    corrupt <- 0.5
  }
  if( missing(batch_size) ){
    batch_size <- 50
  }
  if( missing(batch_size_val) ){
    batch_size_val <- 10
  }
  if( missing(epochs) ){
    epochs <- 1000
  }
  if( missing(verbose) ){
    verbose <- 1
  }
  if( missing(Type) ){
    Type <- "zero"
  }
  if( missing(activation) ){
    activation <- "tanh"
  }
  test.ratio <- 1 - train.ratio - val.ratio
  
  
  ##Generate training, validation and test data
  id.reorder <- sample(1:nrow(X.RNAseq), nrow(X.RNAseq), replace = F)
  id.train <- id.reorder[1:floor(train.ratio*length(id.reorder))]
  id.test.val <- setdiff(id.reorder, id.train)
  id.val <- id.test.val[sample(1:length(id.test.val), floor(val.ratio*length(id.reorder)), 
                               replace = F)]
  id.test <- setdiff(id.test.val, id.val)
  
  mn.train <- colMeans(X.RNAseq[id.train, ])
  sd.train <- apply(X.RNAseq[id.train, ], 2, sd)
  X.train <- (X.RNAseq[id.train, ] - matrix(1, length(id.train), 1)%*%mn.train)/
    matrix(1, length(id.train), 1)%*%sd.train
  X.train.full <- X.train
  X.test <- (X.RNAseq[id.test, ] - matrix(1, length(id.test), 1)%*%mn.train)/
    matrix(1, length(id.test), 1)%*%sd.train
  X.test.full <- X.test
  X.val <- (X.RNAseq[id.val, ] - matrix(1, length(id.val), 1)%*%mn.train)/
    matrix(1, length(id.val), 1)%*%sd.train
  X.val.full <- X.val
  
  ##Generator functions for corrupting data during training and validation of the 
  ##neural networks
  denoising_generator <- function(X_data, Y_data, batch_size, corrupt) {
    function() {
      rows <- sample(1:nrow(X_data), batch_size, replace = F)
      id <- sample(1:ncol(X_data), floor(corrupt*ncol(X_data)), replace = F)
      X_train <- X_data[rows, ]
      Y_train <- Y_data[rows, ]
      X_train[, id] <- 0
      list(X_train, Y_train)
    }
  }
  
  validation_generator <- function(X_data, Y_data, batch_size, corrupt){
    function() {
      rows <- sample(1:nrow(X_data), batch_size, replace = F)
      id <- sample(1:ncol(X_data), floor(corrupt*ncol(X_data)), replace = F)
      X_train <- X_data[rows, ]
      Y_train <- Y_data[rows, ]
      X_train[, id] <- 0
      list(X_train, Y_train)
    }
  }
  
  
  ##Define network architecture
  if( number_layers == 1 ){
    model <- keras_model_sequential()
    model %>%
      layer_dense(units = hidden_neurons[1], activation = activation, 
                  input_shape = ncol(X.RNAseq)) %>%
      layer_dropout(rate = dropout) %>% 
      layer_dense(units = ncol(X.RNAseq))
  }else{
    if( number_layers == 2 ){
      model <- keras_model_sequential()
      model %>%
        layer_dense(units = hidden_neurons[1], activation = activation, 
                    input_shape = ncol(X.RNAseq)) %>%
        layer_dropout(rate = dropout) %>% 
        layer_dense(units = hidden_neurons[2], activation = activation) %>%
        layer_dropout(rate = dropout) %>%
        layer_dense(units = hidden_neurons[1], activation = activation) %>%
        layer_dropout(rate = dropout) %>%
        layer_dense(units = ncol(X.RNAseq))
    }else{
      model <- keras_model_sequential()
      model %>%
        layer_dense(units = hidden_neurons[1], activation = activation, 
                    input_shape = ncol(X.RNAseq)) %>%
        layer_dropout(rate = dropout) %>% 
        layer_dense(units = hidden_neurons[2], activation = activation) %>%
        layer_dropout(rate = dropout) %>%
        layer_dense(units = hidden_neurons[3], activation = activation) %>%
        layer_dropout(rate = dropout) %>%
        layer_dense(units = hidden_neurons[2], activation = activation) %>%
        layer_dropout(rate = dropout) %>%
        layer_dense(units = hidden_neurons[1], activation = activation) %>%
        layer_dropout(rate = dropout) %>%
        layer_dense(units = ncol(X.RNAseq))
    }
  }
  
  
  summary(model)
  
  model %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam", 
    metrics = list("mean_absolute_error")
  )
  k_set_value(model$optimizer$lr, learn_rate)
  
  print_dot_callback <- callback_lambda(
    on_epoch_end = function(epoch, logs) {
      if (epoch %% 20 == 0) cat("\n")
      cat(".")
    }
  )    
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = 10)
  Mae.gene <- rep(0, ncol(X.RNAseq))
  for( i in 1:ncol(X.RNAseq) ){
    if( Type == "zero" ){
      X.train <- X.train.full
      X.train[, i] <- 0
      X.test <- X.test.full
      X.test[, i] <- 0
      X.val <- X.val.full
      X.val[, i] <- 0
    }else{
      X.train <- X.train.full
      X.train[, i] <- X.train.full[sample(1:nrow(X.train.full, 
                                                 nrow(X.train.full))), i]
      X.test <- X.test.full
      X.test[, i] <- X.test.full[sample(1:nrow(X.test.full, 
                                                nrow(X.test.full))), i]
      X.val <- X.val.full
      X.val[, i] <- X.val.full[sample(1:nrow(X.val.full, 
                                              nrow(X.val.full))), i]
    }
    
    model.gene <- model
    history <- model.gene %>% fit_generator(denoising_generator(X.train, 
                                                           X.train.full, 
                                                           batch_size, corrupt),
                                       epochs = epochs,
                                       steps_per_epoch = floor(nrow(X.train)/batch_size),
                                       validation_data = validation_generator(X.val, 
                                                                              X.val.full, 
                                                                              batch_size_val, 
                                                                              corrupt),
                                       validation_steps = floor(nrow(X.val)/batch_size_val),
                                       verbose = verbose,
                                       callbacks = list(print_dot_callback)
    )
    c(loss, mae) %<-% (model.gene %>% evaluate(X.test, X.test.full, 
                                          verbose = verbose))
    Mae.gene[i] <- mae
  }
  genes.rank.DA <- rank(Mae.gene)
  
  return(genes.rank.DA)
}


#Denoising.Autoencoder.Selection.TF
##Denoising.Autoencoder.Selection.TF: Ranks gene by training denoising autoencoder models
##for all the genes. For gene i, we train an autoencoder, while keeping the input for
##gene i either clamped to zero or randomly shuffling the input. The error on
##test set then gives an indication of the importance of gene i. Prior 
##knowledge about transcriptional regulatory network is also included. 
##The hidden layer also has some TF nodes

##Inputs: X.RNaseq, train.ratio, val.ratio, hidden_neurons, dropout, 
##        learn_rate, corrupt, batch_size, batch_size_val, epochs, verbose, 
##        activation, Type

###X.RNAseq - single cell RNAseq data. This is a matrix with cells on the rows and 
###          genes on teh columns.
##train.ratio - percentage of data to use for training
##val.ratio - percentage of data to use for validation
##hidden_neurons - number of neurons in each hidden layer.
##dropout - dropout rate for regularization.
##learn_rate - learning rate
## corrupt - corruption level for input genes.
## batch_size - training batch size.
##batch_size_val - validation batch size.
##epochs - number of epochs to train.
##verbose - either '0' or '1', Whether to print progress or not.
##activation - activation function for th hidden layers
##Type - Type of gene importance technique. If this is 'zero', then for each gene,
##       its neural network has the expression for gene i clamped to zero, and if this 
##       is set to 'random', gene i's expression is randomly shuffled.
##adj.TRN - adjacency matrix for transcriptional regulatory network 
##dir.data - path to data directory which has the file for name conversion between official gene symbol and 
##         flybase ids
##id.TF - position of driver genes in the adjacency matrix

##Output: genes.rank.DA.TF (gene ranks)

Denoising.Autoencoder.Selection.TF <- function( X.RNAseq, train.ratio, 
                                             val.ratio, 
                                             dropout, corrupt, 
                                             batch_size, batch_size_val, 
                                             epochs, verbose, 
                                             activation, 
                                             Type, adj.TRN, dir.data){
  library(keras)
  if( missing(train.ratio) ){
    train.ratio <- 0.7
  }
  if( missing(val.ratio) ){
    val.ratio <- 0.2
  }
  if( missing(dropout) ){
    dropout <- 0.1
  }
  if( missing(corrupt) ){
    corrupt <- 0.5
  }
  if( missing(batch_size) ){
    batch_size <- 50
  }
  if( missing(batch_size_val) ){
    batch_size_val <- 10
  }
  if( missing(epochs) ){
    epochs <- 1000
  }
  if( missing(verbose) ){
    verbose <- 1
  }
  if( missing(Type) ){
    Type <- "zero"
  }
  if( missing(activation) ){
    activation <- "tanh"
  }
  test.ratio <- 1 - train.ratio - val.ratio
  
  # dir.TRN <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.data, "/DAVID_gene_DREAM_84_name_to_flybase.txt", sep = "")
  gene.to.flybase <- read.table(file, sep = "\t", header = T)
  
  names.flybase <- gene.to.flybase[which(gene.to.flybase$From != "D" & 
                                           gene.to.flybase$From != "h" & 
                                           gene.to.flybase$From != "trn"), ]
  names.flybase <- rbind(names.flybase, gene.to.flybase[
    which(gene.to.flybase == "D")[1], ], gene.to.flybase[
      which(gene.to.flybase == "h")[2], ], 
    gene.to.flybase[
      which(gene.to.flybase == "trn")[2], ])
  
  id.DREAM.fly <- match(colnames(X.RNAseq), names.flybase$From)
  DREAM.flybase <- as.character(names.flybase$To[id.DREAM.fly[!is.na(id.DREAM.fly)]])
  
  id.match <- match(DREAM.flybase, colnames(adj.TRN))
  
  adj.DA <- adj.TRN[id.match[!is.na(id.match)], 
                    id.match[!is.na(id.match)]]
  
  id.TF <- which(rowSums(adj.DA) != 0)
  
  
  
  ##Generate training, validation and test data
  id.reorder <- sample(1:nrow(X.RNAseq), nrow(X.RNAseq), replace = F)
  id.train <- id.reorder[1:floor(train.ratio*length(id.reorder))]
  id.test.val <- setdiff(id.reorder, id.train)
  id.val <- id.test.val[sample(1:length(id.test.val), floor(val.ratio*length(id.reorder)), 
                               replace = F)]
  id.test <- setdiff(id.test.val, id.val)
  
  mn.train <- colMeans(X.RNAseq[id.train, ])
  sd.train <- apply(X.RNAseq[id.train, ], 2, sd)
  X.train <- (X.RNAseq[id.train, ] - matrix(1, length(id.train), 1)%*%mn.train)/
    matrix(1, length(id.train), 1)%*%sd.train
  X.train.full <- X.train
  X.test <- (X.RNAseq[id.test, ] - matrix(1, length(id.test), 1)%*%mn.train)/
    matrix(1, length(id.test), 1)%*%sd.train
  X.test.full <- X.test
  X.val <- (X.RNAseq[id.val, ] - matrix(1, length(id.val), 1)%*%mn.train)/
    matrix(1, length(id.val), 1)%*%sd.train
  X.val.full <- X.val
  
  ##Training, validation and test data for network with TF nodes
  X.train <- X.test <- X.val <- list()
  for( i in 1:length(id.TF) ){
    X.train[[i]] <- cbind(X.train.full[, which(adj.DA[id.TF[i], ] != 0)])
    X.test[[i]] <- cbind(X.test.full[, which(adj.DA[id.TF[i], ] != 0)])
    X.val[[i]] <- cbind(X.val.full[, which(adj.DA[id.TF[i], ] != 0)])
  }
  X.train[[i + 1]] <- X.train.full
  X.test[[i + 1]] <- X.test.full
  X.val[[i + 1]] <- X.val.full
  
  ##Generator functions for corrupting data during training and validation of the 
  ##neural networks
  denoising_generator <- function(X_data, Y_data, batch_size, corrupt, 
                                  adj.DA, id.TF) {
    function() {
      X_train <- list()
      rows <- sample(1:nrow(X_data[[length(X_data)]]), batch_size, replace = F)
      id <- sample(1:ncol(X_data[[length(X_data)]]), 
                   floor(corrupt*ncol(X_data[[length(X_data)]])), replace = F)
      for( i in 1:(length(X_data) - 1) ){
        X_train[[i]] <- cbind(X_data[[i]][rows, ])
        id.match <- match(id, which(adj.DA[id.TF[i], ] != 0))
        if( length(id.match[!is.na(id.match)]) != 0 ){
          X_train[[i]][, id.match[!is.na(id.match)]] <- 0
        }
      }
      X_train[[i + 1]] <- X_data[[i + 1]][rows, ]
      X_train[[i + 1]][, id] <- 0
      Y_train <- Y_data[rows, ]
      list(X_train, Y_train)
    }
  }
  
  validation_generator <- function(X_data, Y_data, batch_size, corrupt, 
                                   adj.DA, id.TF){
    function() {
      X_train <- list()
      rows <- sample(1:nrow(X_data[[length(X_data)]]), batch_size, replace = F)
      id <- sample(1:ncol(X_data[[length(X_data)]]), 
                   floor(corrupt*ncol(X_data[[length(X_data)]])), replace = F)
      for( i in 1:(length(X_data) - 1) ){
        X_train[[i]] <- cbind(X_data[[i]][rows, ])
        id.match <- match(id, which(adj.DA[id.TF[i], ] != 0))
        if( length(id.match[!is.na(id.match)]) != 0 ){
          X_train[[i]][, id.match[!is.na(id.match)]] <- 0
        }
      }
      X_train[[i + 1]] <- X_data[[i + 1]][rows, ]
      X_train[[i + 1]][, id] <- 0
      Y_train <- Y_data[rows, ]
      list(X_train, Y_train)
    }
  }
  
  
  
  ##Define network architecture
  ##model
  input <- list()
  for( i in 1:length(id.TF) ){
    cat("Loop = ", i, sep = "\n")
    input[[i]] <- layer_input(shape = c(nnzero(adj.DA[id.TF[i], ])), 
                              name = paste("input_", i, sep = ""))
    
  }
  input[[length(id.TF) + 1]] <- layer_input(shape = ncol(adj.DA), 
                                            name = paste("input_", i + 1, sep = ""))
  output_hidden <- list()
  for( i in 1:length(id.TF) ){
    output_hidden[[i]] <- input[[i]] %>%
      layer_dense(units = 1, activation = activation, 
                  kernel_regularizer = regularizer_l2(l = 0.001), 
                  name = paste("output_", i, sep = ""))
  }
  output_hidden[[i + 1]] <- input[[i + 1]] %>%
    layer_dense(units = 234, activation = activation, 
                kernel_regularizer = regularizer_l2(l = 0.001), 
                name = paste("output_", i + 1, sep = ""))
  output <- layer_concatenate(output_hidden) %>% layer_dense(units = ncol(X.RNAseq))
  
  for( i in 1:length(input) ){
    if( i == 1 ){
      input.list <- c(input[[i]])
    }else{
      input.list <- c(input.list, input[[i]])
    }
  }
  model <- keras_model(
    inputs = c(input.list),
    outputs = output
  )
  summary(model)
  
  model %>% compile(
    loss = "mean_squared_error", 
    optimizer = "adam", 
    metrics = list("mean_absolute_error")
  )
  
  print_dot_callback <- callback_lambda(
    on_epoch_end = function(epoch, logs) {
      if (epoch %% 20 == 0) cat("\n")
      cat(".")
    }
  )    
  
  early_stop <- callback_early_stopping(monitor = "val_loss", patience = 10)
  Mae.gene <- rep(0, ncol(X.RNAseq))
  for( i in 1:ncol(X.RNAseq) ){
    if( Type == "zero" ){
      for( j in 1:(length(X.train) - 1) ){
        id.match <- match(i, which(adj.DA[id.TF[j], ] != 0))
        if( length(id.match[!is.na(id.match)]) != 0 ){
          X.train[[j]][, id.match[!is.na(id.match)]] <- 0
          X.test[[j]][, id.match[!is.na(id.match)]] <- 0
          X.val[[j]][, id.match[!is.na(id.match)]] <- 0
        }
      }
      X.train[[j + 1]][, i] <- 0
      X.test[[j + 1]][, i] <- 0
      X.val[[j + 1]][, i] <- 0
    }else{
      for( j in 1:(length(X.train) - 1) ){
        id.match <- match(i, which(adj.DA[id.TF[j], ] != 0))
        if( length(id.match[!is.na(id.match)]) != 0 ){
          X.train[[j]][, id.match[!is.na(id.match)]] <- 
            sample(X.train[[j]][, id.match[!is.na(id.match)]], 
                   lenth(X.train[[j]][, id.match[!is.na(id.match)]]))
          X.test[[j]][, id.match[!is.na(id.match)]] <- 
            sample(X.test[[j]][, id.match[!is.na(id.match)]], 
                   length(X.test[[j]][, id.match[!is.na(id.match)]]))
          X.val[[j]][, id.match[!is.na(id.match)]] <- 
            sample(X.val[[j]][, id.match[!is.na(id.match)]], 
                   length(X.val[[j]][, id.match[!is.na(id.match)]]))
        }
      }
      X.train[[j + 1]][, i] <- 
        sample(X.train[[j + 1]][, i], length(X.train[[j + 1]][, i]))
      X.test[[j + 1]][, i] <- 
        sample(X.test[[j + 1]][, i], length(X.test[[j + 1]][, i]))
      X.val[[j + 1]][, i] <- 
        sample(X.val[[j + 1]][, i], length(X.val[[j + 1]][, i]))
    }
    
    model.gene <- model
    
    history <- model.gene %>% fit_generator(denoising_generator(X.train, 
                                                           X.train.full, 
                                                           batch_size, corrupt, 
                                                           adj.DA, id.TF),
                                       epochs = epochs,
                                       steps_per_epoch = floor(nrow(X.train[[1]])/batch_size),
                                       validation_data = validation_generator(X.val, X.val.full, 
                                                                              batch_size_val, 
                                                                              corrupt, 
                                                                              adj.DA, id.TF),
                                       validation_steps = floor(nrow(X.val[[1]])/batch_size_val),
                                       verbose = verbose,
                                       callbacks = list(print_dot_callback)
    )
    
    
    c(loss, mae) %<-% (model.gene %>% evaluate(X.test, X.test.full, 
                                          verbose = verbose))
    Mae.gene[i] <- mae
  }
  genes.rank.DA.TF <- rank(Mae.gene)
  
  return(genes.rank.DA.TF)
}


