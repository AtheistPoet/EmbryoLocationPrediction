##Test

##current directory 
dir.cur <- getwd()

##Name of Method files
Method <- c("Lasso", "PCA", "Noise", "RF", "Denoising_Autoencoder", 
                      "DA_one_Layer", "DA_two_Layer", 
                      "DA_three_Layer", "Denoising_Autoencoder_TF", "scPROACTIVE",
            "MLB", "All Genes Top 10")

##Name of the Method in the paper
Method_Name <- c("Lasso", "PCA", "Noise", "RF", "DA", 
                 "DA1", "DA2", "DA3", "DA TF", "scPROACTIVE", "MLB", 
                 "All Genes Top 10")



##loading data
###dir.data <- (add data directory path here)
dm <- load.data(dir.data, RNAseq.full = F)

###RNAseq data
X.RNAseq <- t(dm@data)

i <- 1
##dir.Methods <- (directory for storing gene rankings from Methods)
##Uncomment to run lasso and get gene rankings
###genes.rank.lasso <- LASSO.Selection(X.RNAseq, lambda, nfolds, clsters)
# id <- order(genes.rank.lasso, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.lasso[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1

##Uncomment to run PCA selection and get gene rankings
#num.PC <- 3
# genes.rank.PCA <- PCA.Selection(X.RNAseq, num.PC)
# id <- order(genes.rank.PCA, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.PCA[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1

##Uncomment to run Noise selection and get gene rankings
# genes.rank.noise <- Noise.Selection(X.RNAseq)
# id <- order(genes.rank.noise, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.noise[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1


##Uncomment to run Random Forest Selection and get gene rankings
# mtry <- 10
# ntrees <- 1000
# nboot <- 100
# num.garbage <- 10
# genes.rank.RF <- Random.Forest.Selection(X.RNAseq, mtry, ntrees, nboot, 
#                                             num.garbage)
# id <- order(genes.rank.RF, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.RF[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1


##Uncomment to run Denoising Autoencoder (one hidden layer, gene i clamped to zero) 
##Selection and get gene rankings
# genes.rank.DA <- Denoising.Autoencoder.Selection(X.RNAseq)
# id <- order(genes.rank.DA, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1

##Uncomment to run Denoising Autoencoder (one hidden layer, gene i randomly shuffled) 
##Selection and get gene rankings
# genes.rank.DA1 <- Denoising.Autoencoder.Selection(X.RNAseq, Type = "random")
# id <- order(genes.rank.DA1, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA1[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1


# ##Uncomment to run Denoising Autoencoder (two hidden layer, gene i randomly shuffled) 
# ##Selection and get gene rankings
# genes.rank.DA2 <- Denoising.Autoencoder.Selection(X.RNAseq, Type = "random", 
#                                                   hidden_neurons = c(256, 128))
# id <- order(genes.rank.DA2, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA2[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1

# ##Uncomment to run Denoising Autoencoder (three hidden layer, gene i randomly shuffled) 
# ##Selection and get gene rankings
# genes.rank.DA3 <- Denoising.Autoencoder.Selection(X.RNAseq, Type = "random", 
#                                                   hidden_neurons = c(256, 128, 64))
# id <- order(genes.rank.DA3, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA3[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1

##Uncomment to run Denoising Autoencoder TF (one hidden layer, gene i clamped to zero) 
##Selection and get gene rankings
#dir.data <- (data directory location add here)
# file <- paste(dir.data, "/adj_fly_TRN.csv", sep = "")
# adj.TRN <- as.matrix(read.table(file, sep = ",", header = F))
# genes.rank.DA.TF <- Denoising.Autoencoder.Selection.TF(X.RNAseq, 
#                                                        Type = "zero", 
#                                                        adj.TRN, dir.data)
# id <- order(genes.rank.DA.TF, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA.TF[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1




###Evaluation of Performance using precomputed gene ranks for Methods in the paper
##Distance, and distance z score calculation
nboot <- 785
num.genes <- c(60, 40, 20)
Method <- c("Lasso", "PCA", "Noise", "RF", 
            "RF_TF_nonTF", "Denoising_Autoencoder", 
            "DA_one_Layer", "DA_two_Layer", 
            "DA_three_Layer", "Denoising_Autoencoder_TF", "scPROACTIVE", 
            "RF_geometry", "MLB", "All Genes Top 10")
Method_Name <- c("Lasso", "PCA", "Noise", "RF", "RF TF", "DA", 
                 "DA1", "DA2", "DA3", "DA TF", "scPROACTIVE", "RF Geometry", "MLB", 
                 "All Genes Top 10")
##dir.results <- (Add directory to save results)
##dir.data <- (data directory location add here)
##dir.Methods <- (Add directory to save Method related Data)
for( j in 1:(length(num.genes)) ){
  mcc.dist.zscore.method <- mcc.dist.ratio.zscore.method <- mcc.dist.method <-
    mcc.dist.ratio.method <- matrix(0, length(Method), 
                                    nrow(X.RNAseq))
  mcc.loc.method <- 
    array(0, c(length(Method), nrow(X.RNAseq), 3))
  MCC.list <- list()
  for( i in 1:length(Method) ){
    cat("Num genes = ", num.genes[j], " Method = ", Method[i], "\n")
    if( Method[i] != "scPROACTIVE" & Method[i] != "MLB" & 
        Method[i] != "All Genes Top 10"){
      if( Method[i] != "Random" ){
        file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
        genes <- (read.table(file, sep =",", header = F))
      }else{
        genes <- sample(colnames(X.RNAseq), num.genes[j], replace = F)
      }
    }else{
      if( Method[i] != "All Genes Top 10" ){
        file <- paste(dir.Methods, "/", Method[i], "/", num.genes[j], "genes.csv", 
                      sep = "")
        genes <- as.character(unlist(c(read.table(file, sep = ",", 
                                                  header = F)[1:(num.genes[j]/10), -1])))
      }else{
        genes <- colnames(X.RNAseq)
      }
      
    }
    if( Method[i] != "All Genes Top 10" ){
      id.match <- match(colnames(X.RNAseq), as.character(unlist(genes[1:num.genes[j]])))
    }else{
      id.match <- match(colnames(X.RNAseq), genes)
    }
    
    dm.tmp = new("DistMap",
                 raw.data=as.matrix(dm@raw.data),
                 data=as.matrix(dm@data),
                 insitu.matrix=as.matrix(dm@insitu.matrix[, !is.na(id.match)]),
                 geometry=as.matrix(dm@geometry))
    dm.tmp <- binarizeSingleCellData(dm.tmp, seq(0.15, 0.5, 0.01))
    dm.tmp <- mapCells(dm.tmp)
    dm.tmp@mcc.scores[is.nan(dm.tmp@mcc.scores)] <- 0
    MCC.list[[i]] <- dm.tmp@mcc.scores
    rm(dm.tmp)
  }
 
  file <- paste(dir.data, "/MCC_distance_random_model_Subchallenge_", j, ".RData", 
                sep = "")
  load(file)
  
  for( i in 1:length(Method) ){
    cat("Num genes = ", num.genes[j], " Method = ", Method[i], "\n")
    for( k in 1:nrow(X.RNAseq) ){
      if( Method[i] != "All Genes Top 10" ){
        id.method <- order(MCC.list[[i]][, k], decreasing = T)[1:10]
        mcc.top10 <- sort(MCC.list[[i]][, k], decreasing = T)[1:10]
      }else{
        id.method <- order(dm@mcc.scores[, k], decreasing = T)[1:10]
        mcc.top10 <- sort(dm@mcc.scores[, k], decreasing = T)[1:10]
      }
      id.maxmcc <- order(dm@mcc.scores[, k], decreasing = T)[1:10]
      mcc.top10.orig <- sort(dm@mcc.scores[, k], decreasing = T)[1:10]
      
      if( nnzero(mcc.top10) != 0 ){
        mcc.dist.zscore.method[i, k] <-
          (log10(1 + sqrt(sum((colSums(dm@geometry[id.method, ]*
                                         (cbind(mcc.top10)%*%matrix(1, 1, 3)))/
                                 sum(mcc.top10) -
                                 (dm@geometry[id.maxmcc[1], ]))^2))) -
             mean(mcc.dist.rand[, k]))/sd(mcc.dist.rand[, k])
        mcc.dist.method[i, k] <- log10(1 + sqrt(sum((colSums(dm@geometry[id.method, ]*
                                                               (cbind(mcc.top10)%*%matrix(1, 1, 3)))/
                                                       sum(mcc.top10) -
                                                       (dm@geometry[id.maxmcc[1], ]))^2)))
        
        
      }else{
        mcc.dist.zscore.method[i, k] <-
          (log10(1 + sqrt(sum((colSums(dm@geometry[id.method, ]*
                                         (cbind(rep(1, length(mcc.top10)))%*%matrix(1, 1, 3)))/
                                 length(mcc.top10) -
                                 (dm@geometry[id.maxmcc[1], ]))^2))) -
             mean(mcc.dist.rand[, k]))/sd(mcc.dist.rand[, k])
        mcc.dist.method[i, k] <- log10(1 + sqrt(sum((colSums(dm@geometry[id.method, ]*
                                                               (cbind(rep(1, length(mcc.top10)))%*%matrix(1, 1, 3)))/
                                                       length(mcc.top10) -
                                                       (dm@geometry[id.maxmcc[1], ]))^2)))
        
        
      }
      mcc.loc.method[i, k, ] <- colMeans(dm@geometry[id.method, ])
    }
    
  }
  file <- paste(dir.results, "/MCC_distance_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(mcc.dist.zscore.method, file = file)
  
  dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.results, "/MCC_distance_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(mcc.dist.method, file = file)
  
  dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.results, "/MCC_location_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(mcc.loc.method, file = file)
  
}


##Original location of cells
mcc.loc.orig <- matrix(0, nrow(X.RNAseq), 3)
for( k in 1:nrow(X.RNAseq) ){
  id.maxmcc <- order(dm@mcc.scores[, k], decreasing = T)[1]
  mcc.loc.orig[k, ] <- (dm@geometry[id.maxmcc, ])
}
##dir.results <- "Add result directory"
file <- paste(dir.results, "/MCC_location_original", 
              ".RData", sep = "")
save(mcc.loc.orig, file = file)


###plot distance z-scores
#dir.results <- "ADD result directory path here"
plot.list <- list()
col <- c("magenta", "steelblue2", "lightgreen")
for( j in 1:(length(num.genes)) ){
  file <- paste(dir.results, "/MCC_distance_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  for( i in 1:length(Method) ){
    if( i == 1 ){
      data_bar <- data.frame(x = rep(Method_Name[i], nrow(X.RNAseq)), 
                             y = mcc.dist.zscore.method[i, ])
    }else{
      data_bar1 <- data.frame(x = rep(Method_Name[i], nrow(X.RNAseq)), 
                              y = mcc.dist.zscore.method[i, ])
      data_bar <- rbind(data_bar, data_bar1)
    } 
  }
  file <- paste(dir.results, "/Methods_barplot_new_Subchallenge_", j, ".png", sep = "")
  p <- ggviolin(data_bar, x = "x", y = "y",
                fill = col[j],           # change fill color by mpg_level
                color = "black",            # Set bar border colors to white           # jco journal color palett. see ?ggpar
                x.text.angle = 90,          # Rotate vertically x axis texts
                ylab = "Distance z-score",
                xlab = "",
                legend.title = "",
                ggtheme = theme_minimal(), 
                add = c("mean_se"), 
                font.x = 16, 
                font.y = 16,
                font.tickslab = 14
  ) +  ggtitle(paste("Subchallenge ", 4 - j, sep = ""))
  plot.list[[j]] <- ggplotGrob(p)
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
}
file <- paste(dir.results, "/zscore_distance_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()


###plot distance z-scores barplot
#dir.results <- "ADD result directory path here"
plot.list <- list()
col <- c("magenta", "steelblue2", "lightgreen")
for( j in 1:(length(num.genes)) ){
  dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.results, "/MCC_distance_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  for( i in 1:length(Method) ){
    if( i == 1 ){
      data_bar <- data.frame(x = rep(Method_Name[i], nrow(X.RNAseq)), 
                             y = mcc.dist.zscore.method[i, ])
    }else{
      data_bar1 <- data.frame(x = rep(Method_Name[i], nrow(X.RNAseq)), 
                              y = mcc.dist.zscore.method[i, ])
      data_bar <- rbind(data_bar, data_bar1)
    } 
  }
  id <- order(rowMeans(mcc.dist.zscore.method), decreasing = F)
  data_bar$x <- factor(data_bar$x, levels = Method_Name[id])
  file <- paste(dir.results, "/Methods_mean_distance_barplot_new_Subchallenge_", j, ".png", sep = "")
  p <- ggbarplot(data_bar, x = "x", y = "y",
                 fill = col[j],           # change fill color by mpg_level
                 color = "black",            # Set bar border colors to white           # jco journal color palett. see ?ggpar
                 x.text.angle = 90,          # Rotate vertically x axis texts
                 ylab = "Mean Distance z-score",
                 xlab = "",
                 legend.title = "",
                 ggtheme = theme_minimal(), 
                 add = c("mean_se"), 
                 font.x = 16, 
                 font.y = 16,
                 font.tickslab = 14
  ) +  ggtitle(paste("Subchallenge ", 4 - j, sep = ""))
  plot.list[[j]] <- ggplotGrob(p)
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
}
file <- paste(dir.results, "/zscore_distance_barplot_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()


##Percentage of significantly predicted locations
#dir.results <- "ADD result directory path here"
plot.list <- list()
col <- c("magenta", "steelblue2", "lightgreen")
for( j in 1:(length(num.genes)) ){
  dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.results, "/MCC_distance_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  for( i in 1:length(Method) ){
    id.sign <- which(mcc.dist.zscore.method[i, ] <= -2)
    if( i == 1 ){
      data_bar <- data.frame(x = rep(Method_Name[i], 1), 
                             y = 100*(length(id.sign)/nrow(X.RNAseq)))
    }else{
      data_bar1 <- data.frame(x = rep(Method_Name[i], 1), 
                              y = 100*(length(id.sign)/nrow(X.RNAseq)))
      data_bar <- rbind(data_bar, data_bar1)
    } 
  }

  file <- paste(dir.results, "/Methods_barplot_significant_percentage_new_Subchallenge_", j, ".png", sep = "")
  # p <- ggbarplot(data_bar[order(data_bar$y, decreasing = T), ], x = "x", y = "y",
  #               fill = col[j],           # change fill color by mpg_level
  #               color = "black",            # Set bar border colors to white           # jco journal color palett. see ?ggpar
  #               x.text.angle = 90,          # Rotate vertically x axis texts
  #               ylab = "Fraction of significant predictions",
  #               xlab = "",
  #               legend.title = "",
  #               ggtheme = theme_minimal(),
  #               font.x = 14, 
  #               font.y = 14,
  #               font.tickslab = 12, 
  #               sort.val = "desc"
  # )
  data_bar <- data_bar[order(data_bar$y, decreasing = T), ]
  data_bar$x <- factor(data_bar$x, levels = unique(data_bar$x))
  p <- ggplot(data_bar, 
              mapping = aes(x = x, y = y)) +
    geom_bar(stat = "identity", fill = col[j], color = "black") + 
    theme(axis.text.x = element_text(size=13, angle = 90),
          axis.title=element_text(size=14,face="bold")) + 
    ylab("% significant predictions") + xlab("") + 
    ggtitle(paste("Subchallenge ", 4 - j, sep = ""))
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
  plot.list[[j]] <- ggplotGrob(p)
}
file <- paste(dir.results, "/Percentage_significant_locations_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()

##Overall ranking significant locations
#dir.results <- "ADD result directory path here"
col <- c("magenta", "steelblue2", "lightgreen")
num.genes <- c(60, 40, 20)
Method.ranking.better <- Method.ranking.lesser <- Method.ranking <- list()
plot.list <- list()
for( j in 1:(length(num.genes)) ){
  file <- paste(dir.results, "/MCC_distance_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  
  dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.results, "/MCC_distance_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  Method.ranking.better[[j]] <- Method.ranking.lesser[[j]] <- 
    matrix(0, length(Method), length(Method))
  Method.ranking[[j]] <- rep(0, length(Method) - 1)
  for( i in 1:(length(Method)) ){
    if( i == 1 ){
      id.better <- which(mcc.dist.zscore.method[i, ] <= -2)
      id.lesser <- which(mcc.dist.zscore.method[i, ] >= 2)
    }else{
      id.better <- union(id.better, which(mcc.dist.zscore.method[i, ] <= -2))
      id.lesser <- union(id.lesser, which(mcc.dist.zscore.method[i, ] >= 2))
    }
  }
  for( i in 1:(length(Method) - 1) ){
    for( k in (i + 1):length(Method) ){
      cat("Subchallenge = ", j, " Method 1 = ", i, " Method 2 = ", k, "\n")
      id <- id.better
      z1 <- mcc.dist.zscore.method[i, id]
      z2 <- mcc.dist.zscore.method[k, id]
      P1 <- which(z1 - z2 < 0)
      P2 <- which(z2 - z1 < 0)
      Method.ranking.better[[j]][i, k] <- length(P1)/(length(P1) + length(P2))
      Method.ranking.better[[j]][k, i] <- length(P2)/(length(P1) + length(P2))
      
      id <- id.lesser
      z1 <- mcc.dist.zscore.method[i, id]
      z2 <- mcc.dist.zscore.method[k, id]
      P1 <- which(z1 - z2 < 0)
      P2 <- which(z2 - z1 < 0)
      Method.ranking.lesser[[j]][i, k] <- length(P1)/(length(P1) + length(P2))
      Method.ranking.lesser[[j]][k, i] <- length(P2)/(length(P1) + length(P2))
      
    }
  }
  rank.better <- rank.lesser <- rep(0, length(Method) - 1)
  rank.better[order(Method.ranking.better[[j]][, length(Method)], 
                    decreasing = T)[-length(Method)]] <- 1:(length(Method) - 1)
  rank.lesser[order(Method.ranking.lesser[[j]][, length(Method)], 
                    decreasing = T)[-length(Method)]] <- 1:(length(Method) - 1)
  Method.ranking[[j]] <- (rank.lesser + rank.better)/2
  
  dir.results <- "/home/atheistpoet/Desktop/Work/Study/UIUC/Semester_1/CS598JP/Project/DREAM_single_cell/Data"
  file <- paste(dir.results, "/Ranking_Methods_all_new_Subchallenge_", 
                j, 
                ".csv", sep = "")
  write.table(cbind(sort(rank(Method.ranking[[j]]), decreasing = F), 
                    Method[order(Method.ranking[[j]], decreasing = F)]), 
              file, sep = ",", row.names = F, col.names = F, 
              quote = F)
  
  data_bar <- data.frame(x = rep(Method_Name[-length(Method)], 2), 
                         y = c(Method.ranking.better[[j]][-length(Method), 
                                                          length(Method)], 
                               Method.ranking.lesser[[j]][-length(Method), 
                                                          length(Method)]), 
                         Significance = c(rep("zscore <= -2", length(Method) - 1), 
                                          rep("zscore >= 2", length(Method) - 1)))
  id1 <- order(Method.ranking.better[[j]][-length(Method), length(Method)], decreasing = T)
  data_bar$x <- factor(data_bar$x, levels = Method_Name[id1])
  p <- ggplot(data_bar, 
              mapping = aes(x = x, y = y)) +
    geom_bar(data = data_bar, 
             aes(x = x, y = y, group = Significance, fill = Significance), stat = "identity", 
             color = "black", position = position_dodge(0.9)) + 
    theme(axis.text.x = element_text(size=13, angle = 90),
          axis.title=element_text(size=14,face="bold")) + 
    ylab("Distance score") + xlab("") + 
    ggtitle(paste("Subchallenge ", 4 - j, sep = ""))
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
  plot.list[[j]] <- ggplotGrob(p)
}
file <- paste(dir.results, "/ranking_barplot_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()