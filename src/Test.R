##Test: Sample script to generate gene rankings for each method.....
## Also, shows how to use pre-computed ranks, given in the data directory for each
## method, to generate plots from the paper.

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
#Uncomment below and give path to data directory
#dir.data <- <add data directory path here>
dm <- load.data(dir.data, RNAseq.full = F)

###RNAseq data
X.RNAseq <- t(dm@data)

i <- 1
#Uncomment below and give path to Methods directory which is a subdirectory in
#Data Directory
#dir.Methods <- <directory for storing gene rankings from Methods>

##Uncomment below to run lasso and get gene rankings
###BEGIN:LASSO Ranking
#genes.rank.lasso <- LASSO.Selection(X.RNAseq, lambda, nfolds, clsters)
# id <- order(genes.rank.lasso, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.lasso[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END:LASSO Ranking

##Uncomment below to run PCA selection and get gene rankings
###BEGIN: PCA Ranking
#num.PC <- 3
# genes.rank.PCA <- PCA.Selection(X.RNAseq, num.PC)
# id <- order(genes.rank.PCA, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.PCA[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: PCA Ranking

##Uncomment below to run Noise selection and get gene rankings
###BEGIN: Noise Ranking
# genes.rank.noise <- Noise.Selection(X.RNAseq)
# id <- order(genes.rank.noise, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.noise[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: PCA Ranking


##Uncomment below to run Random Forest Selection and get gene rankings
###BEGIN: Random Forest Ranking
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
###END: PCA Ranking


##Uncomment Below to run Denoising Autoencoder (one hidden layer, gene i clamped to zero) 
##Selection and get gene rankings
###BEGIN: Denoising Autoencoder Ranking
# genes.rank.DA <- Denoising.Autoencoder.Selection(X.RNAseq)
# id <- order(genes.rank.DA, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: Denoising Autoencoder Ranking

##Uncomment to run Denoising Autoencoder (one hidden layer, gene i randomly shuffled) 
##Selection and get gene rankings
###BEGIN: Denoising Autoencoder Ranking
# genes.rank.DA1 <- Denoising.Autoencoder.Selection(X.RNAseq, Type = "random")
# id <- order(genes.rank.DA1, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA1[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: Denoising Autoencoder Ranking


# ##Uncomment to run Denoising Autoencoder (two hidden layer, gene i randomly shuffled) 
# ##Selection and get gene rankings
###BEGIN: Denoising Autoencoder Ranking
# genes.rank.DA2 <- Denoising.Autoencoder.Selection(X.RNAseq, Type = "random", 
#                                                   hidden_neurons = c(256, 128))
# id <- order(genes.rank.DA2, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA2[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: Denoising Autoencoder Ranking

# ##Uncomment to run Denoising Autoencoder (three hidden layer, gene i randomly shuffled) 
# ##Selection and get gene rankings
###BEGIN: Denoising Autoencoder Ranking
# genes.rank.DA3 <- Denoising.Autoencoder.Selection(X.RNAseq, Type = "random", 
#                                                   hidden_neurons = c(256, 128, 64))
# id <- order(genes.rank.DA3, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA3[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: Denoising Autoencoder Ranking

##Uncomment to run Denoising Autoencoder TF (one hidden layer, gene i clamped to zero) 
##Selection and get gene rankings
###BEGIN: Denoising Autoencoder Ranking
#dir.data <- <data directory path add here>
# file <- paste(dir.data, "/adj_fly_TRN.RData", sep = "")
#load(file)
# genes.rank.DA.TF <- Denoising.Autoencoder.Selection.TF(X.RNAseq, 
#                                                        Type = "zero", 
#                                                        adj.TRN, dir.data)
# id <- order(genes.rank.DA.TF, decreasing = F)
# file <- paste(dir.Methods, "/", Method[i], "_", "Genes.csv", sep = "")
# write.table(c(genes.rank.DA.TF[id]), file, sep = ",", row.names = F,
#             col.names = F, quote = F)
# i <- i + 1
###END: Denoising Autoencoder Ranking


##Uncomment below to generate random null model....
###BEGIN:Random Null Model
# num.genes <- c(60, 40, 20) #number of genes in subset in each subchallenge
# nboot <- 1000 #number of gene subset samples in random null model
# dir.data <- <add path to data directory>
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
# tmp <- Random.null.model(X.RNAseq, dm, nboot, num.genes, dir.random)
# tmp <- Random.MCC.dist(X.RNAseq, dir.random, dir.data, num.genes, nboot)
###END:Random Null Model



###Evaluation of Performance using precomputed gene ranks for Methods in the paper
##Distance, and distance z score calculation
nboot <- 785  #number of gene ranking samples in random null models
num.genes <- c(60, 40, 20)
Method <- c("Lasso", "PCA", "Noise", "RF", 
            "RF_TF_nonTF", "Denoising_Autoencoder", 
            "DA_one_Layer", "DA_two_Layer", 
            "DA_three_Layer", "Denoising_Autoencoder_TF", "scPROACTIVE", 
            "RF_geometry", "MLB", "All Genes Top 10")
Method_Name <- c("Lasso", "PCA", "Noise", "RF", "RF TF", "DA", 
                 "DA1", "DA2", "DA3", "DA TF", "scPROACTIVE", "RF Geometry", "MLB", 
                 "All Genes Top 10")
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
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
 
  file <- paste(dir.random, "/MCC_distance_random_model_Subchallenge_", j, ".RData", 
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
#Uncomment belwo and add result directory path
##dir.results <- "Add result directory"
file <- paste(dir.results, "/MCC_location_original", 
              ".RData", sep = "")
save(mcc.loc.orig, file = file)


###plot distance z-scores'
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
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
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
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
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
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
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
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




###vISH results and plots
##Significance for vISH for all genes
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
file <- paste(dir.raandom, "/mismatch_VFISH_random.RData", sep = "")
load(file)

nboot <- 785
num.genes <- c(60, 40, 20)
# Method <- c("Lasso", "PCA", "Noise", "RF", "RF_all", 
#             "RF_TF_nonTF", 
#             "Denoising_Autoencoder", "Denoising_Autoencoder_TF", "scPROACTIVE", 
#             "RF_geometry", "MLB", "All Genes Top 10")
# Method_Name <- c("Lasso", "PCA", "Noise", "RF sqrt", "RF all", "RF TF prior", 
#                  "DA", "DA prior", "scPROACTIVE", "RF Geometry", "MLB", 
#                  "All Genes Top 10")

mcc.mismatch.method <- 
  array(0, c((length(num.genes)), length(Method), ncol(X.RNAseq), 
             nrow(dm@geometry)))
mcc.mismatch.zscore.method <- array(0, c((length(num.genes)), length(Method), 
                                         nrow(dm@geometry)))
for( j in 1:(length(num.genes)) ){
  pha.orig <- matrix(0, ncol(X.RNAseq), nrow(dm@geometry))
  pha.method <- array(0, c(length(Method), ncol(X.RNAseq), nrow(dm@geometry)))
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
    
    for( k in 1:ncol(X.RNAseq) ){
      cat("Gene = ", k, "\n")
      pha.orig[k, ] <- computeVISH(dm, colnames(X.RNAseq)[k], threshold = 0.75)
      pha.method[i, k, ] <- computeVISH(dm.tmp, colnames(X.RNAseq)[k], threshold = 0.75)
      VISH.diff <- abs(pha.orig[k, ] - pha.method[i, k, ])/max(pha.orig[k, ])
      id.nnz <- which(VISH.diff >= 1e-1)
      mcc.mismatch.method[j, i, k, id.nnz] <- 1
    }
    mcc.mismatch.zscore.method[j, i, ] <- 
      (colMeans(mcc.mismatch.method[j, i, , ]) - colMeans(mismatch.rnd[j, , ]))/
      apply(mismatch.rnd[j, , ], 2, sd)
    rm(dm.tmp)
  }
  
  file <- paste(dir.results, "/vISH_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(pha.method, file = file)
  
  file <- paste(dir.results, "/vISH_all_genes_top10_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(pha.orig, file = file)
  
  file <- paste(dir.results, "/vISH_mismatch_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(mcc.mismatch.method, file = file)
  
  file <- paste(dir.results, "/vISH_mismatch_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  save(mcc.mismatch.zscore.method, file = file)
  
  
}



##mismatch z-score barplot
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
plot.list <- list()
col <- c("magenta", "steelblue2", "lightgreen")
for( j in 1:(length(num.genes)) ){
  file <- paste(dir.results, "/vISH_mismatch_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  for( i in 1:length(Method) ){
    if( i == 1 ){
      data_bar <- data.frame(x = rep(Method_Name[i], nrow(dm@geometry)), 
                             y = mcc.mismatch.zscore.method[j, i, ])
    }else{
      data_bar1 <- data.frame(x = rep(Method_Name[i], nrow(dm@geometry)), 
                              y = mcc.mismatch.zscore.method[j, i, ])
      data_bar <- rbind(data_bar, data_bar1)
    } 
  }
  id <- order(rowMeans(mcc.mismatch.zscore.method[j, , ]), decreasing = F)
  data_bar$x <- factor(data_bar$x, levels = Method_Name[id])
  file <- paste(dir.results, "/Methods_barplot_mismatch_zcore_new_Subchallenge_", j, ".png", sep = "")
  p <- ggbarplot(data_bar, x = "x", y = "y",
                 fill = col[j],           # change fill color by mpg_level
                 color = "black",            # Set bar border colors to white           # jco journal color palett. see ?ggpar
                 x.text.angle = 90,          # Rotate vertically x axis texts
                 ylab = "vISH z-score",
                 xlab = "",
                 legend.title = "",
                 ggtheme = theme_minimal(), 
                 add = c("mean_se"), 
                 font.x = 14, 
                 font.y = 14,
                 font.tickslab = 12
  ) + ggtitle(paste("Subchallenge ", 4 - j, sep = ""))
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
  plot.list[[j]] <- ggplotGrob(p)
}
file <- paste(dir.results, "/zscore_vISH_mismatch_barplot_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()



##Percentage significant locations vISH
plot.list <- list()
break.y <- seq(-1, 1.5, by = 0.5)
label.y <- rep("NA", length(break.y))
for(  i in 1:length(break.y) ){
  label.y[i] <-paste("10^", break.y[i], sep = "")
}
label.y <- parse(text = label.y)
col <- c("magenta", "steelblue2", "lightgreen")
for( j in 1:(length(num.genes)  - 1) ){
  file <- paste(dir.results, "/vISH_mismatch_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  for( i in 1:length(Method) ){
    id.sign <- which(mcc.mismatch.zscore.method[j, i, ] <= -2)
    if( i == 1 ){
      data_bar <- data.frame(x = rep(Method_Name[i], 1), 
                             y = length(id.sign)/nrow(dm@geometry))
    }else{
      data_bar1 <- data.frame(x = rep(Method_Name[i], 1), 
                              y = length(id.sign)/nrow(dm@geometry))
      data_bar <- rbind(data_bar, data_bar1)
    } 
  }
  data_bar$y <- log10(100*data_bar$y)
  file <- paste(dir.results, "/Methods_barplot_significant_percentage_vISH_new_Subchallenge_", j, ".png", sep = "")
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
    ggtitle(paste("Subchallenge ", 4 - j, sep = "")) + 
    scale_y_continuous(breaks = break.y, 
                       labels = label.y)
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
  plot.list[[j]] <- ggplotGrob(p)
}
file <- paste(dir.results, "/Percentage_significant_locations_vISH_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()



##Ranking significant bins vISH
#Uncomment below and add result directory path
##dir.results <- <Add path of directory to save results>
#Uncomment below and add data directory path
##dir.data <- <data directory pats add here)
#Uncomment below and add Methods directory path
##dir.Methods <- <Add path to directory to save Method related evaluation Data: This
##                is a subdirectory under dir.data>
#Uncomment below and add Random directory path
# dir.random <- <add path to 'random' directory, which is a subdirectory under 
#                 data directory>
col <- c("magenta", "steelblue2", "lightgreen")
Method.ranking.better <- Method.ranking.lesser <- Method.ranking <- list()
plot.list <- plot.list1 <- list()
for( j in 1:(length(num.genes)) ){
  file <- paste(dir.results, "/vISH_mismatch_zscore_Methods_all_Subchallenge_", 
                j, 
                ".RData", sep = "")
  load(file)
  
  Method.ranking.better[[j]] <- Method.ranking.lesser[[j]] <- 
    matrix(0, length(Method), length(Method))
  Method.ranking[[j]] <- rep(0, length(Method) - 1)
  for( i in 1:(length(Method)) ){
    if( i == 1 ){
      id.better <- which(mcc.mismatch.zscore.method[j, i, ] <= -2)
      id.lesser <- which(mcc.mismatch.zscore.method[j, i, ] >= 2)
    }else{
      id.better <- union(id.better, which(mcc.mismatch.zscore.method[j, i, ] <= -2))
      id.lesser <- union(id.lesser, which(mcc.mismatch.zscore.method[j, i, ] >= 2))
    }
  }
  for( i in 1:(length(Method) - 1) ){
    for( k in (i + 1):length(Method) ){
      cat("Subchallenge = ", j, " Method 1 = ", i, " Method 2 = ", k, "\n")
      id <- id.better
      z1 <- mcc.mismatch.zscore.method[j, i, id]
      z2 <- mcc.mismatch.zscore.method[j, k, id]
      P1 <- which(z1 - z2 < 0)
      P2 <- which(z2 - z1 < 0)
      Method.ranking.better[[j]][i, k] <- length(P1)/(length(P1) + length(P2))
      Method.ranking.better[[j]][k, i] <- length(P2)/(length(P1) + length(P2))
      
      id <- id.lesser
      z1 <- mcc.mismatch.zscore.method[j, i, id]
      z2 <- mcc.mismatch.zscore.method[j, k, id]
      P1 <- which(z1 - z2 < 0)
      P2 <- which(z2 - z1 < 0)
      Method.ranking.lesser[[j]][i, k] <- length(P1)/(length(P1) + length(P2))
      Method.ranking.lesser[[j]][k, i] <- length(P2)/(length(P1) + length(P2))
      
    }
  }
  rank.better <- rank.lesser <- rep(0, length(Method) - 1)
  rank.better[order(Method.ranking.better[[j]][, 12], 
                    decreasing = T)[-length(Method)]] <- 1:(length(Method) - 1)
  rank.lesser[order(Method.ranking.lesser[[j]][, 12], 
                    decreasing = T)[-length(Method)]] <- 1:(length(Method) - 1)
  Method.ranking[[j]] <- (rank.lesser + rank.better)/2
  
  file <- paste(dir.results, "/Ranking_Methods_vISH_all_new_Subchallenge_", 
                j, 
                ".csv", sep = "")
  write.table(cbind(sort(rank(Method.ranking[[j]]), decreasing = F), 
                    Method[order(Method.ranking[[j]], decreasing = F)]), 
              file, sep = ",", row.names = F, col.names = F, 
              quote = F)
  
  data_bar <- data.frame(x = rep(Method_Name[-length(Method)], 2), 
                         y = c(Method.ranking.better[[j]][-12, 12], 
                               Method.ranking.lesser[[j]][-12, 12]), 
                         Significance = c(rep("zscore <= -2", 11), rep("zscore >= 2", 11)))
  id1 <- order(Method.ranking.better[[j]][-12, 12], decreasing = T)
  data_bar$x <- factor(data_bar$x, levels = Method_Name[id1])
  p <- ggplot(data_bar, 
              mapping = aes(x = x, y = y)) +
    geom_bar(data = data_bar, 
             aes(x = x, y = y, group = Significance, fill = Significance), stat = "identity", 
             color = "black", position = position_dodge(0.9)) + 
    theme(axis.text.x = element_text(size=13, angle = 90),
          axis.title=element_text(size=14,face="bold")) + 
    ylab("vISH score") + xlab("") + 
    ggtitle(paste("Subchallenge ", 4 - j, sep = ""))
  png(file, width = 1200, height = 720, res = 120)
  print(p)
  dev.off()
  
  
  library(reshape2)
  mat.better <- Method.ranking.better[[j]][-12, -12]
  for( i in 1:10 ){
    for( k in 2:11 ){
      mat.better[k, i] <- mat.better[i, k]
    }
  }
  colnames(mat.better) <- rownames(mat.better) <- Method_Name[-12]
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist(cormat)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  # Reorder the correlation matrix
  cormat <- mat.better
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                         name="vISH score (z <= -2)") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+ xlab("") + 
    ylab("")
  coord_fixed()
  # Print the heatmap
  plot.list[[j]] <- ggplotGrob(ggheatmap)
  
  
  
  library(reshape2)
  mat.better <- Method.ranking.lesser[[j]][-12, -12]
  for( i in 1:10 ){
    for( k in 2:11 ){
      mat.better[k, i] <- mat.better[i, k]
    }
  }
  colnames(mat.better) <- rownames(mat.better) <- Method_Name[-12]
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist(cormat)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  # Reorder the correlation matrix
  cormat <- mat.better
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                         name="vISH score (z >= 2)") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+ xlab("") + 
    ylab("")
  coord_fixed()
  plot.list1[[j]] <- ggplotGrob(ggheatmap)
  # Print the heatmap
}
file <- paste(dir.results, "/ranking_vISH_better_barplot_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()

file <- paste(dir.results, "/ranking_vISH_lesser_barplot_new.png", sep = "")
png(file, width = 1500, height = 1000, res = 130)
grid.arrange(
  grobs = plot.list1,
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(3:4))
)
dev.off()

