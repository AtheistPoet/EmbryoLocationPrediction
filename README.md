# EmbryoLocationPrediction
This repository contains code for predicting cellular locations on the embryo of fruit fly D. melanogatser using single cell RNAseq and insitu hybridization data. The main task is to identify redundancy present in the driver genes (susbet of genes for which insitu hybridization data is available). Further, we identify a subset of the driver genes only using RNAseq data and prior biological knowledge. We compare the performance of different methods in this feature selection task. We also introduce distance and virtual In-Situ Hybridization (vISH) metrics to compare the performance.

The directory is divided in three sub-directories--'Data', 'Results', 'src'.

1) 'Data' has all the data required for reproducing our results.
2) 'Results' is for saving all the generated results.
3)  'src' has all the source code.
