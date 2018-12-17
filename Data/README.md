#Details of the data files
The details of all the data files is as follows:

1) DAVID_gene_DREAM_84_name_to_flybase.txt -- Contains list for conversion between official gene names and flybase id for fruit fly D.melanogatser. This is used for matching genes in RNAseq and insitu data with the genes in the transcriptional regulatory network (TRN) while using Denoisng Autoencoder with TRN.

2) MCC_distance_model_Subchallenge_<1, 2, 3>.RData -- R data object that contains distance score for cellular location predictions generated from random null model for each sub-challenge. Distance is calculated from the original cellular positions. For each subchallenge (1, 2, or 3), this object conatins a matrix 'mcc.dist.rand', the size of which is nboot (number of samples of subset of genes from random null mode) X num_cells (number of cells in the embryo).

3) adj_fly_TRN.RData -- R data object that conatins the adjacency matrix for transcriptional regulatory network for fruit fly D. melanogaster.

4) bdtnp.txt -- The reference insitu hybridization database from the Berkeley Drosophila Transcription Network Project. This is matrix of in-situ hybridization expression data witjh 84 genes (columns) and 3039 rows (locations on embryo). The order in the rows follows that in the file 'geometry.txt' (which is the x, y, z coordinates of the locations on the embryo.).

5) binarized_bdtnp.csv -- Binarized insitu hybridization expression data.

6) dge_normalized_driver_genes.txt -- Normalized single-cell RNAseq data for driver genes. There are 84 driver genes (columns), and 1297 rows (cells).

7) dge_raw.txt -- raw single-cell RNAseq data for all genes (8924).

8) geometry.txt -- geometry file containing the x, y, z coordinates for locations on the embryo.

The directory also contains a subdirectory 'Methods' that has pre-computed gene list for all methods in (decreasing order of importance for cellular location prediction).
