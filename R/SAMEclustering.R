# !/usr/bin/env Rscript
# 11/10/2017

#' @importFrom SC3 sc3_estimate_k sc3 sc3_run_svm
#' @importFrom e1071 svm
#' @importFrom SingleCellExperiment SingleCellExperiment normcounts normcounts<- logcounts logcounts<-
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<- assayNames
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
sc3_SAME <- function(inputTags, gene_filter, svm_num_cells, SEED){
  exp_cell_exprs <- NULL
  sc3OUTPUT <- NULL
  
  # cell expression
  ### For count data, it would be normalized by the total cound number and then log2 transformed
  exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
  normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
  logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
  
  rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
  exp_cell_exprs <- exp_cell_exprs[!duplicated(rowData(exp_cell_exprs)$feature_symbol), ]
  
  ### Estimating optimal number of clustering
  exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
  optimal_K <- metadata(exp_cell_exprs)$sc3$k_estimation
  
  ### Clustering by SC3 at the optimal K
  if (ncol(inputTags) < svm_num_cells){
    #print(optimal_K)
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter, n_cores = 1, rand_seed = SEED)
  } else if (ncol(inputTags) >= svm_num_cells){
    ### Runing SVM
    exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = gene_filter,
                          svm_max = svm_num_cells, svm_num_cells = svm_num_cells, n_cores = 1, rand_seed = SEED)
    exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
  }
  
  ### Exporting SC3 results
  p_Data <- colData(exp_cell_exprs)
  col_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
  sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
  return(sc3OUTPUT)
}

#' @importFrom cidr scDataConstructor determineDropoutCandidates wThreshold scDissim scPCA nPC scCluster
cidr_SAME <- function(inputTags, percent_dropout, nPC.cidr, SEED){
  set.seed(SEED)
  
  cidrOUTPUT <- NULL
  if(is.null(percent_dropout)){
    inputTags_cidr <- inputTags
  } else{
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags)*100
    inputTags_cidr <- inputTags[-c(which(dropouts <= percent_dropout), which(dropouts >= 100 - percent_dropout)),]
  }
  
  cidrOUTPUT <- scDataConstructor(inputTags_cidr, tagType = "raw")
  cidrOUTPUT <- determineDropoutCandidates(cidrOUTPUT)
  cidrOUTPUT <- wThreshold(cidrOUTPUT)
  cidrOUTPUT <- scDissim(cidrOUTPUT)
  cidrOUTPUT <- scPCA(cidrOUTPUT)
  cidrOUTPUT <- nPC(cidrOUTPUT)
  
  ### Define nPC
  if(!is.null(nPC.cidr)) {
    cidrOUTPUT@nPC <- nPC.cidr
  } else {
    nPC.cidr <- cidrOUTPUT@nPC
  }
  
  ### Clustering by CIDR.
  # The optimal clustering number is determined automatically
  cidrOUTPUT <- scCluster(cidrOUTPUT, nPC = nPC.cidr)
  
  return(cidrOUTPUT)
}

#' @import Seurat
#' @importFrom methods .hasSlot
seurat_SAME <- function(inputTags, nGene_filter = TRUE, low.genes, high.genes, nPC.seurat, resolution, SEED){
  seuratOUTPUT <- NULL
  
  # Initialize the Seurat object with the raw data (non-normalized data)
  # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
  seuratOUTPUT <- CreateSeuratObject(counts = inputTags, min.cells = 0, min.genes = low.genes, project = "single-cell clustering")
  
  # Filter out the cells expressing too few or too many genes
  if (nGene_filter == TRUE){
     seuratOUTPUT <- subset(object = seuratOUTPUT, subset = nFeature_RNA > low.genes & nFeature_RNA < high.genes)
  }
  
  # Perform log-normalization, first scaling each cell to a total of 1e4 molecules (as in Macosko et al. Cell 2015)
  seuratOUTPUT = NormalizeData(object = seuratOUTPUT, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Detection of variable genes across the single cells
  seuratOUTPUT = FindVariableFeatures(object = seuratOUTPUT, selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  all.genes <- rownames(seuratOUTPUT)
  seuratOUTPUT <- ScaleData(object = seuratOUTPUT, features = all.genes)
  
  ### Perform linear dimensional reduction
  if (nPC.seurat <= 20){
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), npcs = 20, seed.use = SEED, verbose = FALSE)
    seuratOUTPUT <- FindNeighbors(seuratOUTPUT, dims = 1:20, verbose = FALSE)
  } else {
    seuratOUTPUT <- RunPCA(object = seuratOUTPUT, features = VariableFeatures(object = seuratOUTPUT), npcs = nPC.seurat, seed.use = SEED, verbose = FALSE)
    seuratOUTPUT <- FindNeighbors(seuratOUTPUT, dims = 1:20, verbose = FALSE)
  }
  
  seuratOUTPUT <- FindClusters(object = seuratOUTPUT, resolution = resolution, verbose = FALSE)
  
  ### Complementing the missing data
  if (length(seuratOUTPUT@active.ident) < ncol(inputTags)){
    seurat_output <- matrix(NA, ncol = ncol(inputTags), byrow = T)
    colnames(seurat_output) <- colnames(inputTags)
    seurat_retained <- t(as.matrix(as.numeric(seuratOUTPUT@active.ident)))
    colnames(seurat_retained) <- names(seuratOUTPUT@active.ident)
    for (i in 1:ncol(seurat_retained)){
      seurat_output[1,colnames(seurat_retained)[i]] <- seurat_retained[1,colnames(seurat_retained)[i]]
    }
  } else {
    seurat_output <- t(as.matrix(as.numeric(seuratOUTPUT@active.ident)))
  }
  
  return(seurat_output)
}

#' @importFrom Rtsne Rtsne
#' @importFrom ADPclust adpclust
#' @importFrom S4Vectors var
#' @importFrom stats kmeans
tSNE_kmeans_SAME <- function(inputTags, percent_dropout, dimensions, perplexity, k.min, k.max, var_genes, SEED){
  input_lcpm <- NULL
  tsne_input <- NULL
  tsne_output <- NULL
  tsne_kmeansOUTPUT <- NULL
  adpOUTPUT <- NULL
  
  if (is.null(percent_dropout)){
    inputTags_tsne <- inputTags
  } else{
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags)*100
    inputTags_tsne <- inputTags[-c(which(dropouts <= percent_dropout), which(dropouts >= 100 - percent_dropout)),]
  }
  
  ### Data tranformation
  ### If the input data is original count data or CPM, it would be tranformed to CPM
  tsne_input <- log2(t(t(inputTags_tsne)/colSums(inputTags_tsne))*1000000+1)
  
  if (is.null(var_genes)){
    set.seed(SEED)
    tsne_output <- Rtsne(t(tsne_input), dims = dimensions, perplexity = perplexity, check_duplicates = FALSE)
  } else{
    se_genes = rep(NA, nrow(tsne_input))
    for (i in 1:nrow(tsne_input)){
      se_genes[i] = sqrt(var(tsne_input[i,])/length(tsne_input[i,]))
    }
    decreasing_rank = order(se_genes, decreasing = TRUE)
    
    set.seed(SEED)
    
    tsne_output <- Rtsne(t(tsne_input[decreasing_rank[1:var_genes],]), dims = dimensions, perplexity = perplexity)
  }
  
  ### Determining the optimal cluster number (k) and centroid by ADPclust
  adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise",centroids="auto", nclust = k.min:k.max)
  
  ### Clustering the cells by kmeans
  tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers,], adpOUTPUT$nclust)
  
  return(tsne_kmeansOUTPUT)
}

#' @importFrom SIMLR SIMLR_Estimate_Number_of_Clusters SIMLR SIMLR_Large_Scale
SIMLR_SAME <- function(inputTags, percent_dropout, k.min, k.max, SEED){
  set.seed(SEED)
  
  #X input is genebycell matrix
  #c is number of clusters
  if(is.null(percent_dropout)){
    inputTags_simlr <- inputTags
  } else{
    dropouts <- rowSums(inputTags == 0)/ncol(inputTags)*100
    inputTags_simlr <- inputTags[-c(which(dropouts <= percent_dropout), which(dropouts >= 100 - percent_dropout)),]
  }
  
  k_range <- k.min:k.max
  
  best_k=SIMLR_Estimate_Number_of_Clusters(log10(t(t(inputTags_simlr)/colSums(inputTags_simlr))*1000000+1), NUMC = k_range, cores.ratio = 1)
  k <- which.min(best_k$K1) + k_range[1] - 1
  
  if (dim(inputTags_simlr)[2] < 1000) {
    simlrOUTPUT <- SIMLR(log10(inputTags_simlr+1), c = k, cores.ratio = 0)
  } else {
    simlrOUTPUT <- SIMLR_Large_Scale(log10(inputTags_simlr+1), c = k)
  }
  
  return(simlrOUTPUT)
}

#' @title SAME
#'
#' @description This function performs single-cell clustering using five state-of-the-art methods,
#' SC3, CIDR, Seurat, tSNE+kmeans and SIMLR.
#'
#' @param inputTags a G*N matrix with G genes and N cells.
#' @param mt_filter is a boolean variable that defines whether to filter outlier cells according to mitochondrial gene percentage.
#' Default is "TRUE".
#' @param mt.pattern defines the pattern of mitochondrial gene names in the data, for example, \code{mt.pattern = "^MT-"} for human and \code{mt.pattern = "^mt-"} for mouse. Only applied when \code{mt_filter = TRUE}
#' Default is \code{mt.pattern = "^MT-"}.
#' @param mt.cutoff defines a high cutoff of mitochondrial percentage (Default is 0.1) that cells having higher percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.
#' @param percent_dropout defines a low cutoff of gene percentage that genes expressed in less than \code{percent_dropout}% or more than (100 - \code{percent_dropout})% of cells are removed for CIDR, tSNE + k-means and SIMLR clustering. 
# Default is 10.
#' @param SC3 is a boolean variable that defines whether to cluster cells using SC3 method.
#' Default is "TRUE".
#' @param gene_filter is a boolean variable that defines whether to perform gene filtering before SC3 clustering, when \code{SC3 = TRUE}.
#' @param svm_num_cells, if \code{SC3 = TRUE}, then defines the mimimum number of cells above which SVM will be run.
#' @param CIDR is a boolean parameter that defines whether to cluster cells using CIDR method.
#' Default is "TRUE".
#' @param nPC.cidr defines the number of principal coordinates used in CIDR clustering, when \code{CIDR = TRUE}.
#' Default value is esimated by \code{nPC} of \code{CIDR}.
#' @param Seurat is a boolean variable that defines whether to cluster cells using Seurat method.
#' Default is "TRUE".
#' @param nGene_filter is a boolean variable that defines whether to filter outlier cells according to unique gene count before Seurat clustering.
#' Default is "TRUE".
#' @param low.genes defines is a low cutoff of unique gene counts (Default is 200) that cells expressing less than 200 genes are filtered out, when \code{nGene_filter = TRUE}. 
#' @param high.genes defines is a high cutoff of unique gene counts (Default is 8000) that cells expressing more than 8000 genes are filtered out, when \code{nGene_filter = TRUE}. 
#' @param nPC.seurat defines the number of principal components used in Seurat clustering, when \code{Seurat = TRUE}.
#' Default is \code{nPC.seurat = nPC.cidr}.
#' @param resolution defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}. Default is \code{resolution = 0.7}.
#' @param tSNE is a boolean variable that defines whether to cluster cells using t-SNE + k-means method.
#' Default is "TRUE".
#' @param dimensions sets the number of dimensions wanted to be retained in t-SNE step. Default is 3.
#' @param perplexity sets the perplexity parameter for t-SNE dimension reduction. Default is 30 when number of cells \code{>=200}.
#' @param tsne_min_cells defines the number of cells in input dataset below which
#' \code{tsne_min_perplexity=10} would be employed for t-SNE step. Default is 200.
#' @param tsne_min_perplexity sets the perplexity parameter of t-SNE step for small datasets (number of cells \code{<200}).
#' @param var_genes defines the number of variable genes used by t-SNE analysis, when \code{tSNE = TRUE}.
#' @param SIMLR is a boolean variable that defines whether to cluster cells using t-SNE + k-means method.
#' Default is "TRUE".
#' @param diverse is a boolean parameter that defines whether to take a subset of 4 out of 5 most diverse sets of clustering. 
#' Can only be implemented when all 5 booleans for individual methods are set to TRUE.
#' Default is "TRUE".
#' @param SEED sets the seed of the random number generator. Setting the seed to a fixed value can
#' produce reproducible clustering results.
#'
#' @return a matrix of indiviudal clustering results is output, where each row represents the cluster results of each method.
#'
#' @author Ruth Huh <rhuh@live.unc.edu>, Yuchen Yang <yyuchen@email.unc.edu>, Yun Li <yunli@med.unc.edu>
#' @references Huh, R., Yang, Y., Jiang, Y., Shen, Y., Li, Y. (2020) SAME-clustering: Single-cell Aggregated Clustering via Mixture Model Ensemble. Nucleic Acids Research, 48: 86-95
#' @examples
#' # Load the example data data_SAME
#' data("data_SAME")
#'
#' # Zheng dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags=data_SAME$Zheng.expr, SEED=123)
#' @importFrom cidr adjustedRandIndex
#' @export
individual_clustering <- function(inputTags, mt_filter = TRUE, mt.pattern = "^MT-", mt.cutoff = 0.1, percent_dropout = 10,
                                  SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL,
                                  Seurat = TRUE, nGene_filter = TRUE, low.genes = 200, high.genes = 8000, nPC.seurat = NULL, resolution = 0.7, 
                                  tSNE = TRUE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                                  SIMLR = TRUE, diverse = TRUE, SEED = 1){
  
  cluster_number <- NULL
  cluster_results <- NULL
  inputTags = as.matrix(inputTags)
  
  # Filter out cells that have mitochondrial genes percentage over 5%
  if (mt_filter == TRUE){
    mito.genes <- grep(pattern = mt.pattern, x = rownames(x = inputTags), value = TRUE)
    percent.mito <- Matrix::colSums(inputTags[mito.genes, ])/Matrix::colSums(inputTags)
    inputTags <- inputTags[,which(percent.mito <= mt.cutoff)]
  }
  
  ##### SC3
  if(SC3 == TRUE){
    message("Performing SC3 clustering...")
    
    sc3OUTPUT <- sc3_SAME(inputTags = inputTags, gene_filter = gene_filter, svm_num_cells = svm_num_cells, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))
  }
  
  
  ##### CIDR
  if(CIDR == TRUE){
    message("Performing CIDR clustering...")
    
    cidrOUTPUT <- cidr_SAME(inputTags = inputTags, percent_dropout = percent_dropout, nPC.cidr = nPC.cidr, SEED = SEED)
    
    if(is.null(nPC.cidr)) {
      nPC.cidr <- cidrOUTPUT@nPC
    }
    
    cluster_results <- rbind(cluster_results, matrix(c(cidrOUTPUT@clusters), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number,  cidrOUTPUT@nCluster)
  }
  
  
  ##### Seurat
  if (Seurat == TRUE){
    message("Performing Seurat clustering...")
    
    if(is.null(nPC.seurat)) {
      nPC.seurat <- nPC.cidr
    }
    
    seurat_output <- seurat_SAME(inputTags = inputTags, nGene_filter = nGene_filter, low.genes = low.genes, high.genes = high.genes, 
                                 nPC.seurat = nPC.seurat, resolution = resolution, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(seurat_output), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(!is.na(seurat_output)))
  }
  
  
  ##### tSNE+kmeans
  if(tSNE == TRUE){
    message("Performing tSNE + k-means clustering...")
    
    ### Dimensionality reduction by Rtsne
    if(length(inputTags[1,]) < tsne_min_cells) {
      perplexity = tsne_min_perplexity
    }
    
    tsne_kmeansOUTPUT <- tSNE_kmeans_SAME(inputTags = inputTags, percent_dropout = percent_dropout, dimensions = dimensions,
                                          perplexity = perplexity, k.min = 2, k.max = max(cluster_number), var_genes = var_genes, SEED = SEED)
    cluster_results <- rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = TRUE))
    cluster_number <- c(cluster_number, max(as.numeric(tsne_kmeansOUTPUT$cluster)))
  }
  
  ##### SIMLR
  if(SIMLR == TRUE){
    message("Performing SIMLR clustering...")
    
    simlrOUTPUT <- SIMLR_SAME(inputTags = inputTags, percent_dropout = percent_dropout, k.min = 2, k.max = max(cluster_number), SEED = SEED)
    cluster_results <- rbind(cluster_results, simlrOUTPUT$y$cluster)
  }
  
  ##### Individual method selection
  if (dim(cluster_results)[1] == 5 && diverse == TRUE){
    message("Selecting clusteirng methods for ensemble...")
    
    rownames(cluster_results) <- c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")

    ARI=matrix(0,5,5)
    rownames(ARI) <- c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")
    colnames(ARI) <- c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")
    
    for(i in c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")){
      for(j in c("SC3","CIDR","Seurat","tSNE+kmeans","SIMLR")){
        ARI[i,j] <- adjustedRandIndex(unlist(cluster_results[i,]), unlist(cluster_results[j,]))
      }
    }
    m1 <- which.min(apply(ARI,1,var))
    cluster_results <- cluster_results[-m1,]
  }
  
  return(cluster_results)
}


#' @title SAMEclustering
#'
#' @description SAME (Single-cell RNA-seq Aggregated clustering via Mixture model Ensemble): Cluster ensemble for single-cell RNA-seq data
#'
#' @param Y a J*N matrix with J individual clustering methods and N cells.
#' @param MAX defines the maximum number of clusters used for cluster ensemble. Default is
#' the maximum cluster number estimated by all the single solutions.
#' @param rep defines how many times wants to run the cluster ensemble step.
#' @param SEED sets the seed of the random number generator. Setting the seed to a fixed value can produce
#' reproducible cluster ensemble results.
#'
#' @return SAFE returns a list object containing:
#' @return \itemize{
#'     \item AICcluster: optimal ensemble clustering result determined by Akaike information criterion (AIC)
#'     \item final_k_AIC: optimal cluster number determined by AIC
#'     \item BICcluster: optimal ensemble clustering result determined by Bayesian information criterion (BIC)
#'     \item final_k_BIC: optimal cluster number determined by BIC
#' }
#' @author Ruth Huh <rhuh@live.unc.edu>, Yuchen Yang <yyuchen@email.unc.edu>, Yun Li <yunli@med.unc.edu>
#' @references Ruth Huh, Yuchen Yang, Yun Li. SAME 2018
#' @examples
#' # Load the example data data_SAME
#' data("data_SAME")
#'
#' # Zheng dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags=data_SAME$Zheng.expr, SEED=123)
#'
#' # Cluster ensemble using SAME clustering:
#' cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED=123)
#'
#' # Biase dataset
#' # Run individual_clustering
#' cluster.result <- individual_clustering(inputTags = data_SAME$Biase.expr, datatype = "FPKM", seurat_min_cell = 200, resolution_min = 1.2, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED=123)
#'
#' # Cluster ensemble using SAME clustering:
#' cluster.ensemble <- SAMEclustering(Y=t(cluster.result), rep = 3, SEED=123)
#' @import Rcpp
#' @useDynLib EM
#' @export
SAMEclustering <- function(Y, MAX = NULL, rep, SEED = 1){
  LL <- NULL
  clusterresults <- NULL
  ensemblecluster <- NULL
  mm_k <- NULL
  ARI <- NULL
  AIC <- NULL
  BIC <- NULL
  numclus <- NULL
  I <- NULL
  N <- dim(Y)[1] #number of subjects
  H <- dim(Y)[2] #number of clustering methods
  set.seed(SEED)
  
  #number of clusters for each of the H clustering methods
  if(is.null(MAX)){
    MAX <- max(Y[!is.na(Y)])
  }
  
  for(M in 2:MAX){
    
    for (r in 1:rep){
      assign(paste('CR', r, sep = ''), EM(Y, M))
      if(r == 1){ CR = get(paste('CR', r, sep = '')) }
      else{ 
        if(get(paste('CR', r, sep = ''))[['Log Likelihood']] > CR[['Log Likelihood']]){
          CR <- get(paste('CR', r, sep = ''))
        }
      }
    }
    
    #if(CR$`Number of parameter`>N){break}
    i <- 1 
    if(length(LL) > 0 && (CR[['Number of clusters']] != M || CR[['Log Likelihood']] < LL[length(LL)])){ 
      
      repeat{
        CR <- EM(Y, M)
        i <- i + 1
        if((CR[['Number of clusters']] == M && CR[['Log Likelihood']] > LL[length(LL)]) || i > 3){break}
      }
    }
    #number of clusters
    I <- c(I, i)
    if(length(LL) > 0 && (CR[['Number of clusters']] != M || CR[['Log Likelihood']] < LL[length(LL)])){break}
    numclus <- c(numclus, CR[['Number of clusters']])
    clusterresults <- c(clusterresults, list(CR[['Cluster results']])) #clusters
    AIC <- c(AIC, CR[['AIC']])
    LL <- c(LL, CR[['Log Likelihood']])
    BIC <- c(BIC, CR[['BIC']])
    final_k_AIC <- numclus[which(AIC == min(AIC))]
    final_k_BIC <- numclus[which(BIC == min(BIC))]
  }
  
  BICresultsummary <- paste("The final optimal k is", numclus[which(BIC == min(BIC))], "which has the lowest BIC of", min(BIC), "resulting in an ARI of", ARI[which(BIC == min(BIC))])
  AICresultsummary <- paste("The final optimal k is", numclus[which(AIC == min(AIC))], "which has the lowest AIC of", min(AIC), "resulting in an ARI of", ARI[which(AIC == min(AIC))])
  BICcluster <- clusterresults[[which(BIC == min(BIC))]]
  AICcluster <- clusterresults[[which(AIC == min(AIC))]]
  result <- list(AIC_result_summary = AICresultsummary, AICcluster = AICcluster, final_k_AIC = final_k_AIC, BIC_result_summary = BICresultsummary, BICcluster = BICcluster, final_k_BIC = final_k_BIC)
  return(result)
}


