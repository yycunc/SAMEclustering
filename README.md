# SAMEclustering
SAME (Single-cell RNA-seq Aggregated clustering via Mixture model Ensemble): Cluster ensemble for single-cell RNA-seq data

Although several methods have been recently developed for clustering cell types using single-cell RNA-seq (scRNA-Seq) data, they utilize different characteristics of data and yield varying results in terms of both the number of clusters and actual cluster assignments. Here, we present SAME-clustering, Single-cell RNA-seq Aggregated clustering via Mixture model Ensemble, a flexible, accurate and robust method for clustering scRNA-Seq data. SAMEclustering takes as input, results from multiple clustering methods, to build one consensus solution. SAME-clustering currently embeds five state-of-the-art methods, SC3, CIDR, Seurat, t-SNE + *k*-means and SIMLR; and ensembles solutions from four diverse sets of five individual methods using EM algorithms.

SAMEclustering is maintained by Yuchen Yang [yyuchen@email.unc.edu]

## News and Updates
Dec 5, 2018
* Version 0.99.0 released
  + First offical release
  + It can work on Windows, Mac and Linux platforms

## Installation
You can install SAMEclustering from github with:
```{r install}
install.packages("devtools)

devtools::install_github("yycunc/SAMEclustering")
```

## SAMEclustering Examples
Here we will provide examples using two datasets: one from Zheng *et al.*, (Nature Communications, 2016) and the other from Biase *et al.*, (Genome Research, 2014). Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells. The original data can be downloaded from [10X GENOMICS website](https://support.10xgenomics.com/single-cell-gene-expression/datasets). The Biase dataset has 49 mouse embryo cells, which were sequenced by SMART-Seq and can be found at [NCBI GEO:GSE57249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249).

### Load the data
```{r setup for Zheng dataset}
library("SAMEclustering")
data("data_SAME")
```

### Zheng dataset
#### Setup the input expression matrix
```{r setup for Zheng dataset}
dim(data_SAME$Zheng.expr)

data_SAME$Zheng.expr[1:5, 1:5]
```

#### Perform individual clustering
Here we perform single-cell clustering using four popular methods, SC3, CIDR, Seurat, t-SNE + *k*-means and SIMLR. Genes expressed in less than 10% or more than 90% of cells are removed for CIDR, tSNE + k-means and SIMLR clustering. To improve the performance of cluster ensemble, we take a subset of 4 out of 5 most diverse sets of clustering.

```{r individual clustering for Baron_human4 dataset, results='hide', fig.show="hide", warning=FALSE}
cluster.result <- individual_clustering(inputTags = data_SAME$Zheng.expr, datatype = "count", 
percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, 
resolution = 0.9, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, SEED = 123)
```

The function *indiviual_clustering* will output a matrix, where each row represents the cluster results of each method, and each column represents a cell. User can also extend SAMEclustering to other scRNA-seq clustering methods, by putting all clustering results into a M * N matrix with M clustering methods and N cells.

```{r, message=FALSE}
cluster.result[1:4, 1:10]
```

#### Cluster ensemble

Using the clustering results generated in last step, we perform cluster ensemble using EM algorithm.

```{r cluster ensemble for Zheng dataset, results='hide'}
cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = 123)
```

Function *SAMEclustering* will output a list for optimal clustering ensemble and clustering number based on AIC and BIC index, respectively.

```{r ensemble results for Zheng dataset, message=FALSE}
cluster.ensemble
```

We can compare the clustering results to the true labels using the Adjusted Rand Index (ARI)

```{r ARI calculation for Zheng dataset}
library(cidr)

# Cell labels of ground truth
head(data_SAME$Zheng.celltype)

# Calculating ARI for cluster ensemble
adjustedRandIndex(cluster.ensemble$optimal_clustering, data_SAME$Zheng.celltype)
```

### Biase dataset

#### Setup the input expression matrix
```{r setup for Biase dataset}
dim(data_SAME$Biase.expr.expr)

data_SAME$Biase.expr[1:5, 1:5]
```

#### Perform individual clustering

Here we perform single-cell clustering using four popular methods, SC3, CIDR, Seurat, t-SNE + *k*-means and SIMLR. Genes expressed in less than 10% or more than 90% of cells are removed for CIDR, tSNE + k-means and SIMLR clustering. Since there are only 49 cells in Biase dataset, the resolution parameter is set to 1.2 according to our benchmarking results. Four out of five most diverse sets of clustering were taken out for downstream cluster ensemble.

```{r individual clustering for Biase dataset, results='hide', fig.show="hide", warning=FALSE}
cluster.result <- individual_clustering(inputTags = data_SAME$Biase.expr, datatype = "FPKM",  
percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, 
seurat_min_cell = 200, resolution_min = 1.2, tSNE = TRUE, dimensions = 2, tsne_min_cells = 200, 
tsne_min_perplexity = 10, SIMLR = TRUE, diverse = TRUE, SEED = 123)
```

#### Cluster ensemble

Using the clustering results, we perform cluster ensemble using EM algorithm.

```{r cluster ensemble for Biase dataset, results='hide', message=FALSE}
cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = 123)
```

```{r ensemble results for Biase dataset, message=FALSE}
cluster.ensemble
```
Compare the cluster ensemble results to the true labels.

```{r ARI calculation for Biase dataset}
# Cell labels of ground truth
head(data_SAME$Biase.celltype)

# Calculating ARI for cluster ensemble
adjustedRandIndex(cluster.ensemble$optimal_clustering, data_SAME$Biase.celltype)
```


