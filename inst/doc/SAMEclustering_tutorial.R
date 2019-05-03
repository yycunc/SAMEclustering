## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE---------------
library(knitr)
opts_chunk$set(echo = TRUE)

## ----setup for Zheng dataset, warning=FALSE-------------------------
# Setup the input expression matrix
library("SAMEclustering")
data("data_SAME")

## ----message=FALSE---------------------------------------------------------
dim(data_SAME$Zheng.expr)
data_SAME$Zheng.expr[1:5, 1:5]

## ----results='hide', fig.show="hide", warning=FALSE------------------------
# Perform individual clustering
cluster.result <- individual_clustering(inputTags = data_SAME$Zheng.expr, datatype = "count", mt_filter = FALSE, percent_dropout = 10, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, SEED = 123)
## ----message=FALSE---------------------------------------------------------
cluster.result[1:4, 1:10]

## ----results='hide'--------------------------------------------------------
# cluster ensemble
cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = 123)

## ----message=FALSE---------------------------------------------------------
# ensemble results,
cluster.ensemble

## ----ARI calculation-------------------------------------------------------
library(cidr)

# Cell labels of ground truth
head(data_SAME$Zheng.celltype)

# Calculating ARI for cluster ensemble
adjustedRandIndex(cluster.ensemble$AICcluster, data_SAME$Zheng.celltype)


# Biase dataset

## ----setup for Biase dataset-----------------------------------------------
# Setup the input expression matrix
dim(data_SAME$Biase.expr.expr)

data_SAME$Biase.expr[1:5, 1:5]

## ----results='hide', fig.show="hide", warning=FALSE------------------------
# Perform individual clustering
cluster.result <- individual_clustering(inputTags = data_SAME$Biase.expr, datatype = "FPKM",  mt_filter = FALSE, percent_dropout = 10, SC3 = TRUE, gene_filter = FALSE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nPC.seurat = NULL, seurat_min_cell = 200, resolution_min = 1.2, tSNE = TRUE, dimensions = 2, tsne_min_cells = 200, tsne_min_perplexity = 10, SIMLR = TRUE, diverse = TRUE, SEED = 123)

## ----results='hide', message=FALSE-----------------------------------------
# cluster ensemble
cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = 123)

## ----message=FALSE--------------------------------------------------------
cluster.ensemble

## ----ARI calculation for Biase dataset------------------------------------
library(cidr)

# Cell labels of ground truth
head(data_SAME$Biase.celltype)

# Calculating ARI for cluster ensemble
adjustedRandIndex(cluster.ensemble$AICcluster, data_SAME$Biase.celltype)
