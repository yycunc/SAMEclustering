\name{individual_clustering}
\alias{individual_clustering}
\alias{sc3_SAME}
\alias{cidr_SAME}
\alias{seurat_SAME}
\alias{tSNE_kmeans_SAME}
\alias{SIMLR_SAME}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
individual_clustering
}
\description{
This function performs single-cell clustering using five state-of-the-art methods, SC3, CIDR, Seurat, tSNE+kmeans and SIMLR.
}
\usage{
individual_clustering(inputTags, datatype = "count", mt_filter = FALSE, low.mt = -Inf, high.mt = 0.05, 
percent_dropout = 10, SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL, 
Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, seurat_min_cell = 200, resolution_min = 1.2, 
tSNE = TRUE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL, 
SIMLR = TRUE, diverse = TRUE, SEED = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputTags}{
a G*N matrix with G genes and N cells.
}
  \item{datatype}{
defines the type of data, which could be "count", "CPM", "RPKM" and "FPKM". Default is "count".
}
  \item{mt_filter}{
is a boolean variable that defines whether to filter outlier cells according to mitochondrial gene percentage. Default is "FALSE".
}
  \item{low.mt}{
defines a low cutoff of mitochondrial percentage (Default is -Inf) that cells having lower percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}.
}
  \item{high.mt}{
defines a high cutoff of mitochondrial percentage (Default is 0.05) that cells having higher percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}
}
  \item{percent_dropout}{
defines a low cutoff of gene percentage that genes expressed in less than \code{percent_dropout}\% or more than (100 - \code{percent_dropout})\% of cells are removed for CIDR, tSNE + k-means, and SIMLR clustering. Default is 10.
}
  \item{SC3}{
a boolean variable that defines whether to cluster cells using SC3 method. Default is "TRUE".
}
  \item{gene_filter}{
is a boolean variable that defines whether to perform gene filtering before SC3 clustering, when \code{SC3 = TRUE}.
}
  \item{svm_num_cells}{
if \code{SC3 = TRUE}, then defines the mimimum number of cells above which SVM will be run.
}
  \item{CIDR}{
is a boolean parameter that defines whether to cluster cells using CIDR method. Default is "TRUE".
}
  \item{nPC.cidr}{
defines the number of principal coordinates used in CIDR clustering, when \code{CIDR = TRUE}. Default value is esimated by \code{nPC} of \code{CIDR}.
}
  \item{Seurat}{
is a boolean variable that defines whether to cluster cells using Seurat method. Default is "TRUE".
}
  \item{nPC.seurat}{
defines the number of principal components used in Seurat clustering, when \code{Seurat = TRUE}. Default is \code{nPC.seurat = nPC.cidr}.
}
  \item{resolution}{
defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}.
}
  \item{seurat_min_cell}{
defines the mimimum number of cells in input dataset below which \code{resolution} is set to 1.2, when \code{Seurat = TRUE}.
}
  \item{resolution_min}{
defines the resolution used in Seurat clustering for small dataset, when \code{Seurat == TRUE} and cell number of input file < \code{seurat_min_cell}.
}
  \item{tSNE}{
is a boolean variable that defines whether to cluster cells using t-SNE + k-means method. Default is "TRUE".
}
  \item{dimensions}{
sets the number of dimensions wanted to be retained in t-SNE step. Default is 3.
}
  \item{perplexity}{
sets the perplexity parameter for t-SNE dimension reduction. Default is 30 when number of cells \code{>=200}.
}
  \item{tsne_min_cells}{
defines the number of cells in input dataset below which \code{tsne_min_perplexity=10} would be employed for t-SNE step. Default is 200.
}
  \item{tsne_min_perplexity}{
sets the perplexity parameter of t-SNE step for small datasets (number of cells \code{<200}).
}
  \item{var_genes}{
defines the number of variable genes used by t-SNE analysis, when \code{tSNE = TRUE}.
}
  \item{SIMLR}{
is a boolean variable that defines whether to cluster cells using t-SNE + k-means method. Default is "TRUE".
}
  \item{diverse}{
is a boolean parameter that defines whether to take a subset of 4 out of 5 most diverse sets of clustering. Can only be implemented when all 5 booleans for individual methods are set to TRUE. Default is "TRUE".
}
  \item{SEED}{
sets the seed of the random number generator. Setting the seed to a fixed value can produce reproducible clustering results.
}
}
\value{
a matrix of indiviudal clustering results is output, where each row represents the cluster results of each method.
}
\references{
Ruth Huh, Yuchen Yang, Yun Li. SAME. 2017
}
\author{
Ruth Huh <rhuh@live.unc.edu>, Yuchen Yang <yyuchen@email.unc.com>, Yun Li <yunli@med.unc.edu>
}
\examples{
# Load the example data data_SAME
data("data_SAME")

# Zheng dataset
# Run individual_clustering
cluster.result <- individual_clustering(inputTags=data_SAME$Zheng.expr, SEED=123)

# Biase dataset
# Run individual_clustering
cluster.result <- individual_clustering(inputTags = data_SAME$Biase.expr, datatype = "FPKM", seurat_min_cell = 200, resolution_min = 1.2, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED=123)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
