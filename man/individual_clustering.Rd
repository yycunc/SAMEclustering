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
individual_clustering(inputTags, mt_filter = TRUE, mt.pattern = "^MT-", mt.cutoff = 0.1, 
percent_dropout = 10, SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, CIDR = TRUE, nPC.cidr = NULL, 
Seurat = TRUE, nGene_filter = TRUE, low.genes = 200, high.genes = 8000, nPC.seurat = NULL, resolution = 0.7, 
tSNE = TRUE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL, 
SIMLR = TRUE, diverse = TRUE, SEED = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputTags}{
a G*N matrix with G genes and N cells.
}
  \item{mt_filter}{
is a boolean variable that defines whether to filter outlier cells according to mitochondrial gene percentage. Default is "TRUE".
}
  \item{mt.pattern}{
defines the pattern of mitochondrial gene names in the data, for example, \code{mt.pattern = "^MT-"} for human and \code{mt.pattern = "^mt-"} for mouse. Only applied when \code{mt_filter = TRUE}. Default is \code{mt.pattern = "^MT-"}.
}
  \item{mt.cutoff}{
defines a high cutoff of mitochondrial percentage (Default is 0.1) that cells having higher percentage of mitochondrial gene are filtered out, when \code{mt_filter = TRUE}
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
  \item{nGene_filter}{
is a boolean variable that defines whether to filter outlier cells according to unique gene count before Seurat clustering. Default is "TRUE".
}
  \item{low.genes}{
defines a low cutoff of unique gene counts (Default is 200) that cells expressing less than 200 genes are filtered out, when \code{nGene_filter = TRUE}.
}
  \item{high.genes}{
defines a high cutoff of unique gene counts (Default is 8000) that cells expressing more than 8000 genes are filtered out, when \code{nGene_filter = TRUE}.
}
  \item{nPC.seurat}{
defines the number of principal components used in Seurat clustering, when \code{Seurat = TRUE}. Default is \code{nPC.seurat = nPC.cidr}.
}
  \item{resolution}{
defines the value of resolution used in Seurat clustering, when \code{Seurat = TRUE}. Default is \code{resolution = 0.7}.
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
Huh, R., Yang, Y., Jiang, Y., Shen, Y., Li, Y. (2020) SAME-clustering: Single-cell Aggregated Clustering via Mixture Model Ensemble. Nucleic Acid Research, 48: 86-95
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
