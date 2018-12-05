\name{SAMEclustering}
\alias{SAMEclustering}
\alias{EM}
\alias{rdirichletC}
\alias{tableC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SAMEclustering
}
\description{
SAME (Single-cell RNA-seq Aggregated clustering via Mixture model Ensemble): Cluster ensemble for single-cell RNA-seq data
}
\usage{
SAMEclustering(Y, MAX = NULL, rep, SEED = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
a J*N matrix with J individual clustering methods and N cells.
}
  \item{MAX}{
defines the maximum number of clusters used for cluster ensemble. Default is the maximum cluster number estimated by all the single solutions.
}
  \item{rep}{
defines how many times wants to run the cluster ensemble step.
}
  \item{SEED}{
sets the seed of the random number generator. Setting the seed to a fixed value can produce reproducible cluster ensemble results.
}
}
\value{
SAMEclustering returns a list object containing:

\itemize{
     \item AICcluster: optimal ensemble clustering result determined by Akaike information criterion (AIC)
     \item final_k_AIC: optimal cluster number determined by AIC
     \item BICcluster: optimal ensemble clustering result determined by Bayesian information criterion (BIC)
     \item final_k_BIC: optimal cluster number determined by BIC
}
}
\references{
Ruth Huh, Yuchen Yang, Yun Li. SAME 2018
}
\author{
Ruth Huh <rhuh@live.unc.edu>, Yuchen Yang <yyuchen@email.unc.edu>, Yun Li <yunli@med.unc.edu>
}
\examples{
# Load the example data data_SAME
data("data_SAME")

# Zheng dataset
# Run individual_clustering
cluster.result <- individual_clustering(inputTags=data_SAME$Zheng.expr, SEED=123)

# Cluster ensemble using SAME clustering:
cluster.ensemble <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED=123)

# Biase dataset
# Run individual_clustering
cluster.result <- individual_clustering(inputTags = data_SAME$Biase.expr, datatype = "FPKM", seurat_min_cell = 200, resolution_min = 1.2, tsne_min_cells = 200, tsne_min_perplexity = 10, SEED=123)

# Cluster ensemble using SAME clustering:
cluster.ensemble <- SAMEclustering(Y=t(cluster.result), rep = 3, SEED=123)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
