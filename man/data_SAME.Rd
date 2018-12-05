\name{data_SAME}
\alias{data_SAME}
\docType{data}
\title{
data_SAME
}
\description{
An expression_matrix.csv (32738 genes and 500 cells of three cell types,
cd56_NK cells, b cells and regulatory T cells from 10X Genomices
(Zheng et al. 2016)), and celltype.info.csv (List of cell type and cluster
information for the 500 cells).
}
\usage{data("data_SAME")}
\format{
  An object of class \code{list} of length 4.
  
}
\source{
data_SAME contains two datasets: one is from Zheng et al., (Nature Communications, 2016). Zheng dataset contains 500 human peripheral blood mononuclear cells (PBMCs) sequenced using GemCode platform, which consists of three cell types, CD56+ natural killer cells, CD19+ B cells and CD4+/CD25+ regulatory T cells. The original data can be downloaded from 10X GENOMICS website. The other is from Biase et al., (Genome Research, 2014). The Biase dataset has 49 mouse embryo cells, which were sequenced by SMART-Seq and can be found at NCBI GEO:GSE57249.
}
\references{
Ruth Huh, Yuchen Yang, Yun Li. SAME 2018
}
\examples{
data(data_SAME)
## maybe str(data_SAME) ; plot(data_SAME) ...
}
\keyword{datasets}
