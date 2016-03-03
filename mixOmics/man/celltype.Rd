\name{celltype}
\encoding{latin1}
\alias{celltype}

\docType{data}

\title{Human Stem Cells Data}

\description{
This data set contains the expression of a random subset of 1,000 genes in 167 samples from 6 independent studies and 3 cell types.
}

\usage{data(celltype)}

\format{
A list containing the following components:
  \describe{
    \item{\code{gene}}{data matrix with 167 rows and 1000 columns. Each row represents
      an experimental sample, and each column a single gene.}
	\item{\code{celltype}}{a factor indicating the cell type of each sample.}
    \item{\code{study}}{a factor indicating the study from which the sample was extracted.}
  }
}

\details{
This data set contains the expression of a random subset of 1,000 genes in 167 samples from 6 independent studies and 3 cell types.
}

\source{
The Human Breast Tumors dataset is a companion resource for the paper of
Perou \emph{et al.} (2000), and was downloaded from the Stanford Genomics Breast 
Cancer Consortium Portal 
\url{http://genome-www.stanford.edu/breast_cancer/molecularportraits/download.shtml}
}

\references{
Perez-Enciso, M. and Tenenhaus, M. (2003). Prediction of clinical outcome with microarray data: 
a partial least squares discriminant analysis (PLS-DA) approach. \emph{Human Genetics} 
\bold{112}, 581-592.

Perou, C. M., Sorlie, T., Eisen, M. B., van de Rijn, M., Jeffrey, S. S., Rees, C. A., 
Pollack, J. R., Ross, D. T., Johnsen, H., Akslen, L. A., Fluge, O.,
Pergamenschikov, A., Williams, C., Zhu, S. X., Lonning, P. E., Borresen-Dale, A. L., 
Brown, P. O. and Botstein, D. (2000). Molecular portraits
of human breast tumours. \emph{Nature} \bold{406}, 747-752.
}

\keyword{datasets}
