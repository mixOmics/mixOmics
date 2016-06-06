\name{plotDiablo}
\encoding{latin1}
\alias{plotDiablo}

\title{Graphical output for the DIABLO framework}

\description{Function to perform Horizontal sparse Partial Least Squares-  Discriminant Analysis (sPLS-DA) regression.

}

\usage{
plotDiablo(x, ncomp = 1, groupOrder, ...)

}

\arguments{
\item{x}{object of class inheriting from \code{"block.splsda"}.}
\item{ncomp}{Which component to plot. Has to be lower than the minimum of \code{object$ncomp}}
\item{groupOrder}{Ordering of the Y-levels on the graphical outputs}
\item{\dots}{not used}
}


\details{

}




\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery.
}

\author{Amrit Singh}

\seealso{\code{\link{block.splsda}} and http://www.mixOmics.org for more details.}

\examples{

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)


nutrimouse.sgccda <- wrapper.sgccda(X=data,
Y = Y,
design = design,
keepX = list(gene=c(10,10), lipid=c(15,15)),
ncomp = c(2, 2),
scheme = "centroid",
verbose = FALSE,
bias = FALSE)

plotDiablo(nutrimouse.sgccda)
plotDiablo(nutrimouse.sgccda, groupOrder = levels(Y)[nlevels(Y):1])


}

\keyword{regression}
\keyword{multivariate}
