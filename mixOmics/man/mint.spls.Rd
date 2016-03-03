\name{mint.spls}
\encoding{latin1}
\alias{mint.spls}

\title{Partial Least Squares (PLS) Regression}

\description{Function to perform Partial Least Squares (PLS) regression.

}

\usage{
mint.spls(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint=NULL,
keepY.constraint=NULL,
keepX=rep(ncol(X), ncomp),
keepY=rep(ncol(Y), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE)

}	

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{Outcome. Numeric vector or matrix of responses (for multi-response models)}
\item{ncomp}{Number of components to include in the model (see Details). Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"} or \code{"canonical"}. See Details.}
\item{study}{grouping factor indicating which samples are from the same study}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components}
\item{keepY.constraint}{A list containing which variables of Y are to be kept on each of the first PLS-components}
\item{keepX}{numeric vector of length \code{ncomp}, the number of variables
to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
\item{keepY}{numeric vector of length \code{ncomp}, the number of variables
to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}

}


\details{
\code{pls} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp} components.
Multi-response models are fully supported. The \code{X} and \code{Y} datasets
can contain missing values.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS 
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis 
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS 
algorithm \code{("classic")} (see References).

The estimation of the missing values can be performed 
by the reconstitution of the data matrix using the \code{nipals} function. Otherwise, missing 
values are handled by casewise deletion in the \code{pls} function without having to 
delete the rows with missing data.  
}

\value{
\code{pls} returns an object of class \code{"pls"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized original response vector or matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{mode}{the algorithm used to fit the model.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the \eqn{X} and \eqn{Y} variates.}
  \item{loadings}{list containing the estimated loadings for the variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
   \item{tol}{the tolerance used in the iterative algorithm, used for subsequent S3 methods}
  \item{max.iter}{the maximum number of iterations, used for subsequent S3 methods}
    \item{iter}{Number of iterations of the algorthm for each component}

}

\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares. 
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.
}

\author{Sebastien Dejean and Ignacio Gonzalez and Kim-Anh Le Cao.}

\seealso{\code{\link{spls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for more details.}

\examples{
data(celltype)
# need to unmap Y for an unsupervised analysis
Y = unmap(celltype)
res=mint.spls(X=gene,Y=Y,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=study)
plotIndiv(res,ind.names=FALSE)
plotIndiv(res,ind.names=FALSE,study=levels(study),col=as.numeric(celltype)) #plot study-specific outputs for all studies
plotIndiv(res,ind.names=TRUE,study=levels(study),col=as.numeric(celltype))
plotIndiv(res,ind.names=FALSE,study="2",col=as.numeric(celltype)) #plot study-specific outputs for study "2"

}

\keyword{regression}
\keyword{multivariate}
