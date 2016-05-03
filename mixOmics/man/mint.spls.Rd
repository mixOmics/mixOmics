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
\item{keepX}{numeric vector of length \code{ncomp -length(keepX.constraint)}, the number of variables
to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
\item{keepY}{numeric vector of length \code{ncomp-length(keepY.constraint)}, the number of variables
to keep in \eqn{Y}-loadings. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}

}


\details{
\code{mint.spls} function fits a vertical PLS models with (\eqn{1, \ldots ,}\code{ncomp} components.

Multi-response models are fully supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by casewise deletion in the \code{block.pls} function without having to delete the rows with missing data.


The type of algorithm to use is specified with the \code{mode} argument. Two PLS
algorithms are available: PLS regression \code{("regression")} and PLS canonical analysis
\code{("canonical")}.


Variable selection is performed on each component of \code{X} via both input parameter \code{keepX.constraint} and \code{keepX}. \code{keepX.constraint[[j]]} is a vector of variable names that indicates which variables to keep on component \code{j}; similarly, \code{keepX[k]} is a single positive value that indicates the number of variables to keep on component \code{j+k}.
Variable selection is performed on \code{Y} via both input parameter \code{keepY.constraint} and \code{keepY}. \code{keepY.constraint[[j]]} is a vector of variables names that indicates which variables to keep on component \code{j}  of block \code{Y}; similarly, \code{keepY[k]} is a single positive value that indicates the number of variables to keep on component \code{j+k} of \code{Y}.

}


\value{
\code{mint.spls} returns an object of class \code{"mint.spls","spls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{indY}{the position of the outcome Y in the output list X.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of X}
\item{keepY}{Number of variables used to build each component of Y}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component of X}
\item{keepY.constraint}{list indicating which variables where constrained to be used to build the first component of Y}
\item{variates}{list containing the variates of X.}
\item{loadings}{list containing the estimated loadings for the variates.}
\item{variates.partial}{list containing the variates of X relative to each study.}
\item{loadings.partial}{list containing the estimated loadings for the partial variates.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information.}
\item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each study}

}



\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares. 
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.
}

\author{Florian Rohart}

\seealso{\code{\link{spls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for more details.}

\examples{
data(stemcells)
# need to unmap Y for an unsupervised analysis
Y = unmap(stemcells$celltype)
res=mint.spls(X=stemcells$gene,Y=Y,ncomp=3,near.zero.var=FALSE,
    keepX=c(10,5,15),study=stemcells$study)

plotIndiv(res)
#plot study-specific outputs for all studies
plotIndiv(res,study=levels(stemcells$study),col=as.numeric(stemcells$celltype))

#plot study-specific outputs for study "2"
plotIndiv(res,study="2",col=as.numeric(stemcells$celltype))
}

\keyword{regression}
\keyword{multivariate}
