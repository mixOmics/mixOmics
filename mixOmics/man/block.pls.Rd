\name{block.pls}
\encoding{latin1}
\alias{block.pls}

\title{Horizontal Partial Least Squares (PLS) Regression}

\description{Function to perform Horizontal Partial Least Squares (PLS) regression.

}

\usage{
block.pls(X,
Y,
indY,
ncomp = rep(2,length(X)),
design,
scheme,
mode,
scale = TRUE,
bias,
init ,
tol = 1e-06,
verbose,
max.iter = 500,
near.zero.var = FALSE)
}

\arguments{
\item{X}{A list of data sets (called 'blocks') matching on the same samples. Data in the list should be arranged in samples x variables, with samples order matching in all data sets.}
\item{Y}{Outcome. Numeric vector or matrix of responses (for multi-response models)}
\item{indY}{To supply if Y is missing, indicates the position of the outcome in the list X}
\item{ncomp}{Numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.}
\item{design}{numeric matrix of size (number of blocks) x (number of blocks) with only 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled. If \code{Y} is provided instead of \code{indY}, the \code{design} matrix is changed to include relationships to \code{Y}. }
\item{scheme}{Either "horst", "factorial" or "centroid" (Default: "centroid"), see reference paper.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{bias}{boleean. A logical value for biaised or unbiaised estimator of the var/cov (defaults to FALSE).}
\item{init}{Mode of initialization use in the algorithm, either by Singular Value Decompostion of the product of each block of X with Y ("svd") or each block independently ("svd.single") . Default to "svd".}
\item{tol}{Convergence stopping value.}
\item{verbose}{if set to \code{TRUE}, reports progress on computing.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}
}


\details{
\code{block.pls} function fits a horizontal PLS models with a specific number of component per block (\eqn{1, \ldots ,}\code{ncomp[i]} components for each block \code{X[i]}).

An outcome needs to be provided, either by \code{Y} or by its position \code{indY} in the list of blocks \code{X}.

Multi-response models are fully supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by casewise deletion in the \code{block.pls} function without having to delete the rows with missing data.


The type of algorithm to use is specified with the \code{mode} argument. Four PLS 
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis 
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS 
algorithm \code{("classic")} (see References).

}

\value{
\code{block.pls} returns an object of class \code{"block.pls"}, a list
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{indY}{the position of the outcome Y in the output list X.}
  \item{ncomp}{the number of components included in the model.}
  \item{mode}{the algorithm used to fit the model.}
  \item{variates}{list containing the variates of each block of X.}
  \item{loadings}{list containing the estimated loadings for the variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
    \item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each block}

}

\references{
Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris: Editions Technic.

Wold H. (1966). Estimation of principal components and related models by iterative least squares. 
In: Krishnaiah, P. R. (editors), \emph{Multivariate Analysis}. Academic Press, N.Y., 391-420.
}

\author{Florian Rohart, etc}

\seealso{\code{\link{spls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for more details.}

\examples{

data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))

# indY indicates where the outcome Y is in the list X
res=block.pls(X=data,indY=3)
plotIndiv(res)
plotVar(res)

# Y is provided
res2=block.pls( list(gene = nutrimouse$gene, lipid = nutrimouse$lipid),Y=Y,ncomp=c(2,2))
plotIndiv(res2)
plotVar(res2)


}

\keyword{regression}
\keyword{multivariate}
