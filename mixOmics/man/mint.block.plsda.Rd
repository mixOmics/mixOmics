\name{mint.block.plsda}
\encoding{latin1}
\alias{mint.block.plsda}

\title{Horizontal and Vertical Partial Least Squares - Discriminant Analysis (PLS-DA) Regression}

\description{Function to perform Horizontal and Vertical Partial Least Squares - Discriminant Analysis (PLS-DA) regression.

}

\usage{
mint.block.plsda(X,
Y,
indY,
study,
ncomp=rep(2,length(X)),
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
\item{study}{grouping factor indicating which samples are from the same study}
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

res=mint.splsda(X=gene,Y=celltype,ncomp=3,near.zero.var=FALSE,keepX=c(10,5,15),study=study)

}

\keyword{regression}
\keyword{multivariate}
