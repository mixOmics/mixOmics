\name{block.spls}
\encoding{latin1}
\alias{block.spls}

\title{Horizontal sparse Partial Least Squares (sPLS) Regression}

\description{Integration of multilple data sets matching on the same samples while selection features in each data set, ie. Horizontal sparse Partial Least Squares (sPLS) regression.}


\usage{
block.spls(X,
Y,
indY,
ncomp = rep(2,length(X)),
keepX.constraint,
keepY.constraint,
keepX,
keepY,
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
\item{indY}{To supply if Y is missing, indicates the position of the outcome in the list \code{X}}
\item{ncomp}{Numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.}
\item{keepX.constraint}{A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components}
\item{keepY.constraint}{Only used if Y is provided, otherwise extracted from keepX.constraint. A list containing which variables of Y are to be kept on each of the first PLS-components}
\item{keepX}{ A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).}
\item{keepY}{Only used if Y is provided. Each entry keepY[i] is the number of Y-variables kept in the model on the last components.}
\item{design}{numeric matrix of size (number of blocks) x (number of blocks) with only 0 or 1 values. A value of 1 (0) indicates a relationship (no relationship) between the blocks to be modelled. If \code{Y} is provided instead of \code{indY}, the \code{design} matrix is changed to include relationships to \code{Y}. }
\item{scheme}{Either "horst", "factorial" or "centroid". Default = \code{centroid}, see reference paper.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details. Default = \code{regression}.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{bias}{boleean. A logical value for biaised or unbiaised estimator of the var/cov. Default = \code{FALSE}.}
\item{init}{Mode of initialization use in the algorithm, either by Singular Value Decompostion of the product of each block of X with Y ("svd") or each block independently ("svd.single"). Default = \code{svd}.}
\item{tol}{Convergence stopping value.}
\item{verbose}{if set to \code{TRUE}, reports progress on computing.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}

}


\details{
\code{block.spls} function fits a horizontal sPLS models with a specific number of component per block (\code{ncomp[i]} components for each block \code{X[i]}).
An outcome needs to be provided, either by \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
Multi-response models are fully supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by casewise deletion in the \code{block.pls} function without having to delete the rows with missing data.

The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).

Variable selection is performed on each component for each block of \code{X} via both input parameter \code{keepX.constraint} and \code{keepX}. \code{keepX.constraint$'blocki'[[j]]} is a vector of variables names that indicates which variables to keep on component \code{j}  of block \code{X$'blocki'}; similarly, \code{keepX[[i]][k]} is a single positive value that indicates the number of variables to keep on component \code{j+k} of block \code{X[[i]]}.
if \code{Y} is given, variable selection is performed on block \code{Y} via both input parameter \code{keepY.constraint} and \code{keepY}. \code{keepY.constraint[[j]]} is a vector of variables names that indicates which variables to keep on component \code{j}  of block \code{Y}; similarly, \code{keepY[k]} is a single positive value that indicates the number of variables to keep on component \code{j+k} of block \code{Y}.
Note that if \code{Y} is missing and \code{indY} is provided, then variable selection on \code{Y} is performed via both input parameter \code{keepX.constraint} and \code{keepX}.

}

\value{
\code{block.spls} returns an object of class \code{"block.spls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{indY}{the position of the outcome Y in the output list X.}
\item{ncomp}{the number of components included in the model.}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of each block}
\item{keepY}{Number of variables used to build each component of Y}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component for each block}
\item{keepY.constraint}{list indicating which variables where constrained to be used to build the first component for Y}
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

Tenenhaus A. and Tenenhaus M., (2011), Regularized
Generalized Canonical Correlation Analysis,
Psychometrika, Vol. 76, Nr 2, pp 257-284.

Schafer J. and Strimmer K., (2005), A shrinkage approach
to large-scale covariance matrix estimation and
implications for functional genomics. Statist. Appl.
Genet. Mol. Biol. 4:32.
}

\author{Florian Rohart}

\seealso{\code{\link{plotIndiv}}, \code{\link{plotArrow}}, \code{\link{plotLoadings}}, \code{\link{plotVar}},
\code{\link{predict}}, \code{\link{perf}}, \code{\link{selectVar}} and http://www.mixOmics.org for more details.}


\examples{

data(nutrimouse)
# need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))

# indY indicates where the outcome Y is in the list X
res = block.spls(X = data, indY = 3, keepX = list(gene = c(10,5,15), lipid = c(3,2)),
    ncomp = c(3,3,3))
res$keepX
selectVar(res, comp = 1)
plotIndiv(res)
plotVar(res)

}

\keyword{regression}
\keyword{multivariate}
