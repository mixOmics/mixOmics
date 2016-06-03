\name{block.splsda}
\encoding{latin1}
\alias{block.splsda}
\alias{wrapper.sgccda}

\title{Horizontal sparse Partial Least Squares - Discriminant Analysis (sPLS-DA) Regression}

\description{Integration of multilple data sets matching on the same samples in the aim of classifying a discrete outcome, while selection features in each data set, ie. Horizontal sparse Partial Least Squares - Discriminant Analysis (sPLS-DA) regression.}


\usage{
block.splsda(X,
Y,
indY,
ncomp=rep(2,length(X)),
keepX.constraint,
keepX,
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
\item{Y}{A factor or a class vector for the discrete outcome.}
\item{indY}{To supply if Y is missing, indicates the position of the outcome in the list \code{X}}
\item{ncomp}{Numeric vector of length the number of blocks in \code{X}. The number of components to include in the model for each block (does not necessarily need to take the same value for each block). By default set to 2 per block.}
\item{keepX.constraint}{A list of same length as X. Each entry keepX.constraint[[i]] is a list containing which variables of X[[i]] are to be kept on each of the first PLS-components}
\item{keepX}{ A vector of same length as X.  Each entry keepX[i] is the number of X[[i]]-variables kept in the model on the last components (once all keepX.constraint[[i]] are used).}
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
\code{block.splsda} function fits a horizontal sPLS-DA models with a specific number of component per block (\code{ncomp[i]} components for each block \code{X[i]}).
A factor for the discrete outcome needs to be provided, either by \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
Multi-response models are fully supported. \code{X} can contain missing values. Missing values are handled by casewise deletion in the \code{block.splsda} function without having to delete the rows with missing data.
The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).

Variable selection is performed on each component for each block of \code{X} via both input parameter \code{keepX.constraint} and \code{keepX}. \code{keepX.constraint$'blocki'[[j]]} is a vector of variables names that indicates which variables to keep on component \code{j}  of block \code{X$'blocki'}; similarly, \code{keepX[[i]][k]} is a single positive value that indicates the number of variables to keep on component \code{j+k} of block \code{X[[i]]}.
}


\value{
\code{block.splsda} returns an object of class \code{"block.splsda", "block.spls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{indY}{the position of the outcome Y in the output list X.}
\item{ncomp}{the number of components included in the model.}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of each block}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component for each block}
\item{variates}{list containing the variates of each block of X.}
\item{loadings}{list containing the estimated loadings for the variates.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information.}
\item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each block}

}


\references{
Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Le Cao K.A. (2016).
DIABLO - multi omics integration for biomarker discovery.

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
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = nutrimouse$diet)
# with this design, all blocks are connected
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
byrow = TRUE, dimnames = list(names(data), names(data)))

# indY indicates where the outcome Y is in the list X
res = block.splsda(X = data, indY = 3, keepX = list(gene = c(10,5,15), lipid = c(3,2)),
    ncomp = c(3,3,3))

res$keepX
selectVar(res, comp = 1)
plotIndiv(res)
plotVar(res)

# -- feature selection with keepX.constraint and keepX
# we force the first PLS-component of X to be build with only the variables "ENSG00000164930" and "ENSG00000044090",
# and the second PLS-component with only "ENSG00000109819"
# we decide to keep 100 variables on the third component
data(stemcells)

res = block.splsda(list(block1 = stemcells$gene), Y = stemcells$celltype, ncomp = 3,
keepX = list(block1 = c(100)),
keepX.constraint = list(block1 = list(comp1 = c("ENSG00000164930", "ENSG00000044090"), comp2 = c("ENSG00000109819"))))

res$keepX
res$keepX.constraint # position of "ENSG00000164930", "ENSG00000044090" and "ENSG00000109819" in the matrix X


}

\keyword{regression}
\keyword{multivariate}
