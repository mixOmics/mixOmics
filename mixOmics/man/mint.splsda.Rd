\name{mint.splsda}
\encoding{latin1}
\alias{mint.splsda}

\title{Partial Least Squares (PLS) Regression}

\description{Function to perform Partial Least Squares (PLS) regression.

}

\usage{
mint.splsda(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint=NULL,
keepX=rep(ncol(X), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE)

}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{a factor or a class vector for the discrete outcome.}
\item{ncomp}{Number of components to include in the model (see Details). Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"} or \code{"canonical"}. See Details.}
\item{study}{grouping factor indicating which samples are from the same study}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components}
\item{keepX}{numeric vector of length \code{ncomp}, the number of variables
to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}

}


\details{
\code{mint.splsda} function fits a vertical sPLS-DA models with \eqn{1, \ldots ,}\code{ncomp} components. The grouping factor is specified with \code{study}. We advise more than 3 samples per study where all outcome categories are represented.

\code{X} can contain missing values. Missing values are handled by casewise deletion in the \code{mint.splsda} function without having to delete the rows with missing data.


The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).


Variable selection is performed on each component of \code{X} via both input parameter \code{keepX.constraint} and \code{keepX} as follows. \code{keepX.constraint[[j]]} is a vector of variable names that indicates which variables to keep on component \code{j}; similarly, \code{keepX[k]} is a single positive value that indicates the number of variables to keep on component \code{length(keepX.constraint)+k}. See examples below.

Useful graphical outputs are available, e.g. \code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}.
}


\value{
\code{mint.splsda} returns an object of class \code{"mint.splsda", "splsda"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ind.mat}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of X}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component of X}
\item{variates}{list containing the variates of X - global variates.}
\item{loadings}{list containing the estimated loadings for the variates - global loadings.}
\item{variates.partial}{list containing the variates of X relative to each study - partial variates.}
\item{loadings.partial}{list containing the estimated loadings for the partial variates - partial loadings.}
\item{names}{list containing the names to be used for individuals and variables.}
\item{nzv}{list containing the zero- or near-zero predictors information.}
\item{iter}{Number of iterations of the algorthm for each component}
\item{explained_variance}{Percentage of explained variance for each component and each study}

}

\references{
Rohart F. et al (2016, in prep). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms.
}

\author{Florian Rohart}

\seealso{\code{\link{spls}},
\code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{tune.mint.splsda}} and http://www.mixOmics.org for more details.}


\examples{
data(stemcells)

# -- feature selection with keepX
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3, keepX = c(10, 5, 15),
study = stemcells$study)

plotIndiv(res)
#plot study-specific outputs for all studies
plotIndiv(res, study = "all.partial")

#plot study-specific outputs for study "2"
plotIndiv(res, study = "2")

# -- feature selection with keepX.constraint and keepX
# we force mint.spls to build the PLS-component 1 and 2 with only specific genes
keepX.constraint = list(comp1 = c("ENSG00000138756", "ENSG00000101470"),
comp2 = c("ENSG00000120049", "ENSG00000138685", "ENSG00000204248", "ENSG00000103260"))

# we select 10 and 15 variables on component 3 and 4
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 4, study = stemcells$study,
keepX.constraint = keepX.constraint, keepX = c(10,15))

}

\keyword{regression}
\keyword{multivariate}
