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
keepX.constraint = NULL,
keepY.constraint = NULL,
keepX = rep(ncol(X), ncomp),
keepY = rep(ncol(Y), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE)

}	

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{Outcome. Numeric vector or matrix of responses (for multi-response models)}
\item{ncomp}{Number of components to include in the model. Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}
\item{study}{grouping factor indicating which samples are from the same study}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components}
\item{keepY.constraint}{A list containing which variables of Y are to be kept on each of the first PLS-components}
\item{keepX}{numeric vector of length \code{ncomp -length(keepX.constraint)}, the number of variables
to keep in the \eqn{X}-loadings. By default all variables are kept in the model.}
\item{keepY}{numeric vector of length \code{ncomp-length(keepY.constraint)}, the number of variables
to keep in the \eqn{Y}-loadings. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, \code{X} and \code{Y} are standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}

}


\details{
\code{mint.spls} function fits a vertical sPLS models with \eqn{1, \ldots ,}\code{ncomp} components. The grouping factor is specified with \code{study}. We advise more than 3 samples per study where all outcome categories are represented.

Multi-response models are fully supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by casewise deletion in the \code{mint.spls} function without having to delete the rows with missing data.


The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).


Variable selection is performed on each component of \code{X} via both input parameter \code{keepX.constraint} and \code{keepX} as follows. \code{keepX.constraint[[j]]} is a vector of variable names that indicates which variables to keep on component \code{j}; similarly, \code{keepX[k]} is a single positive value that indicates the number of variables to keep on component \code{length(keepX.constraint)+k}. See examples below.

Variable selection is performed on each component of \code{Y} via both input parameter \code{keepY.constraint} and \code{keepY} as follows. \code{keepY.constraint[[j]]} is a vector of variables names that indicates which variables to keep on component \code{j}  of block \code{Y}; similarly, \code{keepY[k]} is a single positive value that indicates the number of variables to keep on component \code{length(keepY.constraint)+k} of \code{Y}. See examples below.

Note that both \code{keepX.constraint} and \code{keepY.constraint} must have same length.

Useful graphical outputs are available, e.g. \code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}.
}


\value{
\code{mint.spls} returns an object of class \code{"mint.spls","spls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
\item{keepX}{Number of variables used to build each component of X}
\item{keepY}{Number of variables used to build each component of Y}
\item{keepX.constraint}{list indicating which variables where constrained to be used to build the first component of X}
\item{keepY.constraint}{list indicating which variables where constrained to be used to build the first component of Y}
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
# need to unmap Y for an unsupervised analysis
Y = unmap(stemcells$celltype)

# -- feature selection with keepX
res = mint.spls(X = stemcells$gene, Y = Y, ncomp = 3, keepX = c(10, 5, 15),
    study = stemcells$study)

plotIndiv(res)
#plot study-specific outputs for all studies
plotIndiv(res, study = levels(stemcells$study), col = as.numeric(stemcells$celltype))

#plot study-specific outputs for study "2"
plotIndiv(res, study = "2", col = as.numeric(stemcells$celltype))

# -- feature selection with keepX.constraint and keepX
# we force mint.spls to build the PLS-component 1 and 2 with only specific genes
keepX.constraint = list(comp1 = c("ENSG00000138756", "ENSG00000101470"),
                        comp2 = c("ENSG00000120049", "ENSG00000138685", "ENSG00000204248", "ENSG00000103260"))

# we select 10 and 15 variables on component 3 and 4
res = mint.spls(X = stemcells$gene, Y = Y, ncomp = 4, study = stemcells$study,
    keepX.constraint = keepX.constraint, keepX = c(10,15))

}

\keyword{regression}
\keyword{multivariate}
