\name{mint.splsda}
\encoding{latin1}
\alias{mint.splsda}

\title{Vertical sparse Partial Least Squares - Discriminant Analysis (sPLS-DA) Regression}

\description{Integration of multilple data sets (or studies) matching on the same variables in the aim of classifying a discrete outcome, while selection features in each data set, ie. vertical sparse Partial Least Squares - Discriminant Analysis (sPLS-DA) regression.}

\usage{
mint.splsda(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
keepX.constraint = NULL,
keepX = rep(ncol(X), ncomp),
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE)

}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed. Concatenation of multiple data sets matching on the same variables.}
\item{Y}{a factor or a class vector for the discrete outcome.}
\item{ncomp}{Number of components to include in the model (see Details). Default to 2}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"} or \code{"canonical"}. See Details.}
\item{study}{grouping factor indicating which samples are from the same study (or data set)}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components}
\item{keepX}{numeric vector of length \code{ncomp}, the number of variables
to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}

}



\details{
\code{mint.splsda} function fits a vertical sPLS-DA models with \code{ncomp} components in which several studies matching on the same variables are integrated together in the aim of classifying a discrete outcome while selecting features. The grouping factor that indicates which samples are from which study is specified with \code{study}. We advise more than 3 samples per study where all outcome categories are represented.

The concatenated data matrix \code{X} can contain missing values. Missing values are handled by casewise deletion in the \code{mint.splsda} function without having to delete the rows with missing data.


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
Rohart F. et al (2016, submitted). MINT: A multivariate integrative approach to identify a reproducible biomarker signature across multiple experiments and platforms.

Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group PLS Regression: Application to Epidemiology.
In New Perspectives in Partial Least Squares and Related Methods, pages 243–255. Springer.
}

\author{Florian Rohart}

\seealso{\code{\link{splsda}},
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
keepX.constraint = keepX.constraint, keepX = c(10, 15))

#plot study-specific outputs for study "2", "3" and "4"
plotIndiv(res, study = c(2, 3, 4)) # bad results because the genes on comp 1 and 2 were randomly selected

}

\keyword{regression}
\keyword{multivariate}
