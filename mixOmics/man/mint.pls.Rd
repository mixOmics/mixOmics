\name{mint.pls}
\encoding{latin1}
\alias{mint.pls}

\title{Vertical sparse Partial Least Squares (PLS) Regression}

\description{Integration of multilple data sets (or studies) matching on the same variables, ie. vertical Partial Least Squares (PLS) regression.}


\usage{
mint.pls(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
study,
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
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances. Default = \code{TRUE}.}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default = \code{FALSE}.}

}


\details{
\code{mint.splsda} function fits a vertical sPLS-DA models with \code{ncomp} components in which several studies matching on the same variables are integrated together in the aim of classifying a discrete outcome while selecting features. The grouping factor that indicates which samples are from which study is specified with \code{study}. We advise more than 3 samples per study where all outcome categories are represented.

\code{mint.pls} function fits a vertical PLS models with \code{ncomp} components. in which several studies matching on the same variables are integrated together in the aim of classifying a discrete outcome while selecting features. The grouping factor that indicates which samples are from which study is specified with \code{study}.

Multi-response models are fully supported. \code{X} and \code{Y} can contain missing values. Missing values are handled by casewise deletion in the \code{mint.pls} function without having to delete the rows with missing data.


The type of algorithm to use is specified with the \code{mode} argument. Four PLS
algorithms are available: PLS regression \code{("regression")}, PLS canonical analysis
\code{("canonical")}, redundancy analysis \code{("invariant")} and the classical PLS
algorithm \code{("classic")} (see References).

Useful graphical outputs are available, e.g. \code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}.
}

\value{
\code{mint.pls} returns an object of class \code{"mint.pls", "pls"}, a list
that contains the following components:

\item{X}{the centered and standardized original predictor matrix.}
\item{Y}{the centered and standardized original response vector or matrix.}
\item{ncomp}{the number of components included in the model.}
\item{study}{The study grouping factor}
\item{mode}{the algorithm used to fit the model.}
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
In New Perspectives in Partial Least Squares and Related Methods, pages 243-255. Springer.
}

\author{Florian Rohart}

\seealso{\code{\link{pls}},
\code{\link{plotIndiv}}, \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{tune.mint.splsda}} and http://www.mixOmics.org for more details.}

\examples{
data(stemcells)
# need to unmap Y for an unsupervised analysis
Y = unmap(stemcells$celltype)
res = mint.pls(X = stemcells$gene, Y = Y, ncomp = 3, study = stemcells$study)

plotIndiv(res)
#plot study-specific outputs for all studies
plotIndiv(res, study = "all.partial", col = as.numeric(stemcells$celltype))

#plot study-specific outputs for study "2"
plotIndiv(res, study = "2", col = as.numeric(stemcells$celltype))

}

\keyword{regression}
\keyword{multivariate}
