\name{tune.spls}
\encoding{UTF-8}
\alias{tune.spls}


\title{Tuning functions for sPLS method}

\description{
Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
grid to determine optimal values for the sparsity parameters in \code{spls}.
}

\usage{
tune.spls(X, Y, ncomp = 1,
test.keepX = c(5, 10, 15), test.keepY = ncol(Y), already.tested.X, already.tested.Y,
validation = "Mfold", folds = 10, measure = "MSE", scale = TRUE,
progressBar = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE,
nrepeat = 1, multilevel = NULL, light.output = TRUE, cpus)

}	

\arguments{
    \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
    \item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
    \item{ncomp}{the number of components to include in the model.}
    \item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
    \item{test.keepY}{numeric vector for the different number of variables to test from the \eqn{Y} data set}
    \item{already.tested.X}{Optional, if \code{ncomp > 1} A numeric vector indicating the number of variables to select from the \eqn{X} data set on the firsts components.}
    \item{already.tested.Y}{Optional, if \code{ncomp > 1} A numeric vector indicating the number of variables to select from the \eqn{Y} data set on the firsts components.}
    \item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or
\code{"loo"} (see below). Default is \code{"Mfold"}.}
    \item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{measure}{One of \code{MSE, MAE, Bias} or \code{R2}. Default to \code{MSE}. See details}
    \item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
    \item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
    \item{tol}{Convergence stopping value.}
    \item{max.iter}{integer, the maximum number of iterations.}
    \item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
    \item{nrepeat}{Number of times the Cross-Validation process is repeated.}
    \item{multilevel}{Design matrix for multilevel analysis (for repeated measurements) that indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}

    \item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}
    \item{cpus}{Number of cpus to use when running the code in parallel.}

}

\details{

This tuning function should be used to tune the parameters in the \code{spls} function (number of components and both the number of variables in \code{keepX} and \code{keepY} to select).

If \code{validation = "loo"}, leave-one-out cross-validation is performed. By default \code{folds} is set to the number of unique individuals.
If \code{validation = "Mfold"}, M-fold cross-validation is performed.
How many folds to generate is selected by specifying the number of folds in \code{folds}.

Four measure of accuracy are available: Mean Square Error(\code{MSE}), Mean Absolute Error (\code{MAE}), \code{Bias} and \code{R2}.


The function outputs the optimal number of components that achieve the best performance based on the chosen measure of accuracy. The assessment is data-driven and similar to the process detailed in (Rohart et al., 2016), where one-sided t-tests assess whether there is a gain in performance when adding a component to the model.

}

\value{
Depending on the type of analysis performed, a list that contains:
\item{error.rate}{returns the prediction error for each \code{test.keepX} on each component, averaged across all repeats and subsampling folds. Standard deviation is also output. All error rates are also available as a list.}
\item{choice.keepX}{returns the number of variables selected (optimal keepX) on each component.}
\item{choice.keepY}{returns the number of variables selected (optimal keepY) on each component.}
\item{choice.ncomp}{returns the optimal number of components for the model fitted with \code{$choice.keepX} and \code{$choice.keepY} }

\item{predict}{Prediction values for each sample, each \code{test.keepX,test.keepY}, each comp and each repeat. Only if light.output=FALSE}

\item{cor.value}{only if multilevel analysis with 2 factors: correlation between latent variables.}

}

\references{
mixOmics article:

Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics feature selection and multiple data integration. PLoS Comput Biol 13(11): e1005752

}



\author{Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart.}

\seealso{\code{\link{splsda}}, \code{\link{predict.splsda}} and http://www.mixOmics.org for more details.}

\examples{
\dontrun{
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

tune = tune.spls(X, Y, ncomp=4, test.keepX = c(5,10,15),
nrepeat=3, test.keepY=c(2,5,8,10), progressBar = TRUE)

tune$choice.ncomp
tune$choice.keepX
tune$choice.keepY

# plot only for a specific keepY
plot(tune, keepY = 10)
plot(tune, keepY = 5)
}







}

\keyword{regression}
\keyword{multivariate}
