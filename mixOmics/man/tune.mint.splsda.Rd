\name{tune.mint.splsda}
\encoding{latin1}
\alias{tune.mint.splsda}

\title{Estimate the parameters of regularization for Regularized CCA}

\description{
Computes leave-one-out or M-fold cross-validation scores on a two-dimensional 
grid to determine optimal values for the parameters of regularization in 
\code{rcc}.
}

\usage{
tune.mint.splsda(X, Y, ncomp = 1, study, test.keepX = c(5, 10, 15), already.tested.X,
dist = "max.dist", measure = "BER", progressBar = TRUE, scale = TRUE, tol = 1e-06,
max.iter = 500, near.zero.var = FALSE, light.output = TRUE )
}

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{Outcome. Numeric vector or matrix of responses (for multi-response models)}
\item{ncomp}{Number of components to include in the model (see Details). Default to 1}
\item{study}{grouping factor indicating which samples are from the same study}
\item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
\item{already.tested.X}{if \code{ncomp > 1} A list containing which variables of X are to be kept on each of the first PLS-components}
\item{dist}{only applies to an object inheriting from \code{"plsda"} or \code{"splsda"} to evaluate the classification performance of the model. Should be a subset of \code{"max.dist"}, \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}. See \code{\link{predict}}.}
\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}
}


\details{
This function performs a Leave-One-Group-Out-Cross-Validation (LOGOCV), where each of \code{study} is left out once. It returns a list of variables of \code{X} that were selected on each of the \code{ncomp} components.
Then, a \code{\link{mint.splsda}} can be performed with \code{keepX.constraint} set as the output \code{choice.keepX}.

All component \eqn{1:\code{ncomp}} are tuned, except the first one for which a \code{already.test.X} is provided. See examples below.


}

\value{
The returned value is a list with components: 
\item{mat.mean.error,}{returns the prediction error for each \code{test.keepX} on each component}
\item{choice.keepX}{returns the variables selected with the optimal keepX on each component.}
\item{error.per.class}{returns the classification accuracy for each level of \code{Y} and for each component computed with the optimal keepX}

  \item{predict}{Prediction values for each sample, each \code{test.keepX} and each comp.}
  \item{class}{Predicted class for each sample, each \code{test.keepX} and each comp.}
}

\author{Florian Rohart}

\seealso{\code{\link{mint.splsda}} and http://www.mixOmics.org for more details.}

\examples{
data = stemcells$gene
type.id = stemcells$celltype
exp = stemcells$study

res = mint.splsda(X=data,Y=type.id,ncomp=3,keepX=c(10,5,15),study=exp)

\dontrun{
out = tune.mint.splsda(X=data,Y=type.id,ncomp=2,near.zero.var=FALSE,
    study=exp,test.keepX=seq(1,10,1))
out$choice.keepX
}

}

\keyword{multivariate}
\keyword{dplot}
