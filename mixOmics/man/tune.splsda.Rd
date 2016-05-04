\name{tune.splsda}
\encoding{latin1}
\alias{tune.splsda}


\title{Tuning functions for multilevel analyses}

\description{
These functions were implemented to help tuning the variable selection parameters in the multilevel analyses.
}

\usage{
tune.splsda(X, Y, ncomp = 1,
test.keepX = c(5, 10, 15), already.tested.X = NULL, validation = "Mfold", folds = 10,
dist = "max.dist", measure = "BER", progressBar = TRUE, near.zero.var = FALSE,
nrepeat = 1, logratio = c('none','CLR'), multilevel = NULL, light.output = TRUE)

}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
  \item{Y}{\code{if(method = 'spls')} numeric vector or matrix of continuous responses (for multi-response models) \code{NA}s are allowed.}
  \item{ncomp}{the number of components to include in the model.}
  \item{test.keepX}{numeric vector for the different number of variables to test from the \eqn{X} data set}
  \item{already.tested.X}{if \code{ncomp > 1} numeric vector indicating the number of variables to select rom the \eqn{X} data set on the previous ncomp-1 components}
\item{validation}{character.  What kind of (internal) validation to use, matching one of \code{"Mfold"} or
\code{"loo"} (see below). Default is \code{"Mfold"}.}
\item{folds}{the folds in the Mfold cross-validation. See Details.}
\item{dist}{distance metric to use for \code{splsda} to estimate the classification error rate,
should be a subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).}
\item{measure}{Two misclassification measure are available: overall misclassification error \code{overall} or the Balanced Error Rate \code{BER}}
\item{progressBar}{by default set to \code{TRUE} to output the progress bar of the computation.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Default value is FALSE}
\item{nrepeat}{Number of times the Cross-Validation process is repeated.}
\item{logratio}{one of ('none','CLR'). Default to 'none'}
\item{multilevel}{Design matrix for multilevel anaylis (for repeated measurements) that indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}

\item{light.output}{if set to FALSE, the prediction/classification of each sample for each of \code{test.keepX} and each comp is returned.}

}

\details{
This tuning function should be used to tune the parameters in the \code{multilevel} function.

For a sPLS-DA one-factor analysis, M-fold cross-validation is performed, internally the training data is decomposed into within-subject variation.

For a sPLS-DA two-factor analysis, the correlation between components from the within-subject variation of X and the \code{cond} matrix is computed on the whole data set. The reason why we cannot obtain a corss-validation error rate as for the spls-DA one-factor analysis is because of the dififculty to decompose and predict the within matrices within each fold.

For a sPLS two-factor analysis a sPLS canonical mode is run, and the correlation between components from the within-subject variation of X and Y is computed on the whole data set.

If \code{validation = "Mfold"}, M-fold cross-validation is performed. 
How many folds to generate is selected by specifying the number of folds in \code{folds}.
The folds also can be supplied as a list of vectors containing the indexes defining each 
fold as produced by \code{split}.

If \code{validation = "loo"}, leave-one-out cross-validation is performed. By default \code{folds} is set to the number of unique individuals.


}

\value{
Depending on the type of analysis performed, a list that contains:
  \item{error}{cross-validation overall error rate when one-factor sPLS-DA analysis is performed.}
  \item{prediction.all}{cross-validation prediction for all samples and for the LAST keepX tested parameter when one-factor sPLS-DA analysis is performed.}
  \item{cor.value}{compute the correlation between latent variables for two-factor sPLS-DA analysis or sPLS.}

}


\author{Kim-Anh Le Cao, Benoit Gautier, Francois Bartolo, Florian Rohart.}

\seealso{\code{\link{splsda}} and http://www.mixOmics.org for more details.}

\examples{
## First example: analysis with sPLS-DA
\dontrun{
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- as.factor(breast.tumors$sample$treatment)
tune= tune.splsda(X,Y,ncomp=1,nrepeat=10,logratio="none",
    test.keepX = c(5, 10, 15),folds=10,dist="max.dist",
    already.tested.X=NULL, progressBar = TRUE)


}


## Second example: multilevel one-factor analysis with sPLS-DA
\dontrun{
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample)

tune= tune.splsda(X, Y = Y, ncomp=3, nrepeat=10, logratio="none",
    test.keepX = c(5,50,100),folds=10, dist="max.dist", multilevel = design)

}


}

\keyword{regression}
\keyword{multivariate}
