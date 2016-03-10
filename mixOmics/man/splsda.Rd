\name{splsda}
\encoding{latin1}
\alias{splsda}

\title{Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)}

\description{Function to perform sparse Partial Least Squares to classify samples (supervised analysis) and select variables.
}

\usage{
splsda(X,
Y,
ncomp = 2,
mode = c("regression", "canonical", "invariant", "classic"),
keepX,
keepX.constraint=NULL,
scale = TRUE,
tol = 1e-06,
max.iter = 500,
near.zero.var = FALSE,
logratio="none", 
multilevel=NULL)
}	

\arguments{
\item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
\item{Y}{numeric vector or matrix of responses (for multi-response models).
\code{NA}s are allowed.}
\item{ncomp}{the number of components to include in the model (see Details).
Default is set to from one to the rank of \code{X}.}
\item{mode}{character string. What type of algorithm to use, (partially) matching
one of \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}.
See Details.}
\item{keepX.constraint}{A list containing which variables of X are to be kept on each of the first PLS-components}
\item{keepX}{numeric vector of length \code{ncomp}, the number of variables
to keep in \eqn{X}-loadings. By default all variables are kept in the model.}
\item{scale}{boleean. If scale = TRUE, each block is standardized
to zero means and unit variances (default: TRUE)}
\item{tol}{Convergence stopping value.}
\item{max.iter}{integer, the maximum number of iterations.}
\item{near.zero.var}{boolean, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE}
\item{logratio}{one of ('none','CLR'). Default to 'none'}
\item{multilevel}{Design matrix for multilevel anaylis (for repeated measurements). A numeric matrix or data frame. The first column
indicates the repeated measures on each individual, i.e. the individuals ID. If \code{  method = 'splsda'}, the 2nd and 3rd columns are factors. If \code{method = 'spls'} then you can choose to only input the repeated measures (column 1) or the 2nd AND 3rd columns to split the variation for a 2 level factor. See Details.}
}
\details{
\code{splsda} function fit sPLS models with \eqn{1, \ldots ,}\code{ncomp} components
to the factor or class vector \code{Y}. The appropriate indicator (dummy)
matrix is created.
}

\value{
\code{splsda} returns an object of class \code{"splsda"}, a list 
that contains the following components:

  \item{X}{the centered and standardized original predictor matrix.}
  \item{Y}{the centered and standardized indicator response vector or matrix.}
  \item{ind.mat}{the indicator matrix.}
  \item{ncomp}{the number of components included in the model.}
  \item{keepX}{number of \eqn{X} variables kept in the model on each component.}
  \item{mat.c}{matrix of coefficients to be used internally by \code{predict}.}
  \item{variates}{list containing the variates.}
  \item{loadings}{list containing the estimated loadings for the \code{X} and 
	\code{Y} variates.}
  \item{names}{list containing the names to be used for individuals and variables.}
  \item{nzv}{list containing the zero- or near-zero predictors information.}
   \item{tol}{the tolerance used in the iterative algorithm, used for subsequent S3 methods}
  \item{max.iter}{the maximum number of iterations, used for subsequent S3 methods}
\item{iter}{Number of iterations of the algorthm for each component}
}

\references{
On sPLS-DA:
Le Cao, K.-A., Boitard, S. and Besse, P. (2011). Sparse PLS Discriminant Analysis: biologically relevant feature selection and graphical displays for multiclass problems. \emph{BMC Bioinformatics} \bold{12}:253.
}

\author{Ignacio Gonzalez, Kim-Anh Le Cao and Pierre Monget.}

\seealso{\code{\link{spls}}, \code{\link{summary}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{cim}}, \code{\link{network}}, \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for more details.}

\examples{
## First example
data(breast.tumors)
X <- breast.tumors$gene.exp
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(breast.tumors$sample$treatment)

res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))


# individual names appear
plotIndiv(res, ind.names = Y, add.legend = TRUE, plot.ellipse =TRUE)

## Second example: one-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample,
stimul = vac18$stimulation)

# multilevel sPLS-DA model
res.1level <- splsda(X, ncomp = 3, multilevel = design,
    keepX = c(30, 137, 123))

# set up colors for plotIndiv
col.stim <- c("darkblue", "purple", "green4","red3")
plotIndiv(res.1level, ind.names = Y, col.per.group = col.stim)

## Third example: two-factor analysis with sPLS-DA, selecting a subset of variables
# as in the paper Liquet et al.
#--------------------------------------------------------------
\dontrun{
data(vac18.simulated) # simulated data

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
stimu = vac18.simulated$stimulation,
time = vac18.simulated$time)

res.2level <- splsda(X, ncomp = 2, multilevel = design,
keepX = c(200, 200))

plotIndiv(res.2level, group = design$stimu, ind.names = vac18.simulated$time,
add.legend = TRUE, style = 'lattice')
}



## Fourth example
\dontrun{
## Second example
data(liver.toxicity)
X <- as.matrix(liver.toxicity$gene)
# Y will be transformed as a factor in the function,
# but we set it as a factor to set up the colors.
Y <- as.factor(liver.toxicity$treatment[, 4])

splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))

# individual name is set to the treatment
plotIndiv(splsda.liver, ind.names = Y, plot.ellipse = TRUE, add.legend = TRUE)
}
}

\keyword{regression}
\keyword{multivariate}
