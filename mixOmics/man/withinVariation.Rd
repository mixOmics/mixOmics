\name{withinVariation}
\encoding{latin1}


\title{Within matrix decomposition for repeated measurements (cross-over design)}

\description{
This function is called by the multilevel function for cross-over design data, but can be called independently for other types of analyses (e.g. PCA multilevel).
}

\usage{
withinVariation <- function(X, design)
}	

\arguments{
  \item{X}{numeric matrix of predictors. \code{NA}s are allowed.}
    
  \item{design}{a numeric matrix or data frame of 2 columns for a one-factor discrete
    outcome, or of 3 columns for two-factor discrete outcome. The first column 
    indicates the repeated measures on each individual, i.e. the individuals ID. See Details.}
}

\details{

\code{withinVariation} function decomposes the Within variation in the \eqn{X} data set. The resulting \eqn{Xw} matrix is then input in the 
\code{multilevel} function. 

One or two-factor analyses are available.
}

\value{
\code{withinVariation} simply returns the \eqn{Xw} within matrix, which can be input in the other multivariate approaches already implemented in mixOmics (i.e. spls or splsda, see \code{multilevel}, but also pca or ipca).

}

\references{
On multilevel analysis:

Liquet, B., L\^e Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel approach for
biomarker selection and the integration of repeated measures experiments from two
platforms. \emph{BMC Bioinformatics} \bold{13}:325.

Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K. (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA. \emph{Metabolomics},
\bold{6}(1), 119-128.

}

\author{Benoit Liquet, Kim-Anh L\^e Cao, Benoit Gautier, Ignacio Gonz\`alez.}

\seealso{\code{\link{multilevel}}, \code{\link{spls}}, \code{\link{splsda}}, 
\code{\link{plotIndiv}}, \code{\link{plotVar}}, 
\code{\link{plot3dIndiv}}, \code{\link{plot3dVar}},
\code{\link{cim}}, \code{\link{network}}.}

\examples{
## Example: one-factor analysis matrix decomposition
#--------------------------------------------------------------
data(vac18)
X <- vac18$genes
# sample indicates the repeated measurements
design <- data.frame(sample = vac18$sample, 
                     stimul = vac18$stimulation)
Xw <- withinVariation(X = X, design = design)
# multilevel PCA
res.pca.1level <- pca(Xw, ncomp = 3)

# compare a normal PCA with a multilevel PCA for repeated measurements.
# note: PCA makes the assumptions that all samples are independent, so this analysis is flawed!
res.pca <- pca(X, ncomp = 3)

# set up colors for plotIndiv
col.stim <- c("darkblue", "purple", "green4","red3")
col.stim <- col.stim[as.numeric(Y)]

# plotIndiv comparing both PCA and PCA multilevel
plotIndiv(res.pca, ind.names = vac18$stimulation, col = col.stim)
title(main = 'PCA ')
plotIndiv(res.pca.1level, ind.names = vac18$stimulation, col = col.stim)
title(main = 'PCA multilevel')
}

\keyword{regression}
\keyword{multivariate}
