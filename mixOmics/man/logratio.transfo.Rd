\name{logratio.transfo}
\encoding{latin1}
\alias{logratio.transfo}

\title{Log-ratio transformation}

\description{
This function applies a log transformation to the data, either CLR or ILR
}

\usage{
logratio.transfo(X, logratio = "none")
}	

\arguments{
\item{X}{numeric matrix of predictors}
\item{logratio}{log-ratio transform to applied, one of "none", "CLR" or "ILR"}

}

\details{

\code{logratio.transfo} applies a log transformation to the data, either CLR or ILR.
}

\value{
\code{explained_variance} simply returns the log-ratio transform data.
}


\author{Florian Rohart}

\seealso{\code{\link{pca}}, \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}}, \code{\link{splsda}}.}

