%% $Id: paral.options.Rd 237 2013-08-05 18:01:18Z bhm $
\encoding{latin1}
\name{paral.options}
\alias{paral.options}
\title{Set or return options for the mixOmics package}
\description{
  A function to set options for the \pkg{mix} package, or to return the
  current options.
}
\usage{
paral.options(\dots)
}
\arguments{
  \item{\dots}{any number
    of named arguments (\var{name = value}).}
}
\details{
  If called with no arguments, or with an empty list as the single
  argument, \code{paral.options} returns the current options.


  Otherwise, \code{paral.options} should be called with one named
  arguments \var{name = value}.  For each argument, the option named
  \var{name} will be given the value \var{value}.

  The recognised options are:
  \describe{
    
    \item{parallel}{Specification of how the cross-validation (CV)  should be performed.  If the specification is
      \code{NULL} (default) or \code{1}, the CV is done serially, otherwise
      it is done in parallel using functionality from the
      \code{\link{parallel}} package.

      If it is an integer greater than 1, the CV is done in parallel
      with the specified number of processes, using
      \code{\link{multicore}}.  

      If it is a cluster object created by \code{\link{makeCluster}},
      the CV is done in parallel on that cluster, using
      \code{\link{snow}}.  The user should stop the cluster herself
      when it is no longer needed, using \code{\link{stopCluster}}.

      Finally, if the specification is an unevaluated call to
      \code{\link{makeCluster}}, the call is evaluated, and
      the CV is done in parallel on the resulting cluster, using
      \code{\link{snow}}.  In this case, the cluster will
      be stopped (with \code{\link{stopCluster}}) after the
      CV.  Thus, in the final case, the cluster is created and
      destroyed for each CV, just like when using
      \code{\link{multicore}}.}
  }
}
\value{
  A list with the (possibly changed) options.  If any named argument (or
  list element) was provided, the list is returned invisibly.
}
\author{Francois Bartolo}
\note{
  The function is a slight modification of the function
  \code{\link[sm]{sm.options}} from the package \pkg{sm}.
}
\examples{


library(parallel)


## Parallelised cross-validation, using persistent cluster:
paral.options(parallel = makeCluster(4, type = "FORK")) # not available on Windows
paral.options(parallel = makeCluster(4, type = "PSOCK"))
## The cluster can be used several times
## The cluster should be stopped manually afterwards:
stopCluster(paral.options()$parallel)

## Parallelised cross-validation, using persistent MPI cluster:
## This requires the package Rmpi to be installed
## This creates the cluster:
paral.options(parallel = makeCluster(4, type = "MPI"))
## The cluster can be used several times
## The cluster should be stopped manually afterwards
stopCluster(paral.options()$parallel)
## It is good practice to call mpi.exit() or mpi.quit() afterwards:
mpi.exit()

## Restore `factory settings':
rm(.paral.options)
paral.options()
}

## End(Not run)
