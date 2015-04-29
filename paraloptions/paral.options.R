### paral.options.R:  Package specific options mechanism.
###
### $Id: paral.options.R
###
### Implements a slightly modified version of the sm.options() as found in
### sm 2.1-0.  The difference is that the option list is stored in an
### environment '.mix.data'.

## The list of initial options:
.mix.data <- new.env(parent = emptyenv())
.mix.data$options <- list(parallel = NULL)

paral.options <- function(...) {
    if (nargs() == 0) return(.mix.data$options)
    current <- .mix.data$options
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.mix.data$options[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
    .mix.data$options <- current
    invisible(current)
}

## Internal function to apply FUN over X, optionally in parallel:
applyFunc <- function(parSpec, X,Y, FUN, nonForkInit) {
  if (is.null(parSpec) || (is.numeric(parSpec) && parSpec == 1)) {
    ## Serially
    results <- apply(X,Y, FUN)
  } else {
    ## Parallel
    require(parallel, warn.conflicts = FALSE)
    stop_cluster <- FALSE           # Whether to kill the workers afterwards
    
    if (is.call(parSpec)) {
      ## Unevaluated call => evaluate it to create the cluster:
      parSpec <- eval(parSpec)
      stop_cluster <- TRUE
    }
      
    if (inherits(parSpec, "cluster")) {
      ## Run library(mixOmics) on cluster if type != FORK
      if (!inherits(parSpec[[1]], "forknode") && !missing(nonForkInit)) { eval(nonForkInit) }
      results <- parApply(parSpec, X,Y, FUN)
        
      if (stop_cluster) {
        stopCluster(parSpec)
      }
    } else {
      stop("Unknown parallelity specification: '", parSpec, "'")
    }
  }
  return(results)
}

## Internal function to lapply FUN over X, optionally in parallel:
lapplyFunc <- function(parSpec, X, FUN, nonForkInit) {
  if (is.null(parSpec) || (is.numeric(parSpec) && parSpec == 1)) {
    ## Serially
    results <- lapply(X, FUN)
  } else {
    ## Parallel
    require(parallel, warn.conflicts = FALSE)
    stop_cluster <- FALSE           # Whether to kill the workers afterwards
    
    if (is.call(parSpec)) {
      ## Unevaluated call => evaluate it to create the cluster:
      parSpec <- eval(parSpec)
      stop_cluster <- TRUE
    }
      
    if (inherits(parSpec, "cluster")) {
      ## Run library(mixOmics) on cluster if type != FORK
      if (!inherits(parSpec[[1]], "forknode") && !missing(nonForkInit)) { eval(nonForkInit) }
      results <- parLapply(parSpec, X, FUN)
        
      if (stop_cluster) {
        stopCluster(parSpec)
      }
    } else {
      stop("Unknown parallelity specification: '", parSpec, "'")
    }
  }
  return(results)
}

## Internal function to sapply FUN over X, optionally in parallel:
sapplyFunc <- function(parSpec, X, FUN, nonForkInit) {
  if (is.null(parSpec) || (is.numeric(parSpec) && parSpec == 1)) {
    ## Serially
    results <- sapply(X, FUN)
  } else {
    ## Parallel
    require(parallel, warn.conflicts = FALSE)
    stop_cluster <- FALSE           # Whether to kill the workers afterwards
    
    if (is.call(parSpec)) {
      ## Unevaluated call => evaluate it to create the cluster:
      parSpec <- eval(parSpec)
      stop_cluster <- TRUE
    }
      
    if (inherits(parSpec, "cluster")) {
      ## Run library(mixOmics) on cluster if type != FORK
      if (!inherits(parSpec[[1]], "forknode") && !missing(nonForkInit)) { eval(nonForkInit) }
      results <- parSapply(parSpec, X, FUN)
        
      if (stop_cluster) {
        stopCluster(parSpec)
      }
    } else {
      stop("Unknown parallelity specification: '", parSpec, "'")
    }
  }
  return(results)
}
