# -----------------------------------------------------------------------------------
# testing_imgCor.R
# Author:    KA Le Cao 
# Date started:  24/07/2015
# Last updated:  
# Objective: testing imgCor now that CIM has been modified
# Latest update: 
# -----------------------------------------------------------------------------------


rm(list=ls())

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# Florian
sourceDir("/Users/florian/Work/git/package-mixOmics/mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

# KA
sourceDir("../../mixOmics/R/",trace=FALSE) #load all the functions inside ixOmics/R

# -------------------------------------
# from help file
# ------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

## 'combine' type plot (default)
imgCor(X, Y)

## 'separate' type plot
## Not run: 
imgCor(X, Y, type = "separate")

## 'separate' type plot without the name of datas
imgCor(X, Y, X.names = FALSE, Y.names = FALSE, type = "separate")
