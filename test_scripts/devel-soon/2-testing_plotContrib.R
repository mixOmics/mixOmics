# -----------------------------------------------------------------------------------
# Testing-contrib.R
# Author:    KA Le Cao 
# Date started:  22/07/2015
# Last updated:  
# Objective:
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


# =================================
# first test: PLSDA
# ================================
data(liver.toxicity)
X <- liver.toxicity$clinic
Y <- liver.toxicity$treatment[, 3]

plsda.liver <- plsda(X, Y, ncomp = 2)


contrib.plsda = plotCOotrib()


