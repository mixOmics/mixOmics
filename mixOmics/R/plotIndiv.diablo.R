################################################
#
# visualizationFunctions.R
# Author: Amrit Singh
# Date: April 01, 2016
#
# functions
#     1) plotIndiv_diablo; 2 sub-functions; splotMatPlot() and panel.ellipses
#     2) circosPlot_diablo
#     3) heatmap_diablo
#     4) enrichPathwayNetwork_diablo
#
################################################

## plotIndiv_diablo
plotIndiv.sgccda = function(object, Y, ncomp = 1, groupOrder = levels(Y)){
    #library(mixOmics)  ## needed for color.mixo
  VarX <- do.call(cbind, lapply(object$variates, function(i) i[, ncomp]))
  datNames <- colnames(VarX)
  
  if (!is.factor(Y))
    stop(gettextf("Y must be a factor!"))
  if (length(ncomp) != 1)
    stop(gettextf("You can only choose one component"))
  
  numberOfCols <- ncol(VarX)
  numberOfRows <- numberOfCols - 1
  
  mat <- matrix(0, nrow = numberOfRows, ncol = numberOfRows)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      mat[i,j] <- paste(i,j, sep="_")
    }
  }
  plotType = list(cor=mat[lower.tri(mat)], scatter=mat[upper.tri(mat)],
                  lab=diag(mat), 
                  bar=paste(1:(numberOfRows-1), numberOfCols, sep="_"),
                  stackedbar=paste(numberOfRows, numberOfCols, sep="_"))
  
  par(mfrow = c(numberOfRows, numberOfCols), mar = rep.int(1/2, 4), oma = c(2,2,2,2))
  for(i in 1:numberOfRows){
    for(j in 1:numberOfCols){
      ptype <- unlist(lapply(plotType, function(x){
        intersect(paste(i,j,sep="_"), x)
      }))
      splotMatPlot(x=VarX[, i], y=VarX[, j], datNames, Y, ptype, groupOrder)
      if(i == 1 & j %in% seq(2, numberOfRows, 1)){Axis(side = 3, x=VarX[, i])}
      if(j == numberOfRows & i %in% seq(1, numberOfRows-1, 1)){Axis(side = 4, x=VarX[, i])}
    }
  }
}


splotMatPlot = function(x, y, datNames, Y, ptype, groupOrder){
    if(names(ptype) == "cor"){
        plot(1, type = "n", axes = FALSE)
        r = round(cor(x, y), 2)
        text(1, 1, labels=r, cex = 0.6/strwidth(r)*r)
        box()
    }
    if(names(ptype) == "scatter"){
        panel.ellipses(x=x, y=y, Y = Y)
    }
    if(names(ptype) == "lab"){
        plot(1, type = "n", axes = FALSE)
        ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), unique)))
        text(x=1, y=1, labels=datNames[ind], cex = 2)
        box()
    }
    if(names(ptype) == "bar"){
        Y2 <- factor(as.character(Y), levels = groupOrder)
        boxplot(x ~ Y2, horizontal=TRUE, axes = FALSE, ylim = c(min(x)-3, max(x)),
        col= color.mixo(match(levels(Y2), levels(Y))))
        axis(4, at=1:nlevels(Y2), labels=levels(Y2))
    }
    if(names(ptype) == "stackedbar"){
        Y2 <- factor(as.character(Y), levels = groupOrder)
        bars <- table(Y2)
        barplot(bars, col= color.mixo(match(levels(Y2), levels(Y))),
        axes = FALSE)
        axis(4, at=seq(0,max(bars),length.out=5), labels=seq(0,max(bars),length.out=5))
    }
}

panel.ellipses = function(x, y, Y = Y, pch = par("pch"), col.lm = "red", axes = FALSE,
...) {
    ind.gp = matrice = cdg = variance = list()
    for(i in 1:nlevels(Y)){
        ind.gp[[i]] = which(as.numeric(Y)==i)
    }
    
    matrice = lapply(ind.gp, function(z){matrix(c(x[z], y[z]), ncol = 2)})
    cdg = lapply(matrice, colMeans)
    variance = lapply(matrice, var)
    
    library(ellipse)
    coord.ellipse = lapply(1:nlevels(Y), function(x){ellipse(variance[[x]], centre = cdg[[x]], ellipse.level = 0.95)})
    max.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, max)})
    min.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, min)})
    ind.names <- names(Y)
    cex = 0.5
    
    plot(x, y, xlab = "X.label", ylab = "Y.label", col=color.mixo(as.numeric(Y)), pch=20, axes=axes,
    xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
    #text(x, y, ind.names, col = col, cex = cex)
    box()
    for (z in 1:nlevels(Y)){
        points(coord.ellipse[[z]], type = "l", col = color.mixo(c(1:nlevels(Y))[z]))
    }
}

