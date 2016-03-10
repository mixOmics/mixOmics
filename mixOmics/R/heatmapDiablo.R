#############################################################################################################
# Author :
#   Amrit Singh
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: ?
# last modified: 08-03-2016
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################



################################################
#
## 2) heatmap_diablo
#
################################################

# This function is a small wrapper of cim. For more customisation, please use cim

heatmapDiablo = function(object,
ncomp=1,
margins = c(2, 15),
pos.legend="topleft",
cex.legend=1.5,
groupOrder = NULL)
{
    
    # check input object
    if (!any(class(object) == "block.splsda"))
    stop("heatmapDiablo is only available for 'block.splsda' objects")

    if(length(object$X)<=1)
    stop("This function is only available when there are more than 3 blocks") # so 2 blocks in X + the outcome Y

    if(ncomp>min(object$ncomp))
    stop("'ncomp' needs to be higher than object$ncomp")

    X <- object$X
    Y <- object$Y
    if(is.null(groupOrder)){
        Y2 <- Y
    } else {
        Y2 <- factor(as.character(Y), levels = groupOrder)
    }
    
    #need to reorder variates and loadings to put 'Y' in last
    indY=object$indY
    object$variates=c(object$variates[-indY],object$variates[indY])
    object$loadings=c(object$loadings[-indY],object$loadings[indY])
    
    #reducing loadings for ncomp
    object$loadings=lapply(object$loadings, function(x){x[,1:ncomp,drop=FALSE]})
    
    keepA = lapply(object$loadings, function(i) apply(abs(i), 1, sum) > 0)
    XDatList <- mapply(function(x, y){
        x[, y]
    }, x=X, y=keepA[-length(keepA)])
    XDat <- do.call(cbind, XDatList)
    XDat[which(XDat > 2)] <- 2
    XDat[which(XDat < -2)] <- -2
    
    dark <- brewer.pal(n = 12, name = 'Paired')[seq(2, 12, by = 2)]
    VarLabels <- factor(rep(names(X), lapply(keepA[-length(keepA)], sum)), levels = names(X)[order(names(X))])
    
    ## Plot heatmap
    cim(XDat, row.names = rep("", nrow(XDat)), col.names = rep("", ncol(XDat)),
    col.sideColors = dark[as.numeric(VarLabels)],
    row.sideColors = color.mixo(match(levels(Y2), levels(Y))[as.numeric(Y)]), margins = margins)
    
    legend("topleft",c("Rows",c(levels(Y)[order(levels(Y))],"",
    "Columns",names(X))),col=c(1,color.mixo(match(levels(Y2), levels(Y)))[order(levels(Y))],1,
    1,dark[1:nlevels(VarLabels)][match(levels(VarLabels), names(X))]),
    pch=c(NA,rep(19,nlevels(Y)),NA,NA,rep(19,nlevels(VarLabels))),bty="n",cex=cex.legend,text.font=c(2,rep(1,nlevels(Y)),NA,2,rep(1,nlevels(VarLabels))))
}

