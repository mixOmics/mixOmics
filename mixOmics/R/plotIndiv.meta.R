# Copyright (C) 2014
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
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

# last modified: 01-03-2016


#--------------------------------------------------------------------#
#-- Includes plotIndiv for mint.pls, mint.spls, mint.plsda, mint.splsda --#
#--------------------------------------------------------------------#


#--------------------- PLS and sPLS ---------------------#
plotIndiv.mint.spls <- plotIndiv.mint.splsda <-
function(object,
comp = 1:2,
ind.names = TRUE,
study="all",
rep.space = "X-variate",
X.label = NULL,
Y.label = NULL,
col,
pch,
cex = 1,
abline.line = FALSE,
layout = NULL,
...)
{
    
    
    
    
    
    # validation des arguments #
    #--------------------------#
    if (length(comp) != 2)
    stop("'comp' must be a numeric vector of length 2.")
    
    if (!is.numeric(comp))
    stop("invalid vector for 'comp'.")
    
    if (length(comp) == 1)
    stop("Need at least 2 components to plot the graph")
    
    if (any(comp > object$ncomp))
    stop("the elements of 'comp' must be smaller or equal than ", object$ncomp, ".")
    
    if (is.logical(ind.names)) {
        if (isTRUE(ind.names)) ind.names = object$names$indiv
    }
    
    if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$X))
        stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
    }
    
    if(any(!study%in%levels(object$study)) & any(study!="all"))
    {
        stop("'study' must be one of object$study or 'all'.")
    }
    
    
    comp1 = round(comp[1])
    comp2 = round(comp[2])
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate"))
    
    

    
    if(any(study=="all"))
    {
        
        if(any(class(object)%in%c("mint.plsda","mint.splsda")) & missing(col))
        {
            if (nlevels(object$Y) < 10) {
                #only 10 colors in color.mixo
                col = color.mixo(as.numeric(object$Y))
            } else {
                #use color.jet
                col = color.jet(as.numeric(object$Y))
            }
        }else if (missing(col))
        {
            col=1#rep(col,ceiling(object$N/length(col)))[1:object$N] # complete col to have a vector of length N, so we can choose col per study
        }
        if(missing(pch))
        {
            pch = as.numeric(object$study)
        }else{
            pch=rep(pch,ceiling(nrow(object$X)/length(col)))[1:nrow(object$X)] # complete pch to have a vector of length N, so we can choose pch per study
            
        }
        
        # l'espace de representation #
        #----------------------------#
        if (rep.space == "X-variate"){
            x = object$variates$X[, comp1]
            y = object$variates$X[, comp2]
            if (is.null(X.label)) X.label = paste("X-variate", comp1)
            if (is.null(Y.label)) Y.label = paste("X-variate", comp2)
        }
        
        if (rep.space == "Y-variate"){
            x = object$variates$Y[, comp1]
            y = object$variates$Y[, comp2]
            if (is.null(X.label)) X.label = paste("Y-variate", comp1)
            if (is.null(Y.label)) Y.label = paste("Y-variate", comp2)
        }
        
        if (rep.space == "XY-variate"){
            x = (object$variates$X[, comp1] + object$variates$Y[, comp1]) / 2
            y = (object$variates$X[, comp2] + object$variates$Y[, comp2]) / 2
            if (is.null(X.label)) X.label = paste("X-variate", comp1)
            if (is.null(Y.label)) Y.label = paste("Y-variate", comp2)
        }
        
        
        # le plot des individus #
        #-----------------------#
        if (length(ind.names) > 1)
        {
            plot(x, y, type = "n", xlab = X.label, ylab = Y.label, ...)
            text(x, y, ind.names, col = col, cex = cex, ...)
            if(abline.line) abline(v = 0, h = 0, lty = 2)
        }else{
            if (isTRUE(ind.names))
            {
                plot(x, y, type = "n", xlab = X.label, ylab = Y.label, ...)
                text(x, y, ind.names, col = col, cex = cex, ...)
                if(abline.line) abline(v = 0, h = 0, lty = 2)
            }else{
                plot(x, y, xlab = X.label, ylab = Y.label,
                col = col, cex = cex, pch = pch, ...)
                if(abline.line) abline(v = 0, h = 0, lty = 2)
            }
        }
        
    }else{ #loop on the levels
        nResp = length(study)  # Number of study
        #print(nResp)
        if (nResp > 1)
        {
            if (is.null(layout))
            {
                nRows = min(c(3, nResp))
                nCols=min(c(3,ceiling(nResp/nRows)))
                layout = c(nRows, nCols)
            }else{
                if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
                stop("'layout' must be a numeric vector of length 2.")
                nRows = layout[1]
                nCols = layout[2]
            }
            if (nRows * nCols < nResp) devAskNewPage(TRUE)
        }
        opar=par("mfrow") # save the current mfrow to load it at the end of the function
        #print(layout)
        par(mfrow=layout)
        
        
        if(any(class(object)%in%c("mint.plsda","mint.splsda")) & missing(col))
        {
            if (nlevels(object$Y) < 10) {
                #only 10 colors in color.mixo
                col = color.mixo(as.numeric(object$Y))
            } else {
                #use color.jet
                col = color.jet(as.numeric(object$Y))
            }
        }else if (missing(col))
        {
            col=rep(1,nrow(object$X)) # complete col to have a vector of length N, so we can choose col per study
        }
        if(missing(pch))
        {
            pch = as.numeric(object$study)
        }else{
            pch=rep(pch,ceiling(nrow(object$X)/length(col)))[1:nrow(object$X)] # complete pch to have a vector of length N, so we can choose pch per study
            
        }

        for (m in 1:length(study))
        {
            

            ind=match(study[m],levels(object$study))
            #print(ind)
            ind.X=which(object$study==study[m])
            #print(ind.X)
            col.temp=col[ind.X]
            pch.temp =pch[ind.X]

    
            # l'espace de representation #
            #----------------------------#
            if (rep.space == "X-variate"){
                x = object$variates.partial$X[[ind]][, comp1]
                y = object$variates.partial$X[[ind]][, comp2]
                if (is.null(X.label)) X.label = paste("X-variate", comp1)
                if (is.null(Y.label)) Y.label = paste("X-variate", comp2)
            }
            #print(dim(x))
            
            if (rep.space == "Y-variate"){
                x = object$variates.partial$Y[[ind]][, comp1]
                y = object$variates.partial$Y[[ind]][, comp2]
                if (is.null(X.label)) X.label = paste("Y-variate", comp1)
                if (is.null(Y.label)) Y.label = paste("Y-variate", comp2)
            }
            
            if (rep.space == "XY-variate"){
                x = (object$variates.partial$X[[ind]][, comp1] + object$variates.partial$Y[[ind]][, comp1]) / 2
                y = (object$variates.partial$X[[ind]][, comp2] + object$variates.partial$Y[[ind]][, comp2]) / 2
                if (is.null(X.label)) X.label = paste("X-variate", comp1)
                if (is.null(Y.label)) Y.label = paste("Y-variate", comp2)
            }
            
            
            
            # le plot des individus #
            #-----------------------#
            if (length(ind.names) > 1)
            {
                plot(x, y, type = "n", xlab = X.label, ylab = Y.label, ...)
                text(x, y, ind.names, col = col.temp, cex = cex, pch = pch.temp, ...)
                if(abline.line) abline(v = 0, h = 0, lty = 2)
            }else{
                if (isTRUE(ind.names))
                {
                    plot(x, y, type = "n", xlab = X.label, ylab = Y.label, ...)
                    text(x, y, ind.names, col = col.temp, cex = cex, pch = pch.temp, ...)
                    if(abline.line) abline(v = 0, h = 0, lty = 2)
                }else{
                    plot(x, y, xlab = X.label, ylab = Y.label,
                    col = col.temp, cex = cex, pch = pch.temp, ...)
                    if(abline.line) abline(v = 0, h = 0, lty = 2)
                }
            }
            
         title(study[m],...)
        }
        
        #par(mfrow=opar) # commented out so that the user can keep plotting graphs next to these ones
        if (nRows * nCols < nResp) devAskNewPage(FALSE)

    }

}




