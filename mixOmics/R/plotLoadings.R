#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 15-04-2016
# last modified: 19-04-2016
#
# Copyright (C) 2016
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


plotLoadings  =
function(object, ...) UseMethod("plotLoadings")


#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA --#
#----------------------------------------------------------------------------------------------------------#


#plotLoadings =

plotLoadings.pls =
#plotLoadings.mlpls =      # because pls too
plotLoadings.spls =
#plotLoadings.mlspls =     # because spls too
plotLoadings.rcc =
plotLoadings.pca =
plotLoadings.sgcca =
plotLoadings.rgcca =


function(object, block, #single value
comp = 1,
show.ties = TRUE,
col = NULL,
ndisplay = NULL,
cex.name = 0.7,
cex.legend = 0.8,
name.var = NULL,
complete.name.var = FALSE,
main = NULL,
subtitle,
layout = NULL,
size.title = rel(2),
size.subtitle = rel(1.5),
...)
{
    
    # -------------------
    # input checks
    # ------------------
    
    if (is.null(object$loadings))
    stop("'plotLoadings' should be used on object for which object$loadings is present.")
    
    # block
    # --
    if (missing(block))
    block = object$names$blocks
    
    if(is.numeric(block))
    {
        if(any(block>length(object$names$blocks)))
        stop("'block' needs to be lower than the number of blocks in the fitted model, which is length(object$names$blocks)")
        
    }else if(is.character(block) & any(is.na(match(block,object$names$blocks)))) {
        stop("Incorrect value for 'block', 'block' should be among the blocks used in your object: ", paste(object$names$blocks,collapse=", "), call. = FALSE)
    }
    
    if (!missing(subtitle))
    {
        if (length(subtitle)!=length(block))
        stop("'subtitle' indicates the subtitle of the plot for each block and it needs to be the same length as 'block'.")
    }



    # layout
    # --
    opar = par(no.readonly = TRUE)
    reset.mfrow = FALSE # if set to TRUE, the algorithm ends up with  par(mfrow=reset.mfrow)
    nResp = length(block)#length(names(object$names$colnames)) #number of blocks
    if (is.null(layout))
    {
        # check if there are enough plots in mfrow
        omfrow = par("mfrow")
        available.plots = prod(omfrow)
        if (available.plots<nResp) # if not enough plots available, we create our new plot
        {
            nRows = min(c(3, ceiling(nResp/2)))
            nCols = min(c(3, ceiling(nResp/nRows)))
            layout = c(nRows, nCols)
            par(mfrow = layout)
            
            if (nRows * nCols < nResp)
            devAskNewPage(TRUE)
            
            reset.mfrow=TRUE # we changed mfrow to suits our needs, so we reset it at the end
        }
    } else {
        if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
        stop("'layout' must be a numeric vector of length 2.")
        
        nRows = layout[1]
        nCols = layout[2]
        par(mfrow = layout)
        
        if (nRows * nCols < nResp)
        devAskNewPage(TRUE)
    }



    # cex
    # --
    if (cex.name <= 0)
    cex.name = 0.7
    
    if (cex.legend <= 0)
    cex.name = 0.8
    
    ncomp = object$ncomp
    
    #col
    #-----


    #title
    #-----
    if (!is.null(main) & !is.character(main))
    warning('main needs to be of type character')
    
    #col
    #-----
    if (!is.null(col) & (length(col) !=  1))
    {
        warning('col must be the of length 1, by default set to default colors')
        col = color.mixo(1)  # by default set to the colors in color.mixo (10 colors)
    }
    if (is.null(col))
    col = color.mixo(1) # by default set to the colors in color.mixo (10 colors)
    
    
    for (i in 1 : length(block))
    {
        ##selectvar
        selected.var = selectVar(object, comp = comp, block = block[i]) # gives name and values of the blocks in 'block'
        name.selected.var = selected.var[[1]]$name
        value.selected.var = selected.var[[1]]$value

        # ndisplay
        # ------
        # if null set by default to all variables from selectVar
        if (is.null(ndisplay))
        {
            ndisplay.temp = length(name.selected.var)
        } else if(ndisplay > length(name.selected.var)) {
            message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var))
            ndisplay.temp = length(name.selected.var)
        } else {
            ndisplay.temp = ndisplay
        }
        
        name.selected.var = name.selected.var[1:ndisplay.temp]
        value.selected.var = value.selected.var[1:ndisplay.temp,]
        
        #comp
        # ----
        if (any(class(object) %in% c("pls","spls", "rcc")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
        {
            ncomp = object$ncomp
            object$X = list(X = object$X, Y = object$Y) # so that the data is in object$X, either it's a pls or block approach
        } else {
            ncomp = object$ncomp[block[i]]
        }
        
        if (any(max(comp) > ncomp))
        stop(paste("Argument 'comp' should be less or equal to ", ncomp))
        
        names.block = as.character(names(selected.var)[1]) #it should be one block and ncomp, so we take the first one
        
        X = object$X[names.block][[1]]
        
        #name.var
        ind.match = match(name.selected.var, colnames(X)) # look at the position of the selected variables in the original data X
        if(!is.null(name.var))
        {
            if(length(name.var)!= ncol(X))
            stop("name.var should be a vector of length ", ncol(X))
            
            colnames.X = as.character(name.var[ind.match]) # get the
        }else{
            colnames.X = as.character(colnames(X))[ind.match]
        }
        X = X[, name.selected.var] #reduce the problem to ndisplay
        
        #completing colnames.X by the original names of the variables when missing
        if (complete.name.var == TRUE)
        {
            ind = which(colnames.X == "")
            if (length(ind) > 0)
            colnames.X[ind] = colnames(X)[ind]
        }
        
        # --------------------
        # end check inputs
        # ---------------------

        contrib = data.frame( importance = value.selected.var) # contribution of the loading
        
        # barplot with contributions
        #layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.1), TRUE)
        if (!is.null(main) & length(block) > 1)
        {
            par(mar = c(4, max(6, max(sapply(colnames.X, nchar))/2), 6, 2))
        } else {
            par(mar = c(4, max(6, max(sapply(colnames.X, nchar))/2), 4, 2))
        }

        mp = barplot(contrib$importance, horiz = T, las = 1, col = col, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(contrib),
        cex.names = cex.name, cex.axis = 0.7, beside = TRUE,border = NA)
        
        if ( (length(block) == 1 & is.null(main)) | (length(block) > 1 & missing(subtitle)))
        {
            title(paste0('Loadings on comp ', comp, "\nBlock '", names.block,"'"), line=1, cex.main = size.subtitle)
        } else if (length(block) == 1) {
            title(paste(main), line=1, cex.main = size.subtitle)
        } else if (length(block) > 1 & !missing(subtitle)) {
            title(paste(subtitle[i]), line=1, cex.main = size.subtitle)
        }
    }
    
    if (length(block) > 1 & !is.null(main))
    title(main, outer=TRUE, line = -2, cex.main = size.title)
    
    if (reset.mfrow)
    par(opar)#par(mfrow = omfrow)
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  =
    # return the contribution matrix
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  =
    return(invisible(list(contrib = contrib)))
}
