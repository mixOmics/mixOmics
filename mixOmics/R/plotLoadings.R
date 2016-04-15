# Copyright (C) 2015
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France

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


# ---------------------------------------------------------

# to do later
# create a S3 function?



plotLoadings = function(object,
                        block = 1, #single value
                        comp = 1,
                        show.ties = TRUE,
                        col = NULL,
                        ndisplay = NULL,
                        cex.name = 0.7,
                        cex.legend = 0.8,
                        name.var = NULL,
                        complete.name.var = FALSE,
                        main = NULL,
                        xlim = NULL)
{
    
    # -------------------
    # input checks
    # ------------------
    
    if (is.null(object$loadings))
    stop("'plotLoadings' should be used on object for which object$loadings is present.")
    

    # cex
    # --
    if (cex.name <=  0)
    cex.name = 0.7
    
    if (cex.legend <=  0)
    cex.name = 0.8
    
    ncomp = object$ncomp
    
    if (length(block) > 1)
    stop("Incorrect value for 'block', a single value (numeric or character) is expected", call. = FALSE)
    #xlim
    #-----
    if (!is.null(xlim))
    {
        if (length(xlim) !=2 & !is.numeric(xlim))
        stop("'xlim' needs to be a numeric vector of length 2")
    }

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


    ##selectvar
    selected.var = selectVar(object, comp = comp, block = block) # gives name and values of the blocks in 'block'
    name.selected.var = selected.var[[1]]$name
    value.selected.var = selected.var[[1]]$value

    # ndisplay
    # ------
    # if null set by default to all variables from selectVar
    if (is.null(ndisplay))
    {
        ndisplay = length(name.selected.var)
    } else if(ndisplay > length(name.selected.var)) {
        message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var))
        ndisplay = length(name.selected.var)
    }
    
    name.selected.var = name.selected.var[1:ndisplay]
    value.selected.var = value.selected.var[1:ndisplay,]

    #comp
    # ----
    if (any(class(object) %in% c("pls","spls")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
    {
        ncomp = object$ncomp
        object$X = list(X = object$X, Y = object$Y) # so that the data is in object$X, either it's a pls or block approach
    } else {
        ncomp = object$ncomp[block]
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


    # in case the user does not want to show the legend, then margins can be reduced

    # barplot with contributions
    layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.1), TRUE)

    par(mar = c(5, min(9, max(sapply(colnames.X, nchar))), 4, 0))
    mp = barplot(contrib$importance, horiz = T, las = 1, col = col, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(contrib),
    cex.names = cex.name, cex.axis = 0.7, beside = TRUE,border = NA, xlim=xlim)
    if (is.null(main))
    {
        title(paste0('Loadings on comp ', comp, "\nBlock '", names.block,"'"))
    } else {
        title(paste(main))
    }
    
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    par(mfrow = c(1,1))
    
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  = 
    # return the contribution matrix
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  = 
    return(invisible(list(contrib = contrib)))
}
