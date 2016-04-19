#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 19-04-2016
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



#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS-DA, SPLS-DA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#


#plotLoadings.mint.pls      =
#plotLoadings.mint.spls     =
plotLoadings.mint.plsda    =
plotLoadings.mint.splsda   =

function(object,
contrib = NULL,  # choose between 'max" or "min", NULL does not color the barplot
method = "mean", # choose between 'mean" or "median"
study = "all",
comp = 1,
show.ties = TRUE,
col.ties = "white",
ndisplay = NULL,
cex.name = 0.7,
cex.legend = 0.8,
name.var = NULL,
complete.name.var = FALSE,
legend = TRUE,
legend.color = NULL,
main = NULL,
subtitle,
legend.title = 'Outcome',
plot = TRUE,
layout = NULL,
size.title = rel(1.8),
size.subtitle = rel(1.4),
...
) {
    
    # -------------------
    # input checks
    # ------------------
    if (is.null(object$loadings))
    stop("'plotLoadings' should be used on object for which object$loadings is present.")
    
    
    
    # what I want is to modify the input and call plotLoadings.pls and plotLoadings.splsda where blocks are now studies
    # do not forget to change object$names$block in levels(object$study) and it should work, see you tomorrow
    
    # block
    # --
    if (missing(study))
    {
        block = "X"
    }
    
    if (any(class(object) %in% c("plsda", "splsda")) & (!all(block %in% c(1,"X")) | length(block) > 1 ))
    stop("'block' can only be 'X' or '1' for plsda and splsda object")
    
    if (any(class(object) %in% c("plsda", "splsda")))
    object$indY=2
    
    
    if(is.numeric(block))
    {
        if(any(block>length(object$names$blocks[-object$indY])))
        stop("'block' needs to be lower than the number of blocks in the fitted model, which is ",length(object$names$blocks)-1)
        
    }else if(is.character(block) & any(is.na(match(block,object$names$blocks[-object$indY])))) {
        stop("Incorrect value for 'block', 'block' should be among the blocks used in your object: ", paste(object$names$blocks[-object$indY],collapse=", "), call. = FALSE)
    }
    
    
    # contrib
    # --
    
    # if contrib is NULL, then we switch to the classical plotLoadings (without contribution/colors)
    if(is.null(contrib))
    {
        if(plot)
        {
            plotLoadings.pls(object = object, block = block, comp = comp, ndisplay = ndisplay,
            cex.name = cex.name,
            cex.legend = cex.legend,
            name.var = name.var,
            complete.name.var = complete.name.var,
            main = main,
            subtitle = subtitle,
            xlim = xlim,
            layout = layout,
            size.title = size.title,
            size.subtitle = size.subtitle)
        } else {
            stop("'contrib' is NULL and 'plot' is FALSE => no results to show", call. = FALSE)
        }
        # stop the script without error message
        # blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
        # stop(simpleError(blankMsg))
    } else {
        
        
        # layout
        # --
        if(plot == TRUE)
        {
            opar = par(no.readonly = TRUE)
            reset.mfrow = FALSE # if set to TRUE, the algorithm ends up with  par(mfrow=reset.mfrow)
            nResp = length(block) + length(block) * legend  #number of blocks *2 if legend is plotted
            if (is.null(layout))
            {
                # check if there are enough plots in mfrow
                omfrow = par("mfrow")
                available.plots = prod(omfrow)
                if (available.plots<nResp) # if not enough plots available, we create our new plot
                {
                    nCols = min(c(3 + legend, ceiling(nResp)))
                    nRows = min(c(3 + legend, ceiling(nResp/nCols)))
                    layout = c(nRows, nCols)
                    
                    layout(matrix(1 : min(nCols * nRows, nResp), nRows, nCols, byrow=TRUE),rep(c(0.7,0.7 -0.4*legend),nCols/(1+legend)))
                    
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
            
        }
        
        
        
        #layout(matrix(c(1, 2,3,4,5,6), 3, 2, byrow = TRUE), c(0.7, 0.3), TRUE)
        
        #} else {
        #layout(matrix(c(1, 2), 1, 2, byrow = TRUE), c(0.7, 0.1), TRUE)
        
        
        class.object = class(object)
        
        # method
        # ----
        if (length(method) !=1 || !method %in% c("mean","median"))
        {
            method = "median"
            warning("'method' should be either 'mean' or 'median', set to 'median' by default")
        }
        # cex
        # --
        if (cex.name <= 0)
        cex.name = 0.7
        
        if (cex.legend <= 0)
        cex.name = 0.8
        
        ncomp = object$ncomp
        
        
        #title
        #-----
        if (!is.null(main) & !is.character(main))
        warning('main needs to be of type character')
        
        
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
            }else if (ndisplay > length(name.selected.var)) {
                message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var))
                ndisplay.temp = length(name.selected.var)
            } else {
                ndisplay.temp = ndisplay
            }
            
            name.selected.var = name.selected.var[1:ndisplay.temp]
            value.selected.var = value.selected.var[1:ndisplay.temp, ]
            
            #comp
            # ----
            if (any(class(object) %in% c("plsda","splsda")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
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
            if (!is.null(name.var))
            {
                if (length(name.var) != ncol(X))
                stop("name.var should be a vector of length ", ncol(X))
                
                colnames.X = as.character(name.var[ind.match]) # get the
            } else {
                colnames.X = as.character(colnames(X))[ind.match]
            }
            X = X[, name.selected.var] #reduce the problem to ndisplay
            
            #completing colnames.X by the original names of the variables when missing
            if (complete.name.var == TRUE)
            {
                ind = which(colnames.X == "")
                if (length(ind)>0)
                colnames.X[ind] = colnames(X)[ind]
            }
            
            Y = object$Y #v6: all $Y are factors for DA methods
            
            #legend.color
            #-----
            if (!is.null(legend.color) & (length(legend.color) != nlevels(Y)))
            {
                warning('legend.color must be the same length than the number of group, by default set to default colors')
                legend.color = color.mixo(1:10)  # by default set to the colors in color.mixo (10 colors)
            }
            if (is.null(legend.color))
            legend.color = color.mixo(1:10)[1:nlevels(Y)] # by default set to the colors in color.mixo (10 colors)
            
            if (col.ties%in%legend.color[1:nlevels(Y)])
            stop("'col.ties' should not be in 'legend.color'")
            
            # --------------------
            # end check inputs
            # ---------------------
            
            # ==================================================
            # First, calculate the contribution of the loading
            # =================================================
            
            # Start: Initialisation
            which.comp = method.group = list()
            which.contrib = data.frame(matrix(FALSE, ncol = nlevels(Y) + 2, nrow = ndisplay.temp,
            dimnames = list(name.selected.var, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
            # End: Initialisation
            
            # calculate the max.method per group for each variable, and identifies which group has the max max.method
            for(k in 1:ncol(X))
            {
                method.group[[k]] = tapply(X[, k], Y, method, na.rm=TRUE) #method is either mean or median
                # determine which group has the highest mean/median
                which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == get(contrib)((method.group[[k]])) # contrib is either min or max
            }
            
            
            # if ties, we set the color to white
            which.contrib$Contrib = apply(which.contrib, 1, function(x)
            {
                if (length(which(x)) > 1)
                {
                    return(col.ties)
                } else { # otherwise we use legend color provided
                    return(legend.color[1 : nlevels(Y)][which(x)])
                }
            })
            
            # we also add an output column indicating the group that is max
            which.contrib$GroupContrib = apply(which.contrib[, 1:(nlevels(Y))], 1, function(x)
            {
                if (length(which(x)) > 1)
                {
                    return("tie")
                } else {
                    return(levels(Y)[which(x)])
                }
            })
            
            method.group = do.call(rbind, method.group)
            
            
            df = data.frame(method.group, which.contrib, importance = value.selected.var)
            # End contribution calculation
            
            # =====================================
            # then determine the colors/groups matching max contribution
            # =======================================
            
            # when working with sparse counts in particular and using the median to measure contribution
            # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
            if (show.ties == FALSE)
            {
                df = df[!df$Contrib %in% col.ties, ]
                colnames.X = rownames(df)
            }
            
            #display barplot with names of variables
            
            # ==================================
            # represent the barplot
            # ==================================
            #added condition if all we need is the contribution stats
            if (plot)
            {    #opar <- par(no.readonly = TRUE)
                
                if (!is.null(main) & length(block) > 1)
                {
                    par(mar = c(4, max(6, max(sapply(colnames.X, nchar))/2), 6, 2))
                } else {
                    par(mar = c(4, max(6, max(sapply(colnames.X, nchar))/2), 4, 2))
                }
                
                
                mp = barplot(df$importance, horiz = T, las = 1, col = df$Contrib, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(df),
                cex.names = cex.name, cex.axis = 0.7, beside = TRUE, border = NA)

                
                if ( length(block) == 1 & is.null(main) )
                {
                    title(paste0('Contribution on comp ', comp), line=1, cex.main = size.subtitle)
                } else if (length(block) == 1) {
                    title(paste(main), line=1, cex.main = size.subtitle)
                } else if ((length(block) > 1 & missing(subtitle))) {
                    title(paste0('Contribution on comp ', comp, "\nStudy '", names.block,"'"), line=1, cex.main = size.subtitle)
                } else if (length(block) > 1 & !missing(subtitle)) {
                    title(paste(subtitle[i]), line=1, cex.main = size.subtitle)
                }
                
                
                if (legend)
                {
                    par(mar = c(5, 0, 4, 3) + 0.1)
                    plot(1,1, type = "n", axes = FALSE, ann = FALSE)
                    legend(0.8, 1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19,
                    title = paste(legend.title),
                    cex = cex.legend)
                    
                    
                }
                
                
                
                
                #par(opar)
                
            } # end if plot
            
        }
        
        
        if(plot)
        {
            # legend
            if (length(block) > 1 & !is.null(main))
            title(main, outer=TRUE, line = -2, cex.main = size.title)
        }
        if (plot == TRUE && reset.mfrow)
        par(opar)#par(mfrow = omfrow)
        
        # ===================================
        # return the contribution matrix
        # ===================================
        return(invisible(list(contrib = df)))
    }# end contrib missing
}
