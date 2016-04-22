#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 19-04-2016
# last modified: 20-04-2016
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
#-- Includes plotLoadings for mint.pls and mint.spls --#
#----------------------------------------------------------------------------------------------------------#


plotLoadings.mint.pls    =
plotLoadings.mint.spls   =

function(object, block, 
study = "all",
comp = 1,
col = NULL,
ndisplay = NULL,
cex.name = 0.7,
name.var = NULL,
complete.name.var = FALSE,
main = NULL,
subtitle,
plot = TRUE,
layout = NULL,
size.title = rel(1.8),
size.subtitle = rel(1.4),
border = NA,
...
) {
    
    # what I want is to modify the input and call plotLoadings.pls and plotLoadings.splsda where blocks are now studies
    # do not forget to change object$names$block in levels(object$study) and it should work, see you tomorrow
    
    if(any(study == "all"))
    {
        # if study == "all" then we plot the results on the concatenated data, thus direct call to plotLoadings.plsda
        plotLoadings.pls(object = object, contrib = contrib, method = method, block = "X", comp = comp, ndisplay = ndisplay,
        cex.name = cex.name,
        name.var = name.var,
        complete.name.var = complete.name.var,
        main = main,
        subtitle = subtitle,
        xlim = xlim,
        layout = layout,
        size.title = size.title,
        size.subtitle = size.subtitle,
        border = border,
        col = "white")
        
    } else {
        # if study != "all" then we plot the results on each study

        # -- input checks
        check = check.input.plotLoadings(object = object, block = "X", subtitle = subtitle, main = main, col = col, cex.name = cex.name)
        
        col = check$col
        cex.name = check$cex.name
        block = check$block # "X"
        
        # swap block for study
        block = study
        
        # -- layout
        res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
        reset.mfrow = res$reset.mfrow
        opar = res$opar
        omar = par("mar") #reset mar at the end
        
        # method
        if (length(method) !=1 || !method %in% c("mean","median"))
        {
            method = "median"
            warning("'method' should be either 'mean' or 'median', set to 'median' by default")
        }
        
        # get the selected variables on the concatenated data
        res = get.loadings.ndisplay(object = object, comp = comp, block = "X", name.var = name.var, complete.name.var = complete.name.var, ndisplay = ndisplay)
        X = res$X
        colnames.X = res$colnames.X
        name.selected.var = res$name.selected.var
        value.selected.var = res$value.selected.var


        # swap loadings partial for loadings
        object$loadings.global = object$loadings
        object$loadings = object$loadings.partial$X
        object$names$block = levels(object$study)
        
        for (i in 1 : length(block))
        {
            value.selected.var = object$loadings.partial$X [[block[i]]][, comp] [name.selected.var]

            df = data.frame(importance = value.selected.var, color = col, stringsAsFactors = FALSE) # contribution of the loading
           
            #display barplot with names of variables
            #added condition if all we need is the contribution stats
            if (!is.null(main) & length(block) > 1)
            {
                par(mar = c(4, max(6, max(sapply(colnames.X, nchar))/2), 6, 2))
            } else {
                par(mar = c(4, max(6, max(sapply(colnames.X, nchar))/2), 4, 2))
            }

            mp = barplot(df$importance, horiz = T, las = 1, col = df$color, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(df),
            cex.names = cex.name, cex.axis = 0.7, beside = TRUE, border = border)
            
            if ( length(block) == 1 & is.null(main) )
            {
                title(paste0('Contribution on comp ', comp), line=1, cex.main = size.subtitle)
            } else if (length(block) == 1) {
                title(paste(main), line=1, cex.main = size.subtitle)
            } else if ((length(block) > 1 & missing(subtitle))) {
                title(paste0('Contribution on comp ', comp, "\nStudy '", block[i],"'"), line=1, cex.main = size.subtitle)
            } else if (length(block) > 1 & !missing(subtitle)) {
                title(paste(subtitle[i]), line=1, cex.main = size.subtitle)
            }
            
        }
        
        if (length(block) > 1 & !is.null(main))
        title(main, outer=TRUE, line = -2, cex.main = size.title)
        
        if (reset.mfrow)
        par(opar)#par(mfrow = omfrow)
        
        par(mar = omar) #reset mar
        
        
    }
}
