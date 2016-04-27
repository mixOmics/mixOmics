#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 15-04-2016
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
#-- Includes plotLoadings for PLS, sPLS, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA --#
#----------------------------------------------------------------------------------------------------------#


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
col = NULL,
ndisplay = NULL,
cex.name = 0.7,
name.var = NULL,
complete.name.var = FALSE,
main = NULL,
subtitle,
layout = NULL,
size.title = rel(2),
size.subtitle = rel(1.5),
border = NA,
...)
{
    
    # -- input checks
    check = check.input.plotLoadings(object = object, block = block, subtitle = subtitle, cex.name = cex.name, main = main, col = col, name.var = name.var)
    
    col = check$col
    cex.name = check$cex.name
    block = check$block


    # -- layout
    res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
    reset.mfrow = res$reset.mfrow
    opar = res$opar
    omar = par("mar") #reset mar at the end

    if (length(block) == 1 & !is.null(name.var))
    name.var = list(name.var = name.var)
    
    for (i in 1 : length(block))
    {
        res = get.loadings.ndisplay(object = object, comp = comp, block = block[i], name.var = name.var[[i]], complete.name.var = complete.name.var, ndisplay = ndisplay)
        X = res$X
        names.block = res$names.block
        colnames.X = res$colnames.X
        value.selected.var = res$value.selected.var

        df = data.frame(importance = value.selected.var) # contribution of the loading
        
        # barplot with contributions
        if (!is.null(main) & length(block) > 1)
        {
            par(mar = c(4, max(7, max(sapply(colnames.X, nchar))/3), 6, 2))
        } else {
            par(mar = c(4, max(7, max(sapply(colnames.X, nchar))/3), 4, 2))
        }

        mp = barplot(df$importance, horiz = T, las = 1, col = col, axisnames = TRUE, names.arg = colnames.X, #names.arg = row.names(df),
        cex.names = cex.name, cex.axis = 0.7, beside = TRUE, border = border)
        
        if ( (length(block) == 1 & is.null(main)) | (length(block) > 1 & missing(subtitle)))
        {
            title(paste0('Loadings on comp ', comp, "\nBlock '", names.block,"'"), line=1, cex.main = size.title)
        } else if (length(block) == 1) {
            title(paste(main), line=1, cex.main = size.title)
        } else if (length(block) > 1 & !missing(subtitle)) {
            title(paste(subtitle[i]), line=1, cex.main = size.subtitle)
        }
    }
    
    if (length(block) > 1 & !is.null(main))
    title(main, outer=TRUE, line = -2, cex.main = size.title)
    
    if (reset.mfrow)
    par(opar)#par(mfrow = c(1,1))

    par(mar = omar) #reset mar

    # return the contribution matrix
    return(invisible(list(contrib = df)))
}
