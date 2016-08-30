#############################################################################################################
# Author:
#   Francois Bartolo, 
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 29-08-2016
# last modified: 30-08-2016
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

# --------------------------------------------------------------------------------------
# Internal helpers functions to run "plot.perf" functions
# --------------------------------------------------------------------------------------


internal_graphic.perf<- function (error.rate, error.rate.sd, overlay, type, measure, dist, legend, xlab, ylab, ...)
{
    # error.rate is a list [[measure]]
    # error.rate[[measure]] is a matrix of dist columns and ncomp rows
    
    # concatenation of the error in a single matrix
    error.rate.concat = matrix(nrow = nrow(error.rate[[1]]), ncol = 0)
    for(mea in measure)
    {
        temp = error.rate[[mea]][, dist, drop = FALSE]
        colnames(temp) = paste(mea, colnames(temp),sep="_")
        error.rate.concat = cbind(error.rate.concat, temp)
    }
    
    if(!is.null(error.rate.sd))
    {
        error.rate.sd.concat = matrix(nrow = nrow(error.rate.sd[[1]]), ncol = 0)
        for(mea in measure)
        {
            temp = error.rate.sd[[mea]][, dist, drop = FALSE]
            colnames(temp) = paste(mea, colnames(temp),sep="_")
            error.rate.sd.concat = cbind(error.rate.sd.concat, temp)
        }
    } else {
        error.rate.sd.concat = NULL
    }
   
    
    if(overlay == "all")
    {
        out<-matplot(error.rate.concat, type = type, lty = rep(c(1:length(measure)), each = length(dist)), col = rep(color.mixo(1:length(dist)), length(measure)), 
        lwd = 2, xlab = xlab, ylab = ylab, axes=FALSE, ylim=c(min(error.rate.concat), max(error.rate.concat)))
        
        axis(1, 1:nrow(error.rate.concat), rownames(error.rate.concat))
        axis(2)
        
        if(legend == "vertical")
        {
            legend('topright', legend = c(measure, dist), lty = c(1:length(measure), rep(NA, length(dist))), 
            pch = c(rep(NA, length(measure)), rep(16, length(dist))), col = c(rep('black', length(measure)), color.mixo(1:length(dist))), ncol = 1, lwd = 2)
        } else if(legend == "horizontal") {
            legend('topright', legend = c(measure, "" , dist), lty = c(1:length(measure), rep(NA, (length(dist)+1))), 
            pch = c(rep(NA, (length(measure)+1)), rep(16, length(dist))), col = c(rep('black', length(measure)), NA, color.mixo(1:length(dist))), ncol = 2, lwd = 2)
            
        }
        if(!is.null(error.rate.sd.concat))
        {
            for(j in 1:ncol(error.rate.concat))
            plot_error_bar(error.rate.concat[, j], uiw=error.rate.sd.concat[, j], add=T, col = rep(color.mixo(rep(j,each=nrow(error.rate.concat))), length(measure)))#, ...)
        }
    } else if(overlay == "measure") {
        
        # overlaying all measure on one graph, a graph per distance
        par(mfrow=c(1, length(dist)))
        for(di in dist)
        {
            new_mat.error=error.rate.concat[, grep(di, colnames(error.rate.concat))]
            
            out<-matplot(new_mat.error, type = type, lty = c(1:length(measure)),
            col = rep(color.mixo(which(di == dist))),
            lwd = 2, xlab = xlab, ylab = ylab, axes=FALSE, ylim=c(min(error.rate.concat), max(error.rate.concat)))
            
            axis(1, 1:nrow(error.rate.concat), rownames(error.rate.sd.concat))
            axis(2)
            
            if(legend == "vertical")
            {
                legend('topright', legend = measure, lty = 1:length(measure), col = rep(color.mixo(which(di == dist))), ncol = 1, lwd = 2)
            } else if(legend == "horizontal") {
                legend('topright', legend = c(measure), lty = 1:length(measure), col = rep(color.mixo(which(di == dist))), ncol = 2, lwd = 2)
            }
            if(!is.null(error.rate.sd.concat))
            {
                if(is.matrix(new_mat.error))
                {
                    for(col in 1:ncol(new_mat.error))
                    {
                        new_sd.error=error.rate.sd.concat[, grep(di, colnames(error.rate.sd.concat))]
                        plot_error_bar(new_mat.error[, col], uiw=new_sd.error[, col], add=T)#, ...)
                    }
                } else {
                    plot_error_bar(new_mat.error, uiw=error.rate.sd.concat, add=T)#, ...)
                }
            }
            #add title for each subgraph
            title(di)
        }
        
    } else if(overlay == "dist") {
        
        # overlaying all distance on one graph, a graph per measure

        par(mfrow=c(1, length(measure)))
        for(mea in measure)
        {
            new_mat.error=error.rate.concat[, grep(mea, colnames(error.rate.concat))]
            
            out<-matplot(new_mat.error, type = type, lty = c(1:length(dist)),
            col = rep(color.mixo(1:length(dist))),
            lwd = 2, xlab = xlab, ylab = ylab, axes=FALSE, ylim=c(min(error.rate.concat), max(error.rate.concat)))
            
            axis(1, 1:nrow(error.rate.concat), rownames(error.rate.concat))
            axis(2)
            
            if(legend == "vertical")
            {
                legend('topright', legend = dist, lty = 1:length(dist), col = rep(color.mixo(1:length(dist))), ncol = 1, lwd = 2)
            } else if(legend == "horizontal") {
                legend('topright', dist, lty = 1:length(dist), col = rep(color.mixo(1:length(dist))), ncol = 2, lwd = 2)
            }
            if(!is.null(error.rate.sd.concat))
            {
                if(is.matrix(new_mat.error))
                {
                    for(col in 1:ncol(new_mat.error))
                    {
                        new_sd.error=error.rate.sd.concat[, grep(mea, colnames(error.rate.sd.concat))]
                        plot_error_bar(new_mat.error[, col], uiw=new_sd.error[, col], add=T)#, ...)
                    }
                } else {
                    plot_error_bar(new_mat.error, uiw=error.rate.sd.concat, add=T)#, ...)
                }
            }
            title(mea)
        }
        
    }
    
}


plot_error_bar <- function (x, y = NULL, uiw, add=FALSE, col = "black", ...)
{
    sfrac = 0.01
    gap = 0
    slty = par("lty")
    pt.bg = par("bg")
    arglist <- list()#...)
    liw = uiw
    y <- as.numeric(x)
    x <- seq(along = x)
    z <- y
    ui <- z + uiw
    li <- z - liw
    arglist$xlab <- deparse(substitute(x))
    arglist$ylab <- deparse(substitute(y))
    arglist$ylim <- range(c(y, ui, li), na.rm = TRUE)
    arglist$xlim <- range(c(x, ui, li), na.rm = TRUE)
    scol <- par("col")
    plotpoints <- TRUE
    ul <- c(li, ui)
    pin <- par("pin")
    usr <- par("usr")
    x.to.in <- pin[1]/diff(usr[1:2])
    y.to.in <- pin[2]/diff(usr[3:4])
    gap <- rep(gap, length(x)) * diff(par("usr")[3:4])
    smidge <- par("fin")[1] * sfrac
    nz <- abs(li - pmax(y - gap, li)) * y.to.in > 0.001
    scols <- col#rep(scol, length.out = length(x))[nz]
    
    arrow.args <- c(list(lty = slty, angle = 90, length = smidge, 
    code = 1, col = scols), clean.args(arglist, arrows, 
    exclude.other = c("col", "lty", "axes")))
    
    do.call("arrows", c(list(x[nz], li[nz], x[nz], pmax(y -
    gap, li)[nz]), arrow.args))
    
    nz <- abs(ui - pmin(y + gap, ui)) * y.to.in > 0.001
    #scols <- rep(scol, length.out = length(x))[nz]
    arrow.args$col <- scols
    
    do.call("arrows", c(list(x[nz], ui[nz], x[nz], pmin(y +
    gap, ui)[nz]), arrow.args))
    
    do.call("points", c(list(x, y, bg = pt.bg, col =col), clean.args(arglist,
    points, exclude.other = c("xlab", "ylab", "xlim", 
    "ylim", "axes"))))
    
    invisible(list(x = x, y = y))
}


# from plotrix: 3.0-6
clean.args<-function(argstr, fn, exclude.repeats=FALSE, exclude.other=NULL, 
dots.ok=TRUE) {
    
    fnargs<-names(formals(fn))
    if(length(argstr) > 0 && !("..." %in% fnargs && dots.ok)) {
        badargs<-names(argstr)[!sapply(names(argstr), "%in%", c(fnargs, ""))]
        for(i in badargs) argstr[[i]]<-NULL
    }
    if(exclude.repeats) {
        ntab<-table(names(argstr))
        badargs<-names(ntab)[ntab > 1 & names(ntab) != ""]
        for (i in badargs) argstr[[i]]<-NULL
    }
    for(i in exclude.other) argstr[[i]]<-NULL
    argstr
}
