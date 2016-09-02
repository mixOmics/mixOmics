#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#
# created: 20-08-2016
# last modified: 25-08-2016
#
# Copyright (C) 2010
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

plot.tune<-function(x,...) NextMethod("plot")


plot.tune.splsda = #plot.tune.mint.splsda
function(x, optimal = TRUE, sd=TRUE,horiz=FALSE,overlay=TRUE,log="", ...)
{
    
    if (!is.logical(optimal))
    stop("'optimal' must be logical.", call. = FALSE)
    
    
    error <- x$error.rate
    if(sd & !is.null(x$error.rate.sd))
    {
        error.rate.sd = x$error.rate.sd
        ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
    } else {
        error.rate.sd = NULL
        ylim = range(error)
    }
    
    select.keepX <- x$choice.keepX[colnames(error)]
    comp.tuned = length(select.keepX)
    
    legend=NULL
    measure = x$measure
    
    if (length(select.keepX) < 10)
    {
        #only 10 colors in color.mixo
        col.per.comp = color.mixo(1:comp.tuned)
    } else {
        #use color.jet
        col.per.comp = color.jet(comp.tuned)
    }
    
    if(measure == "overall")
    {
        ylab = "Classification error rate"
    } else if (measure == "BER")
    {
        ylab = "Balanced error rate"
    }
    
    if(overlay)
    {
        matplot(rownames(error),error, type = "l", axes = TRUE, lwd = 2, lty = 1, log = log,
        xlab = "Number of selected genes", ylab = ylab,
        col = col.per.comp, ylim = ylim)
        
        if(optimal)
        {
            for(i in 1:comp.tuned)
            {
                # store coordinates of chosen keepX
                index = which(rownames(error) == select.keepX[i])
                # choseen keepX:
                points(rownames(error)[index], error[index,i], col = col.per.comp[i], lwd=2, cex=3, pch = 18)
            }
        }
        
        if(!is.null(error.rate.sd))
        {
            for(j in 1:ncol(error))
            plot_error_bar(x = as.numeric(names(error[, j])), y =error[, j] , uiw=error.rate.sd[, j], add=T, col = color.mixo(rep(j,nrow(error))))#, ...)
        }
        
        
        for(test in 1:comp.tuned)
        {
            if(j==1 && strsplit(colnames(error),"p")[[1]][2]==1)
            legend=c(legend,strsplit(colnames(error),"p")[[1]][2])
            else
            legend=c(legend,paste("1 to",strsplit(colnames(error),"p")[[test]][2]))
            
        }
        
        legend("topright", lty = 1,  horiz = horiz, col = col.per.comp,
        legend = legend,title="Component :",seg.len = 1)
    }
    else
    {
        def.par <- par(no.readonly = TRUE)
        par(mfrow=c(1,comp.tuned))
        for(i in 1:comp.tuned)
        {
            matplot(rownames(error),error[,i], type = "l", axes = TRUE, lwd = 2, lty = 1, log = log,
            xlab = "Number of selected genes", ylab = ylab,
            col = col.per.comp[i], ylim = ylim)
            
            if(optimal)
            {
                # store coordinates of chosen keepX
                index = which(rownames(error) == select.keepX[i])
                # choseen keepX:
                points(rownames(error)[index], error[index,i], col = col.per.comp[i], lwd=2, cex=3, pch = 18)
                
            }
            
            if(!is.null(error.rate.sd))
            {
                plot_error_bar(x = as.numeric(rownames(error)), y =error[, i] , uiw=error.rate.sd[, i], add=T, col = color.mixo(rep(i,nrow(error))))#, ...)
            }
            
            if(strsplit(colnames(error),"p")[[i]][2]==1)
            title(paste("Component 1"))
            else
            title(paste("Component 1 to",strsplit(colnames(error),"p")[[i]][2]))
            
        }
        
        par(def.par)
    }
    
}


