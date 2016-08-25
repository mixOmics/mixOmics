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


plot.tune.splsda = #plot.spca <- plot.ipca <- plot.sipca <-
function(x, ...)
{
    
    error <- x$error.rate
    select.keepX <- x$choice.keepX[colnames(error)]
    first.keepX <- names(x$choice.keepX[1])
    
    legend=NULL
    measure = x$measure
    
    if (length(select.keepX) < 10)
    {
        #only 10 colors in color.mixo
        col.per.comp = color.mixo(1:length(select.keepX))
    } else {
        #use color.jet
        col.per.comp = color.jet(length(select.keepX))
    }
    
    if(measure == "overall")
    {
         ylab = "Classification error rate"
    } else if (measure == "BER")
    {
        ylab = "Balanced error rate"
    }
   

    matplot(rownames(error),error, type = "l", axes = TRUE, lwd = 2, lty = 1, log = "x",
    xlab = "Number of selected genes", ylab = ylab,
    col = col.per.comp)
    
    for(i in 1:length(select.keepX))
    {
        # store coordinates of chosen keepX
        index = which(rownames(error) == select.keepX[i])
        
        # choseen keepX:
        points(rownames(error)[index], error[index,i], col = col.per.comp[i], lwd=2, cex=1.5)
        if(names(select.keepX[i])==first.keepX)
        legend=first.keepX
        else
        legend=c(legend,paste(first.keepX," to ",names(select.keepX[i])))
    }
    #axis(1, c(1:length(rownames(error))),labels=rownames(error))
    #axis(2)
    legend("topright", lty = 1, lwd = 2, horiz = FALSE, col = col.per.comp,
    legend = legend)
    
}


