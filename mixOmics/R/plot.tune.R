#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#
# created: 2010
# last modified: 19-04-2016
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
function(object)
{
    
    error <- object$error.rate
    select.keepX <- object$choice.keepX[colnames(error)]
    first.keepX <- names(object$choice.keepX[1])
    
    legend=NULL
    
    matplot(error, type = "l", axes = FALSE, lwd = 2, lty = 1,
    xlab = "number of selected genes", ylab = "error rate",
    col = 1:length(select.keepX)    )
    
    for(i in 1:length(select.keepX))
    {
        # store coordinates of chosen keepX
        index = which(rownames(error) == select.keepX[i])
        
        # choseen keepX:
        points(index, error[index,i], col = i, lwd=2, cex=1.5)
        if(names(select.keepX[i])==first.keepX)
        legend=first.keepX
        else
        legend=c(legend,paste(first.keepX," to ",names(select.keepX[i])))
    }
    axis(1, c(1:length(rownames(error))),labels=rownames(error))
    axis(2)
    legend("topright", lty = 1, lwd = 2, horiz = FALSE, col = 1:length(select.keepX),
    legend = legend)
    
}


