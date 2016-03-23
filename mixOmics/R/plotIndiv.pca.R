#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 16-03-2016
# last modified: 23-03-2016
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


plotIndiv <-
function(object, ...) UseMethod("plotIndiv")


#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#


plotIndiv.pca=

function(object,
comp = NULL,
ind.names = TRUE,
group,  # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
col.per.group,
style="ggplot2", # can choose between graphics,3d, lattice or ggplot2
plot.ellipse = FALSE,
ellipse.level = 0.95,
plot.centroid=FALSE,
plot.star=FALSE,
main=NULL,
add.legend=FALSE,
X.label = NULL,
Y.label = NULL,
Z.label = NULL,
abline.line = FALSE,
xlim = NULL,
ylim = NULL,
col,
cex,
pch,
alpha=0.2,
axes.box = "box",
layout=NULL,
...
)
{
    
    
    blocks = "X"
    rep.space= "X-variate"
    
    check = check.input.plotIndiv(object=object,comp = comp,rep.space = rep.space,blocks = blocks, ind.names = ind.names,
    style=style, plot.ellipse = plot.ellipse, ellipse.level = ellipse.level, plot.centroid=plot.centroid,
    plot.star=plot.star, add.legend=add.legend,X.label = X.label,Y.label = Y.label,Z.label = Z.label,abline.line = abline.line,
    xlim = xlim,ylim = ylim,alpha=alpha,axes.box = axes.box)
    
    # retrieve outputs from the checks
    axes.box=check$axes.box
    comp=check$comp
    rep.space=check$rep.space
    xlim=check$xlim
    ylim=check$ylim
    ind.names=check$ind.names
    display.names=check$display.names


    #-- Get variates
    x=y=z=list()
    x[[1]] = object$x[, comp[1]]
    y[[1]] = object$x[, comp[2]]
    if(style=="3d") z[[1]] = object$x[, comp[3]]
    
    
    #-- Variance explained on X, Y and Z labels

    if (style == "3d")
    {
        inf=object$explained_variance[c(comp[1],comp[2],comp[3])]#c((object$sdev[comp[1]])^2/object$var.tot,(object$sdev[comp[2]])^2/object$var.tot,(object$sdev[comp[3]]^2)/object$var.tot)
        inf = round(inf,2)
    } else {
        inf = object$explained_variance[c(comp[1],comp[2])]#c((object$sdev[comp[1]])^2/object$var.tot,(object$sdev[comp[2]])^2/object$var.tot)
        inf = round(inf,2)}
    

    if (is.null(X.label))
    {
        X.label = paste("PC", comp[1],sep='')
        percentage=paste0(inf[1]*100,"% expl. var")
        X.label = paste(X.label,percentage,sep=": ")
    }
    if (is.null(Y.label))
    {
        Y.label = paste("PC", comp[2],sep='')
        percentage=paste0(inf[2]*100,"% expl. var")
        Y.label = paste(Y.label,percentage,sep=": ")
    }
    if (is.null(Z.label)&&style=="3d")
    {
        Z.label = paste("PC", comp[3],sep='')
        percentage=paste0(inf[3]*100,"% expl. var")
        Z.label = paste(Z.label,percentage,sep=": ")
    }


    n=nrow(object$X)

    # create data frame df that contains (almost) all the ploting information
    out=shape.input.plotIndiv(object=object,n=n,blocks = blocks,x=x,y=y,z=z,ind.names = ind.names,group,col.per.group=col.per.group,
    style=style,study="all",plot.ellipse = plot.ellipse,ellipse.level = ellipse.level,
    plot.centroid=plot.centroid,plot.star=plot.star,main=main,xlim = xlim,ylim = ylim,
    col=col,cex=cex,pch=pch,display.names=display.names)
    #-- retrieve outputs
    df=out$df
    df.ellipse=out$df.ellipse
    col.per.group=out$col.per.group
    main=out$main
    display.names=out$display.names
    xlim=out$xlim
    ylim=out$ylim
    missing.col=out$missing.col
    plot.ellipse=out$plot.ellipse
    plot.centroid=out$plot.centroid
    plot.star=out$plot.star
    
    save(list=ls(),file="temp.Rdata")
    
    #call plot module (ggplot2,lattice,graphics,3d)
    internal_graphicModule(df=df,plot.centroid=plot.centroid,col.per.group=col.per.group,main=main,X.label=X.label,Y.label=Y.label,Z.label=Z.label,
    xlim=xlim,ylim=ylim,class.object=class(object),display.names=display.names,add.legend=add.legend,
    abline.line=abline.line,plot.star=plot.star,plot.ellipse=plot.ellipse,df.ellipse=df.ellipse,style=style,layout=layout,missing.col=missing.col,
    axes.box=axes.box)
    
    return(invisible(list(df=df,df.ellipse=df.ellipse)))

}
