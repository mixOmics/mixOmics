#############################################################################################################
# Authors:
#   Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
#   Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2011
# last modified: 19-04-2016
#
# Copyright (C) 2011
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

plot.perf<-function(x,...) NextMethod("plot")


# PLS object
# ----------------------

plot.perf.pls.mthd <-plot.perf.spls.mthd <-
function (x,
criterion = c("MSEP", "RMSEP", "R2", "Q2"),
xlab = "number of components",
ylab = NULL,
LimQ2 = 0.0975,
LimQ2.col = "darkgrey",
cTicks = NULL,
layout = NULL,
...)
{
    
        if (!any(criterion %in% c("MSEP", "RMSEP", "R2", "Q2")) || missing(criterion))
        stop("Choose a validation criterion: MSEP, RMSEP, R2 or Q2.")
        y = switch(criterion, MSEP = x$MSEP, RMSEP = sqrt(x$MSEP), R2 = x$R2, Q2 = x$Q2)
        
        Q2.total = NULL
        if ((criterion == "Q2") & is.list(y)) {
            Q2.total = y$Q2.total
            y = y$variables
        }
        
        if (is.null(ylab))
        ylab = switch(criterion, MSEP = "MSEP", RMSEP = "RMSEP",
        R2 = expression(R^~2), Q2 = expression(Q^~2))
        
        nResp = nrow(y)  # Number of response variables
        nComp = ncol(y)  # Number of components
        
        if (nResp > 1) {
            if (is.null(layout)) {
                nRows = min(c(3, nResp))
                nCols = min(c(3, ceiling(nResp / nRows)))
                layout = c(nRows, nCols)
            }
            else {
                if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
                stop("'layout' must be a numeric vector of length 2.")
                nRows = layout[1]
                nCols = layout[2]
            }
            
            if (nRows * nCols < nResp) devAskNewPage(TRUE)
            ynames = rownames(y)
        }
        else {
            ynames = "Y"
        }
        
        val = comps = vector("numeric")
        varName = vector("character")
        
        for (i in 1:nResp) {
            val = c(val, y[i, ])
            comps = c(comps, 1:nComp)
            varName = c(varName, rep(ynames[i], nComp))
        }
        
        df = data.frame(val = val, comps = comps, varName = varName)
        if (is.null(cTicks)) cTicks = 1:ncol(y)
        yList = list(relation = "free")
        
    
    if (criterion == "Q2")
    {
        plt = xyplot(val ~ comps | varName, data = df, xlab = xlab, ylab = ylab,
        scales = list(y = yList, x = list(at = cTicks)),
        as.table = TRUE, layout = layout,
        panel = function(x, y) {
            if (LimQ2.col != "none") panel.abline(h = LimQ2, col = LimQ2.col)
            panel.xyplot(x, y, ...)})
        plot(plt)
        
        if (!is.null(Q2.total)) {
            devAskNewPage(TRUE)
            Q2.df = data.frame(Q2 = Q2.total, comps = 1:nComp, varName = rep("Total", nComp))
            xyplot(Q2 ~ comps | varName, data = Q2.df, xlab = xlab, ylab = ylab,
            scales = list(y = yList, x = list(at = cTicks)), as.table = TRUE,
            panel = function(x, y) {
                if (LimQ2.col != "none") panel.abline(h = LimQ2, col = LimQ2.col)
                panel.xyplot(x, y, ...)})
        }
    }
    else {
        xyplot(val ~ comps | varName, data = df, xlab = xlab, ylab = ylab,
        scales = list(y = yList, x = list(at = cTicks)),
        as.table = TRUE, layout = layout, ...)
    }
  

  
}

# PLSDA object
# ----------------------

plot.perf.plsda.mthd <-plot.perf.splsda.mthd <-
  function (x,
            dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
            measure = c("all","overall","BER"),
            type="l",
            xlab = NULL,
            ylab = NULL,
            overlay=c("all","measure"),
            legend=c("vertical","horizontal"),
            sd=TRUE,
            ...)
  {
      
      if (hasArg(pred.method))
        stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
      pred.method = NULL # to pass R CMD check
      
      
      if (any(measure == "all"))
        measure = names(x$error.rate)
      
      if (is.null(measure) || !any(measure %in% names(x$error.rate)))
        stop("'measure' should be among the ones used in your call to 'perf': ", paste(names(x$error.rate),collapse = ", "),".")
      
      if (any(dist == "all"))
      dist = colnames(x$error.rate[[1]])
      
      
      if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
        stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
      
      # KA changed
      #x = matrix(x[, pred.method], ncol = length(pred.method))
      mat.error.plsda=matrix(nrow=nrow(x$error.rate[[1]]),ncol=0)
      for(mea in measure)
        {
        mat.error.plsda=cbind(mat.error.plsda,x$error.rate[[mea]][, dist])
      }
    colnames(mat.error.plsda)=rep(dist,length(measure))
    
    sd.error.plsda=matrix(nrow=nrow(x$error.rate.sd[[1]]),ncol=0)
    for(mea in measure)
    {
      sd.error.plsda=cbind(sd.error.plsda,x$error.rate.sd[[mea]][, dist])
    }
    colnames(sd.error.plsda)=rep(dist,length(measure))
      
      if (is.null(ylab))
      {
          ylab = 'Classification error rate'
        
      }
      
      if (is.null(xlab))
      {
        xlab = 'PLSDA components'
        
      }
    def.par <- par(no.readonly = TRUE) 
    internal_graph_plot.perf(mat.error.plsda,sd.error.plsda, overlay, type,measure,dist,legend,xlab,ylab,sd=sd,  ...)
    par(def.par)
    # error.bar(out,as.vector(mat.error.plsda),as.vector(cbind(x$error.rate.sd$overall,x$error.rate.sd$BER)))
    
return(invisible(out))
    
  }

# mint.PLSDA object
# ----------------------

plot.perf.mint.plsda.mthd <-plot.perf.mint.splsda.mthd <-
  function (x,
            dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
            measure = c("all","overall","BER"),
            type="l",
            xlab = NULL,
            ylab = NULL,
            study="all",
            error.rate=c("average","study"),
            overlay=c("all","measure"),
            legend=c("vertical","horizontal"),
            ...)
  {
    
    if (hasArg(pred.method))
      stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
    pred.method = NULL # to pass R CMD check
    
    if(any(error.rate=="average"))
  {  if (any(measure == "all"))
      measure = c("BER","overall")
    
    if (is.null(measure) || !any(measure %in% c("BER","overall")))
      stop("'measure' should be among the ones used in your call to 'perf': ", paste(c("BER","overall"),collapse = ", "),".")
    
    
    if (any(dist == "all"))
      dist = colnames(x$global.error[[1]])
    
    
    if (is.null(dist) || !any(dist %in% colnames(x$global.error[[1]])))
      stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$global.error[[1]]),collapse = ", "),".")
    
    # KA changed
    #x = matrix(x[, pred.method], ncol = length(pred.method))
    mat.error.plsda=matrix(nrow=nrow(x$global.error[[1]]),ncol=0)
    for(mea in measure)
    {
      mat.error.plsda=cbind(mat.error.plsda,x$global.error[[mea]][, dist])
    }
    
    if (is.null(ylab))
    {
      ylab = 'Classification error rate'
      
    }
    
    if (is.null(xlab))
    {
      xlab = 'PLSDA components'
      
    }
    def.par <- par(no.readonly = TRUE) 
    internal_graph_plot.perf(mat.error.plsda,sd.error.plsda=NULL, overlay, type,measure,dist,legend,xlab,ylab,sd=FALSE,  ...)
    par(def.par)
    }
    else if(any(error.rate=="study"))
    {  
      
      def.par <- par(no.readonly = TRUE) 
      
      if (any(measure == "all"))
      measure = c("BER","overall")
    
    if (is.null(measure) || !any(measure %in% c("BER","overall")))
      stop("'measure' should be among the ones used in your call to 'perf': ", paste(c("BER","overall"),collapse = ", "),".")
    
      if (any(study == "all"))
        study = 1:length(x$study.specific.error)
      
      
      if (any(dist == "all"))
        dist = colnames(x$study.specific.error[[1]][[1]])
      
      
      if (is.null(dist) || !any(dist %in% colnames(x$study.specific.error[[1]][[1]])))
        stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$study.specific.error[[1]][[1]]),collapse = ", "),".")
      
      
      if(any(overlay=="all"))
        par(mfrow=c(1,length(study)))
      else if(any(overlay=="measure"))
        par(mfrow=c(length(study),length(dist)))
      
    for(stu in study)
    {
    # KA changed
    #x = matrix(x[, pred.method], ncol = length(pred.method))
    mat.error.plsda=matrix(nrow=nrow(x$study.specific.error[[stu]][[1]]),ncol=0)
    for(mea in measure)
    {
      mat.error.plsda=cbind(mat.error.plsda,x$study.specific.error[[stu]][[mea]][, dist])
    }
    
    if (is.null(ylab))
    {
      ylab = 'Classification error rate'
      
    }
    
    if (is.null(xlab))
    {
      xlab = 'PLSDA components'
      
    }
    
    
    internal_graph_plot.perf(mat.error.plsda,sd.error.plsda=NULL, overlay, type,measure,dist,legend,xlab,ylab,sd=FALSE,  ...)
    }
      
      par(def.par)
      
    }
    return(invisible(out))
    
  }

# SGCCDA object
# ----------------------

plot.perf.sgccda.mthd <-
  function (x,
            dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
            measure = "all",
            type="l",
            xlab = NULL,
            ylab = NULL,
            overlay=c("all","measure"),
            legend=c("vertical","horizontal"),
            ...)
  {
    
    if (hasArg(pred.method))
      stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
    pred.method = NULL # to pass R CMD check
    
    
    if (any(measure == "all"))
      measure = names(x$error.rate)
    
    if (is.null(measure) || !any(measure %in% names(x$error.rate)))
      stop("'measure' should be among the ones used in your call to 'perf': ", paste(names(x$error.rate),collapse = ", "),".")
    
    if (any(dist == "all"))
      dist = colnames(x$error.rate[[1]])
    
    
    if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
      stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
    
    # KA changed
    #x = matrix(x[, pred.method], ncol = length(pred.method))
    mat.error.plsda=matrix(nrow=nrow(x$error.rate[[1]]),ncol=0)
    for(mea in measure)
    {
      mat.error.plsda=cbind(mat.error.plsda,x$error.rate[[mea]][, dist])
    }
    
    if (is.null(ylab))
    {
      ylab = 'Classification error rate'
      
    }
    
    if (is.null(xlab))
    {
      xlab = 'PLSDA components'
      
    }
    def.par <- par(no.readonly = TRUE) 
    
    internal_graph_plot.perf(mat.error.plsda,sd.error.plsda=NULL, overlay, type,measure,dist,legend,xlab,ylab,sd=FALSE,  ...)
   
       par(def.par)
    return(invisible(out))
    
  }

