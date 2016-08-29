
internal_graph_plot.perf<- function (mat.error.plsda,sd.error.plsda, overlay, type,measure,dist,legend,xlab,ylab,sd,  ...)
{
if(any(overlay=="all"))
{
  out<-matplot(mat.error.plsda, type = type, lty = rep(c(1:length(measure)), each = length(dist)), col = rep(color.mixo(1:length(dist)), length(measure)), 
               lwd = 2, xlab = xlab, ylab = ylab,axes=FALSE,ylim=c(min(mat.error.plsda),max(mat.error.plsda)))
  axis(1,1:nrow(mat.error.plsda),rownames(mat.error.plsda))
  axis(2)
  if(any(legend=="vertical"))
  {legend('topright', legend = c(measure, dist), lty = c(1:length(measure), rep(NA, length(dist))), 
          pch = c(rep(NA, length(measure)), rep(16, length(dist))), col = c(rep('black',length(measure)), color.mixo(1:length(dist))), ncol = 1, lwd = 2)}
  else if(legend=="horizontal")
  {
    legend('topright', legend = c(measure,"" ,dist), lty = c(1:length(measure), rep(NA, (length(dist)+1))), 
           pch = c(rep(NA, (length(measure)+1)), rep(16, length(dist))), col = c(rep('black',length(measure)), NA, color.mixo(1:length(dist))), ncol = 2, lwd = 2)
    
  }
  if(sd)
 { for(col in 1:ncol(mat.error.plsda))
  {
    plot_error_bar(mat.error.plsda[,col],uiw=sd.error.plsda[,col],add=T,...)
  }}
}
else if(overlay=="measure")
{
  par(mfrow=c(1,length(dist)))
  for(di in dist)
  {
    new_mat.error=mat.error.plsda[,which(colnames(mat.error.plsda)==di)]
    out<-matplot(new_mat.error, type = type, lty = c(1:length(measure)), col ="black", 
                 lwd = 2, xlab = xlab, ylab = ylab,axes=FALSE,ylim=c(min(mat.error.plsda),max(mat.error.plsda)))
    axis(1,1:nrow(mat.error.plsda),rownames(mat.error.plsda))
    axis(2)
    if(any(legend=="vertical"))
    {legend('topright', legend = measure, lty = 1:length(measure),  col = rep('black',length(measure)), ncol = 1, lwd = 2)}
    
    else if(legend=="horizontal")
    {
      legend('topright', legend = c(measure), lty = 1:length(measure),  col = rep('black',length(measure)), ncol = 2, lwd = 2)
    }
    if(sd)
    {for(col in 1:ncol(new_mat.error))
    {
      new_sd.error=sd.error.plsda[,which(colnames(sd.error.plsda)==di)]
      plot_error_bar(new_mat.error[,col],uiw=new_sd.error[,col],add=T,...)
    }}
  }
  
}

}

plot_error_bar <- function (x, y = NULL, uiw,add=FALSE,  ...) 
{
  sfrac = 0.01
  gap = 0
  slty = par("lty")
  pt.bg = par("bg")
  arglist <- list(...)
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
  scols <- rep(scol, length.out = length(x))[nz]
  arrow.args <- c(list(lty = slty, angle = 90, length = smidge, 
                       code = 1, col = scols), clean.args(arglist, arrows, 
                                                          exclude.other = c("col", "lty", "axes")))
  do.call("arrows", c(list(x[nz], li[nz], x[nz], pmax(y - 
                                                        gap, li)[nz]), arrow.args))
  nz <- abs(ui - pmin(y + gap, ui)) * y.to.in > 0.001
  scols <- rep(scol, length.out = length(x))[nz]
  arrow.args$col <- scols
  do.call("arrows", c(list(x[nz], ui[nz], x[nz], pmin(y + 
                                                        gap, ui)[nz]), arrow.args))
  do.call("points", c(list(x, y, bg = pt.bg), clean.args(arglist, 
                                                         points, exclude.other = c("xlab", "ylab", "xlim", 
                                                                                   "ylim", "axes"))))
  invisible(list(x = x, y = y))
}
