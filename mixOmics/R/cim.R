# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
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


# -------------------------------------
# CIM S3 method
# -------------------------------------
cim <-
  function(...) UseMethod("cim")


# -------------------------------------
# CIM default
# -------------------------------------

cim.default <- 
  function(mat,
           color = NULL,
           row.names = TRUE,
           col.names = TRUE,
           row.sideColors = NULL,
           col.sideColors = NULL,
           row.cex = NULL,
           col.cex = NULL,
           cluster = "both",
           dist.method = c("euclidean", "euclidean"),
           clust.method = c("complete", "complete"),
           cut.tree = c(0, 0),
           transpose = FALSE,
           symkey = TRUE, 
           keysize = c(1, 1),            
           zoom = FALSE, 
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           margins = c(5, 5),
           lhei = NULL,
           lwid = NULL,
           ...)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                   error = function(e) e)
    
    if ("simpleError" %in% class(err))
      stop(err[[1]], ".", call. = FALSE)
    
    function.arg = names(mget(names(formals()), sys.frame(sys.nframe())))
    not.arg = !(user.arg %in% function.arg)
    
    if (any(not.arg)) {
      unused.arg = user.arg[not.arg]
      not.arg = which(not.arg) + 1
      output = rep("", length(not.arg))
      
      for (i in 1:length(not.arg)) {
        output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
      }
      
      output = paste0("(", paste(output, collapse = ", "), ").")
      msg = "unused argument "
      if (length(not.arg) > 1) msg = "unused arguments "  
      stop(msg, output, call. = FALSE)
    }
    
    #-- mat
    isMat = tryCatch(is.matrix(mat), error = function(e) e)
    
    if ("simpleError" %in% class(isMat))
      stop(isMat[[1]], ".", call. = FALSE)
    
    if (!is.matrix(mat) || !is.numeric(mat))
      stop("'mat' must be a numeric matrix.", call. = FALSE)
    
    p = nrow(mat)
    q = ncol(mat)
    
    #-- color
    if (is.null(color)) 
      color = spectral.colors(25)
    
    #-- row.names
    if (is.logical(row.names)) {
      if(isTRUE(row.names)) row.names = rownames(mat) else row.names = rep("", p)
    }
    else {
      row.names = as.vector(row.names)
      if (length(row.names) != p)
        stop("'row.names' must be a character vector of length ", p, ".", call. = FALSE)
    }
    
    #-- col.names
    if (is.logical(col.names)) {
      if(isTRUE(col.names)) col.names = colnames(mat) else col.names = rep("", q)
    }
    else {
      col.names = as.vector(col.names)
      if (length(col.names) != q)
        stop("'col.names' must be a character vector of length ", q, ".", call. = FALSE)
    }
    
    #-- row.sideColors
     if (!is.null(row.sideColors)) {
       row.sideColors = as.matrix(row.sideColors)
       if (nrow(row.sideColors) != p)
         stop("'row.sideColors' must be a colors character vector (matrix) of length (nrow) ", p, ".", 
              call. = FALSE)
     }

    #-- col.sideColors
     if (!is.null(col.sideColors)) {
       col.sideColors = as.matrix(col.sideColors)
       if (nrow(col.sideColors) != q)
         stop("'col.sideColors' must be a colors character vector (matrix) of length (nrow) ", q, ".", 
              call. = FALSE)
     }
    
    #-- cluster
    choices = c("both", "row", "column", "none")
    cluster = choices[pmatch(cluster, choices)]
    
    if (is.na(cluster)) 
      stop("'cluster' should be one of 'both', 'row', 'column' or 'none'.", 
           call. = FALSE)
    
    #-- cluster method
    if (!is.character(clust.method) | length(as.vector(clust.method)) != 2)
      stop("'clust.method' must be a character vector of length 2.", call. = FALSE)
    
    choices = c("ward", "single", "complete", "average", "mcquitty",
                "median", "centroid")
    clust.method = choices[c(pmatch(clust.method[1], choices),
                             pmatch(clust.method[2], choices))]
    
    if (any(is.na(clust.method))) 
      stop("invalid clustering method.", call. = FALSE)
    
    #-- distance method
    if (!is.character(dist.method) | length(as.vector(dist.method)) != 2)
      stop("'dist.method' must be a character vector of length 2.", call. = FALSE)
    
    choices = c("euclidean", "correlation", "maximum", "manhattan", 
                "canberra", "binary", "minkowski")
    dist.method = choices[c(pmatch(dist.method[1], choices),
                            pmatch(dist.method[2], choices))]
    
    if (any(is.na(dist.method))) 
      stop("invalid distance method.", call. = FALSE)
    
    #-- checking other arguments
    check.cim.arg(color = color,
                  row.cex = row.cex,
                  col.cex = col.cex,
                  cut.tree = cut.tree,
                  transpose = transpose,
                  symkey = symkey, 
                  keysize = keysize,            
                  zoom = zoom, 
                  xlab = xlab, 
                  ylab = ylab,
                  main = main,
                  margins = margins,
                  row.sideColors = row.sideColors,
                  col.sideColors = col.sideColors,
                  lhei = lhei,
                  lwid = lwid)
    
    #-- end checking --#
    #------------------#
    
    
    #-- clustering -------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    if ((cluster == "both") || (cluster == "row")) {
      Rowv = rowMeans(mat)
      
      if (dist.method[1] == "correlation") 
        dist.mat = as.dist(1 - cor(t(as.matrix(mat)), method = "pearson"))
      else
        dist.mat = dist(mat, method = dist.method[1])
      
      hcr = hclust(dist.mat, method = clust.method[1])
      ddr = as.dendrogram(hcr)
      ddr = reorder(ddr, Rowv)
      rowInd = order.dendrogram(ddr)
      mat = mat[rowInd, ]
      row.names = row.names[rowInd]
      
      if (!is.null(row.sideColors)) 
        row.sideColors = as.matrix(row.sideColors[rowInd, ])
    }      
   
    if ((cluster == "both") || (cluster == "column")) {
      Colv = colMeans(mat)
      
      if (dist.method[2] == "correlation") 
        dist.mat = as.dist(1 - cor(as.matrix(mat), method = "pearson"))
      else
        dist.mat = dist(t(mat), method = dist.method[2])
      
      hcc = hclust(dist.mat, method = clust.method[2])
      ddc = as.dendrogram(hcc)
      ddc = reorder(ddc, Colv)
      colInd = order.dendrogram(ddc)
      mat = mat[, colInd]
      col.names = col.names[colInd]
      
      if (!is.null(col.sideColors)) 
        col.sideColors = as.matrix(col.sideColors[colInd, ])
    }
    
    #-- calling the image.map function -----------------------------------------#
    #---------------------------------------------------------------------------#
    imageMap(mat,
             color = color,
             row.names = row.names,
             col.names = col.names,
             row.sideColors = row.sideColors,
             col.sideColors = col.sideColors,             
             row.cex = row.cex,
             col.cex = col.cex,
             cluster = cluster,
             ddr = ddr,
             ddc = ddc,
             cut.tree = cut.tree,
             transpose = transpose,
             symkey = symkey, 
             keysize = keysize,            
             zoom = zoom, 
             main = main,
             xlab = xlab,
             ylab = ylab,
             margins = margins,
             lhei = lhei,
             lwid = lwid)          
    
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    res = list(mat = mat, row.names = row.names, col.names = col.names, 
               row.sideColors = row.sideColors, col.sideColors = col.sideColors)
    
    if ((cluster == "both") || (cluster == "row")) {
      res$rowInd = rowInd
      res$ddr = ddr    
    }
    
    if ((cluster == "both") || (cluster == "column")) {
      res$colInd = colInd
      res$ddc = ddc
    }
    
    class(res) = "cim_default"
    return(invisible(res))
  }


# -------------------------------------
# check CIM arguments
# -------------------------------------

check.cim.arg <-
  function(color = NULL,
           row.cex = NULL,
           col.cex = NULL,
           cut.tree = c(0, 0),
           transpose = FALSE,
           symkey = TRUE, 
           keysize = c(1, 1),            
           zoom = FALSE, 
           xlab = NULL, 
           ylab = NULL,
           main = NULL,
           margins = c(5, 5),
           row.sideColors = NULL,
           col.sideColors = NULL,
           lhei = NULL,
           lwid = NULL) 
{
    #-- checking general input arguments ---------------------------------------#
    #---------------------------------------------------------------------------#  
    
    #-- color
    if (any(!isColor(color))) 
      stop("'color' must be a character vector of recognized colors.", 
           call. = FALSE)
    
    #-- row.sideColors
    if (any(!isColor(row.sideColors))) 
      stop("color names for vertical side bar must be a character vector of recognized colors.", 
           call. = FALSE)
    
    #-- col.sideColors
    if (any(!isColor(col.sideColors))) 
      stop("color names for horizontal side bar must be a character vector of recognized colors.", 
           call. = FALSE)
    
    #-- row.cex
    if (!is.null(row.cex)) {
      if (!is.numeric(row.cex) || length(row.cex) != 1)
        stop("'row.cex' must be a numerical value.", call. = FALSE)
    }
    
    #-- col.cex
    if (!is.null(col.cex)) {
      if (!is.numeric(col.cex) || length(col.cex) != 1)
        stop("'col.cex' must be a numerical value.", call. = FALSE)
    }
      
    #-- transpose
    if (!is.logical(transpose))
      stop("'transpose' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- cut.tree
    if (!is.numeric(cut.tree) || length(cut.tree) != 2) 
      stop("'cut.tree' must be a numeric vector of length 2.",
           call. = FALSE)
    else {
      if (!(all(0 <= cut.tree & cut.tree <= 1)))
        stop("Components of 'cut.tree' must be between 0 and 1.",
             call. = FALSE) 
    }
    
    #-- keysize
    if (length(keysize) != 2 || any(!is.finite(keysize))) 
      stop("'keysize' must be a numeric vector of length 2.",
           call. = FALSE)
    
    #-- zoom
    if (!is.logical(zoom))
      stop("'zoom' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- margins
    if (!is.numeric(margins) || length(margins) != 2) 
      stop("'margins' must be a numeric vector of length 2.",
           call. = FALSE)
    
    #-- symkey
    if (!is.logical(symkey))
      stop("'symkey' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
   
    #-- lhei
    if (!is.null(lhei)) {
      if (is.null(col.sideColors)) {
        if (length(lhei) != 2 | !is.numeric(lhei) | any(is.na(lhei))) 
          stop("'lhei' must be a numeric vector of length 2.",
               call. = FALSE)
      }
      else {
        if (length(lhei) != 3 | !is.numeric(lhei) | any(is.na(lhei))) 
          stop("'lhei' must be a numeric vector of length 3.",
               call. = FALSE)
      }
    }
     
    #-- lwid
    if (!is.null(lwid)) {
      if (is.null(row.sideColors)) {
        if (length(lwid) != 2 | !is.numeric(lwid) | any(is.na(lwid))) 
          stop("'lwid' must be a numeric vector of length 2.",
               call. = FALSE)
      }
      else {
        if (length(lwid) != 3 | !is.numeric(lwid) | any(is.na(lwid))) 
          stop("'lwid' must be a numeric vector of length 3.",
               call. = FALSE)
      }
    }
    
    #-- xlab
    xlab = as.graphicsAnnot(xlab)
    
    #-- ylab
    ylab = as.graphicsAnnot(ylab)
    
    #-- main
    main = as.graphicsAnnot(main)
  }
