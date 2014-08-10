# Copyright (C) 2009 
# S?bastien D?jean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonz?lez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh L? Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
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


pls <-
function(X, 
         Y, 
         ncomp = 2, 
         mode = c("regression", "canonical", "invariant", "classic"),
         max.iter = 500, 
         tol = 1e-06,
         near.zero.var = TRUE,
         ...)
{

   result=spls(X,Y,ncomp,mode,max.iter,tol,near.zero.var,keepX=rep(ncol(X), ncomp),keepY=rep(ncol(Y), ncomp))
	
    class(result) = "pls"
    return(invisible(result))
}

