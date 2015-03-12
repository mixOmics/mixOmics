# Copyright (C) 2015 
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
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


#-- jet colors --#
#----------------#
jet.colors <-
  function (n, alpha = 1) 
  {
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- n
    if (length(n) > 1 || !is.finite(n))
      stop("'n' must be an integer positive value.", 
           call. = FALSE)
    
    if (n < 1)
      stop("'n' must be an integer positive value.", 
           call. = FALSE)
    
    #-- alpha
    if (length(alpha) > 1 || !is.finite(alpha))
      stop("'alpha' must be an numeric value in the range [0, 1].", 
           call. = FALSE)
    
    if (alpha < 0 || alpha > 1)
      stop("'alpha' must be an numeric value in the range [0, 1].", 
           call. = FALSE)
    
    alpha = round(255 * alpha)
    
    #-- end checking --#
    #------------------#
    
    ramp = colorRamp(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
                       "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
                       "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
                       "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
                       "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
                       "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
                       "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
                       "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
                       "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
                       "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
                       "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
                       "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
                       "#AF0000", "#9F0000", "#8F0000", "#800000"), 
                     space = "Lab")
    
    rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
  }


#-- spectral colors --#
#---------------------#
spectral.colors <-
  function (n, alpha = 1) 
  {
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- n
    if (length(n) > 1 || !is.finite(n))
      stop("'n' must be an integer positive value.", 
           call. = FALSE)
    
    if (n < 1)
      stop("'n' must be an integer positive value.", 
           call. = FALSE)
    
    #-- alpha
    if (length(alpha) > 1 || !is.finite(alpha))
      stop("'alpha' must be an numeric value in the range [0, 1].", 
           call. = FALSE)
    
    if (alpha < 0 || alpha > 1)
      stop("'alpha' must be an numeric value in the range [0, 1].", 
           call. = FALSE)
    
    alpha = round(255 * alpha)
    
    #-- end checking --#
    #------------------#
    
    ramp = colorRamp(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", 
                       "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", 
                       "#9E0142"), space = "Lab")
    
    rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
  }


#-- green-black-red gradient colors --#
#-------------------------------------#
GreenRed.colors <-
  function (n, alpha = 1) 
  {
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- n
    if (length(n) > 1 || !is.finite(n))
      stop("'n' must be an integer positive value.", 
           call. = FALSE)
    
    if (n < 1)
      stop("'n' must be an integer positive value.", 
           call. = FALSE)
    
    #-- alpha
    if (length(alpha) > 1 || !is.finite(alpha))
      stop("'alpha' must be an numeric value in the range [0, 1].", 
           call. = FALSE)
    
    if (alpha < 0 || alpha > 1)
      stop("'alpha' must be an numeric value in the range [0, 1].", 
           call. = FALSE)
    
    alpha = round(255 * alpha)
    
    #-- end checking --#
    #------------------#
    
    ramp = colorRampPalette(c("green", "darkgreen", "black", "darkred", "red"))
    ramp = ramp(101)
    green = ramp[1:43]
    red = ramp[59:101]
    ramp = colorRamp(c(green, "black", red), space = "Lab")
    
    rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
  }





#-- mixOmics colors --#
#---------------------#
mixo.colors <- function(num.vector){
  
  if(is.factor(num.vector)) num.vector=as.numeric(num.vector)
  
  if(!is.numeric(num.vector)) stop(paste("num.vector has to be numeric", call. = FALSE))
  
  # these are the colors in the logo (the first 3) and used in the Shiny web interface
  mixo.gray = gray.colors(1, start = 0.76, gamma = 1)
  
  mixo.col = c('#388ECC', # mixOmics logo blue
               '#F68B33', # mixOmics logo orange
               mixo.gray, # mixOmics logo grey
               '#009E73', # shiny dark green
               '#CC79A7', # shiny purple/pink
               '#F0E442', #shiny yellow
               'black',
               '#D55E00', #shiny dark orange
               '#0072B2', #shiny dark blue
               '#999999'  # shiny grey
               #'#E69F00', # shiny orange
               #'#56B4E9' #Shiny blue
  )
  
  #-- checking general input parameters --------------------------------------#
  #---------------------------------------------------------------------------#
  
  n = length(num.vector)
  #-- n: check that there are more colors available than requested
  if (isTRUE(num.vector) > length(mixo.col)){
    stop(paste("We only have a few mix.colors available, n <= ", length(mixo.col)), call. = FALSE)
  }
  if (isTRUE(!is.finite((num.vector))) ||  (n < 1)){
    stop("'num.vector' must be an integer vector with positive values.", call. = FALSE)
  }
  #-- end checking --#
  #------------------#
  
  return(mixo.col[num.vector])
}


