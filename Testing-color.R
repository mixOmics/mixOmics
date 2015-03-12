# -----------------------------------------------------------------------------------
# Testing-color.R
# Author:    KA Le Cao 
# Date started:  26/02/2015
# Last updated:  12/03/2015
# Objective: testing an update of the new mixOmics palette
# Latest update: 
# -----------------------------------------------------------------------------------

library(mixOmics)

# -----------------------
# existing jet colors
# ----------------------
par(mfrow = c(3, 1))
z <- seq(-1, 1, length = 125)
for (n in c(11, 33, 125)) {
  image(matrix(z, ncol = 1), col = jet.colors(n), 
        xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
  box()
  par(usr = c(-1, 1, -1, 1))  	
  axis(1, at = c(-1, 0, 1))
}

# 
# # --------------------------------
# # including the mixOmics colors
# # -------------------------------
# mixo.gray = gray.colors(1, start = 0.76, gamma = 1)
# 
# mixo.col = c('#388ECC', # mixOmics logo blue
#              '#F68B33', # mixOmics logo orange
#              mixo.gray, # mixOmics logo grey
#              '#009E73', # shiny dark green
#              '#CC79A7', # shiny purple/pink
#              '#F0E442' #shiny yellow
#              #'#D55E00', #shiny dark orange
#              #'#0072B2', #shiny dark blue
#              #'#999999',  # shiny grey
#              #'#E69F00', # shiny orange
#              #'#56B4E9' #Shiny blue
#              )
# 
# par(mfrow=c(1,1))            
# plot( 1:length(mixo.col), pch = 18, col = mixo.col, cex = 5)            


# -------------------
# testing the function mixo.colors
# -------------------
source('color.palettes.R')
help(pca)
data(multidrug)

pca.res <- pca(multidrug$ABC.trans, ncomp = 4, scale = TRUE)

my.colors = mixo.colors(as.numeric(factor(multidrug$cell.line$Class)))
# samples representation
plotIndiv(pca.res, ind.names = multidrug$cell.line$Class, cex = 0.5, 
          col = my.colors)

mixo.colors(3)
