# ----------------------------------------
# After compile, changing the workshop material
# date: 15/03/2015
# latest update: 15/03/2015.
# ----------------------------------------

# ================================================================
# Testing for workshop material (starting from version 07/10/2014)
# ================================================================
detach("package:mixOmics", unload=TRUE)

library(mixOmics, lib.loc = 'MyR/')
sessionInfo() # ok mixOmics 5.0-4



###################################################
### code chunk number 141: cs_multilevel.Rnw:33-39
###################################################
data(vac18)
?vac18
X <- vac18$genes
trt <- vac18$stimulation
# summarises the treatments
summary(trt)
# sumamrises the number of repeated measurements per individual
summary(as.factor(vac18$sample))


###################################################
### code chunk number 142: cs_multilevel.Rnw:48-52
###################################################

# setup the design matrix by indicating the repeated measurements
design <- data.frame(sample = vac18$sample)
# calculate the within matrix Xw
Xw <- withinVariation(X = X, design = design)
dim(Xw)
pca.multilevel.vac18 <- pca(Xw, ncomp = 3, scale = TRUE, center = TRUE)
pca.multilevel.vac18


###################################################
### code chunk number 143: cs_multilevel.Rnw:57-59
###################################################
# compare with a 'normal' PCA (analysis flawed since PCA assumes that all samples are independent!)

pca.classic.vac18 <-pca(X, ncomp = 3, scale = TRUE, center = TRUE)
pca.classic.vac18


###################################################
### code chunk number 144: pca-vac18
###################################################

# set up colors for plotIndiv
col.stimu <- c("darkblue", "purple", "green4","red3")
col.stimu <- col.stimu[as.numeric(trt)]

plotIndiv(pca.classic.vac18, ind.names = trt, col = col.stimu)
title(main = 'VAC18, classical PCA, comp 1 - 2')


###################################################
### code chunk number 145: pca-multilevel-vac18
###################################################
plotIndiv(pca.multilevel.vac18, ind.names = trt, col = col.stimu)
title(main = 'VAC18, PCA multilevel, comp 1 - 2')


###################################################
### code chunk number 146: cs_multilevel.Rnw:93-99
###################################################
# classical PLS-DA
vac18.plsda <-plsda(X, Y=trt, ncomp = 3)

# multilevel PLS-DA
# we redefine the design matrix. 
#Since PLS-DA is supervised, we need to mention the outcome of interest
design <- data.frame(sample = vac18$sample,
                     stimu = trt)
# note that here we dont need to specify Y as the information is contained in 'design'
vac18.plsda.multilevel <- multilevel(X, design = design, ncomp = 3, 
                                     method = "splsda",
                                     keepX = c(ncol(X), ncol(X), ncol(X)))


###################################################
### code chunk number 147: plsda-vac18
###################################################
# numbers here indicate the ID of the individuals
plotIndiv(vac18.plsda, col = col.stimu, ind.names = vac18$sample)
title('VAC18, classical PLS-DA, comp 1 - 2')


###################################################
### code chunk number 148: plsda-multilevel-vac18
###################################################
# numbers here indicate the ID of the individuals
plotIndiv(vac18.plsda.multilevel, col = col.stimu, ind.names = vac18$sample)
title('VAC18, multilevel PLS-DA, comp 1 - 2')


###################################################
### code chunk number 149: cs_multilevel.Rnw:123-127
###################################################
vac18.splsda.multilevel <- multilevel(X, design = design, ncomp = 3, 
                                    method = "splsda", 
                                    keepX = c(30, 137, 123))


###################################################
### code chunk number 150: multilevel-plotIndiv3d (eval = FALSE)
###################################################
plot3dIndiv(vac18.splsda.multilevel , ind.names = trt, col = col.stimu, 
            axes.box = "both")


###################################################
### code chunk number 151: heatmap-vac18
###################################################

# set up colors for pheatmap
col.samp <- c("lightgreen", "red", "lightblue", "darkorange",
              "purple", "maroon", "blue", "chocolate", "turquoise",
              "tomato1", "pink2", "aquamarine")
# define the colors for the conditions for the heatmap
col.stimu.htm <- c("darkblue", "purple", "green4","red3")
col.stimu.htm <- col.stimu.htm[as.numeric(trt)]
# ! make sure the number of colors matches the number of levels in Y
col.stimu.htm <- unique(col.stimu.htm)

pheatmap.multilevel(vac18.splsda.multilevel, 
                    # colors:
                    col_sample = col.samp, 
                    col_stimulation = col.stimu.htm,
                    #labels:
                    label_annotation = c("Subject", "Stimulus"),
                    # scaling:
                    scale = 'row',
                    # distances and clutering
                    clustering_distance_rows = "euclidean", 
                    clustering_distance_cols = "euclidean", 
                    clustering_method = "complete",
                    #  show col/row names and font
                    show_colnames = FALSE,
                    show_rownames = FALSE)


###################################################
### code chunk number 152: cs_multilevel.Rnw:169-175
###################################################
data(liver.toxicity)
# note: we made up those data, pretending they are repeated measurements
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each


###################################################
### code chunk number 153: cs_multilevel.Rnw:180-186
###################################################

# here we want to integrate genes and clinical measurements with a sPLS.
# the design matrix only include the repeated measurement information since spls is unsupervised
# (no information about the treatement will be taken into account)

design <- data.frame(sample = repeat.indiv)
vac18.spls.multilevel <- multilevel(X = liver.toxicity$gene,
                              Y=liver.toxicity$clinic,
                              design = design,
                              ncomp = 3,
                              keepX = c(50, 50, 50), keepY = c(5, 5, 5),
                              method = 'spls', mode = 'canonical')


###################################################
### code chunk number 154: cs_multilevel.Rnw:190-207 (eval = FALSE)
###################################################
# variable plots
plotVar(vac18.spls.multilevel, comp = 1:2, X.label = TRUE, Y.label = TRUE, 
        cex = c(0.5, 0.9)) 

jpeg('graphics/CIM_spls_multilevel.jpeg')
CIM <- cim(vac18.spls.multilevel, comp = 1:2, xlab = "genes", 
           ylab = "clinic var", 
           margins = c(5, 6), zoom = FALSE)
dev.off()

jpeg('graphics/network_spls_multilevel.jpeg')
network(vac18.spls.multilevel, comp = 1:2, threshold = 0.8, 
        Y.names = NULL, keep.var = TRUE,
        color.node = c( "lightcyan","mistyrose"),
        shape.node = c("circle", "rectangle"),
        color.edge = jet.colors(100),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
        show.edge.labels = FALSE, interactive = FALSE)
dev.off()


###################################################
### code chunk number 155: cs_multilevel.Rnw:215-239 (eval = FALSE)
###################################################
data(vac18.simulated)  # simulated data from VAC18 with 2 levels

X <- vac18.simulated$genes
design <- data.frame(sample = vac18.simulated$sample,
                     stimu = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

vac18.splsda2.multilevel <- multilevel(X, ncomp = 2, design = design,
                                       keepX = c(200, 200), method = 'splsda')



# set up color for plotIndiv
col.stimu <- as.numeric(vac18.simulated$stimulation)
# pch for plots
pch.time <- c(20, 4)
pch.time[as.numeric(vac18.simulated$time)] 

plotIndiv(vac18.splsda2.multilevel, col = col.stimu, pch = pch.time, ind.names = FALSE)
# combined legend!
legend('bottomright', col = c(unique(col.stimu), 'black', 'black'), 
       legend = c(levels(design$stimu), levels(design$time)),
       lty = c(1,1,1,1, NA, NA), lwd = 2,
       pch = c(NA, NA, NA, NA, unique(pch.time)), cex = 0.8)

# set up colors for each level of pheatmap 
col.sample <- c("lightgreen", "red","lightblue","darkorange","purple","maroon") # 6 samples
col.time <- c("pink","lightblue1") # two time points
# here the coloc vector should be the length of the # of stimulations
col.stimu.htm <- c('green', 'black', 'red', 'blue') 
# set up labels for the 2 levels in design matrix
label.stimu <- unique(design[, 2])
label.time <- unique(design$time)



#jpeg('pheatmap_splsda_multilevel2.jpeg')
pheatmap.multilevel(vac18.splsda2.multilevel,
                    # colors:
                    col_sample=col.sample, 
                    col_stimulation=col.stimu.htm, 
                    col_time=col.time,
                    #labels for each level
                    label_color_stimulation=label.stimu,
                    label_color_time=label.time, 
                    #clustering method
                    clustering_method="ward",
                    #show col/row names and font size
                    show_colnames = FALSE,
                    show_rownames = TRUE,
                    fontsize_row=2)
#dev.off()



