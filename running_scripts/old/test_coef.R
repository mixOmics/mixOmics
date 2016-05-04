data(linnerud)
X <- linnerud$exercise
Y <- linnerud$physiological

#test ncomp
linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")
coef.pls=coef(linn.pls)
coef.pls$comp1
linn.pls <- pls(X, Y, ncomp = 3, mode = "classic")
coef.pls=coef(linn.pls)
coef.pls$comp1
coef.pls

#test comp

linn.pls <- pls(X, Y, ncomp = 3)

coef.pls=coef(linn.pls, comp = c(1,3))
coef.pls
coef.pls=coef(linn.pls, comp = 3)
coef.pls


X <- pentosys$gene
Y <- pentosys$metabolom
pento.spls <- spls(X, Y, ncomp = 3,
                   keepX = c(50, 50, 50), keepY = c(10, 10, 10))
coef.spls=coef(pento.spls)
coef.spls

#test keep.coef
coef.spls=coef(pento.spls,comp=1,keep.coef = FALSE)
length(which(coef.spls$comp1==0))
coef.spls=coef(pento.spls,comp=1,keep.coef = TRUE)
length(which(coef.spls$comp1==0)) #nb différent de 0
coef.spls=coef(pento.spls,comp=2,keep.coef = TRUE)
length(which(coef.spls$comp2==0)) #nb différent de 0


data(vac18.simulated)
X <- vac18.simulated$genes
design <- data.frame(samp = vac18.simulated$sample,
                     stim = vac18.simulated$stimulation,
                     time = vac18.simulated$time)

res.2level <- multilevel(X, ncomp = 2, design = design,
                         keepX = c(120, 10), method = 'splsda')
coef.ml=coef(res.2level,keep.coef = FALSE)
length(which(coef.spls$comp1==0))
coef.ml=coef(res.2level,comp=c(2,1))
coef.ml$comp2
coef.ml=coef(res.2level,comp=c(1,2))
coef.ml$comp2
