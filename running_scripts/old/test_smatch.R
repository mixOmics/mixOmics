data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, mode = "canonical", ncomp = 3, 
                      keepX = c(50, 50, 50), keepY = c(10, 10, 10))

color.toxicity <- as.numeric(liver.toxicity$treatment[, 2])
label.toxicity <- liver.toxicity$treatment[, 1]
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)], 
        cex.lab = 0.5, label = label.toxicity, col = color.toxicity)


#color
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        col="red")
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        col=rep(c("red","blue"),each=32))

#pch
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        pch=1)
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        pch=5)

#cex
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        cex=1.2)
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        cex=2)

#label
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        label=TRUE)

#cex.lab
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        label=TRUE,cex.lab=0.8)
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        label=FALSE,cex.lab=0.8)

#arrows
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        arrows=FALSE)

#xlim, ylim
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        xlim=c(0,5))
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        ylim=c(0,5))
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        xlim=c(0,5),ylim=c(-5,5))

#add.axes
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        add.axes=FALSE)

#origin
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        origin = c(5,5))
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        origin = c(20,20))

#include.origin
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        origin = c(50,50),include.origin=F)

#main
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        main="plot")
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        main="plot plot plot plot plot plot")

#cex.main
s.match(toxicity.spls$variates$X[, c(1, 2)], 
        toxicity.spls$variates$Y[, c(1, 2)],
        main="plot",cex.main=0.5)



