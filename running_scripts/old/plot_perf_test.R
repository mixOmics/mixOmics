data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

liver.pls <- pls(X, Y, ncomp = 3)
liver.perf <- perf(liver.pls, validation = "Mfold")

plot(liver.perf, criterion = "R2", type = "l", layout = c(2, 2))


#test criterion
plot(liver.perf,criterion = "R2")
plot(liver.perf,criterion = "Q2")
plot(liver.perf,criterion = "Q2.total")
plot(liver.perf,criterion = "MSEP") 



#test xlab,ylab
plot(liver.perf,criterion = "R2",xlab="Xlab")
plot(liver.perf,criterion = "Q2",ylab="Y")
plot(liver.perf,criterion = "MSEP",xlab="xlab",ylab="ylab")

#test cTicks
plot(liver.perf,criterion = "R2",cTicks=1:100)
plot(liver.perf,criterion = "R2",cTicks=5:1)
plot(liver.perf,criterion = "R2",1:5)

#test limQ2
plot(liver.perf,criterion = "Q2",limQ2=0.01)
plot(liver.perf,criterion = "Q2",limQ2=1)
plot(liver.perf,criterion = "Q2",limQ2=5)

#test limQ2.col
plot(liver.perf,criterion = "Q2",limQ2=1,limQ2.col=7)
plot(liver.perf,criterion = "Q2",limQ2=1,limQ2.col="green")
plot(liver.perf,criterion = "Q2",limQ2.col="green")

#test layout
plot(liver.perf,criterion = "MSEP",layout=c(2,2),style="lattice")
plot(liver.perf,criterion = "MSEP",layout=c(4,2),style="lattice")


#test type
plot(liver.perf,criterion = "Q2",type="p")
plot(liver.perf,criterion = "Q2",type="o")
plot(liver.perf,criterion = "Q2",type="b")
plot(liver.perf,criterion = "Q2",type="p",style="lattice")
plot(liver.perf,criterion = "Q2",type="o",style="lattice")
plot(liver.perf,criterion = "Q2",type="b",style="lattice")