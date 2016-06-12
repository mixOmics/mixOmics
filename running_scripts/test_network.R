#######################################################################################################
#######################################################################################################
#                                           from help file
#######################################################################################################
#######################################################################################################
opar <- par(no.readonly = TRUE)

## network representation for objects of class 'rcc'
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

network(nutri.res, comp = 1:3, threshold = 0.6)


## Changing the attributes of the network

network(nutri.res, comp = 1:3, threshold = 0.7,
color.node = c("mistyrose", "lightcyan"),
shape.node = c("circle", "rectangle"),
color.edge = color.jet(100),
lty.edge = "solid", lwd.edge = 2,
show.edge.labels = FALSE)


## interactive 'threshold'
#network(nutri.res, comp = 1:3, threshold = 0.55, interactive = TRUE)
## select the 'threshold' and "see" the new network


## network representation for objects of class 'spls'
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
keepY = c(10, 10, 10))

network(toxicity.spls, comp = 1:3, threshold = 0.8,
color.node = c("mistyrose", "lightcyan"),
shape.node = c("rectangle", "circle"),
color.edge = color.spectral(100),
lty.edge = "solid", lwd.edge =  1,
show.edge.labels = FALSE, interactive = FALSE)



#######################################################################################################
#######################################################################################################
#                                           additional tests
#######################################################################################################
#######################################################################################################

if(additional.test==TRUE)
{
    
    data(liver.toxicity)
    
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    toxicity.spls <- spls(X, Y, ncomp = 3,
    keepX = c(50, 50, 50), keepY = c(10, 10, 10))
    
    network(toxicity.spls)
    
    # test comp
    network(toxicity.spls,comp=1)
    network(toxicity.spls,comp=c(1,2))
    network(toxicity.spls,comp=c(3,2))
    network(toxicity.spls,comp=1:3)
    
    # test threshold
    network(toxicity.spls, threshold = 0.1)
    network(toxicity.spls, threshold = 0.8)
    network(toxicity.spls, threshold = 0.2)
    
    # test row.names, col.names
    network(toxicity.spls,row.names = FALSE)
    network(toxicity.spls,col.names = FALSE)
    network(toxicity.spls,row.names = FALSE, col.names = FALSE)
    
    # test color.node
    network(toxicity.spls,color.node = c("orange","blue"))
    network(toxicity.spls,color.node = c("orange","blue"),row.names = FALSE, col.names = FALSE)
    
    # test shape.node
    network(toxicity.spls,shape.node = c('none','circle'))
    network(toxicity.spls,shape.node = c('rectangle','none'))
    network(toxicity.spls,shape.node = c('circle','rectangle'))
    
    # test shape.node + color.node
    network(toxicity.spls,shape.node = c('circle','rectangle'),color.node = c("orange","blue"))
    network(toxicity.spls,shape.node = c('rectangle','none'),color.node = c("orange","blue"))
    network(toxicity.spls,shape.node = c('none','circle'),color.node = c("orange","blue"))
    
    # test cex.node.name
    network(toxicity.spls,cex.node.name = 0.5)
    network(toxicity.spls,cex.node.name = 1.5)
    ##(rajouter length=2)
    
    # test color.edge
    network(toxicity.spls,color.edge = c("blue","red"))
    network(toxicity.spls,color.edge = c("red","orange","cyan","blue"),row.names=FALSE)
    
    # test lty.edge, lwd.edge
    network(toxicity.spls,lty.edge='dashed')
    network(toxicity.spls,lty.edge='dotted',row.names = FALSE)
    network(toxicity.spls,lwd.edge=3,row.names = FALSE)
    network(toxicity.spls,lwd.edge=2,lty.edge = "twodash",row.names = FALSE)
    
    # test show.edge.label
    network(toxicity.spls,show.edge.labels = TRUE,row.names=FALSE)
    # test cex.edge.label
    network(toxicity.spls,show.edge.labels = TRUE,cex.edge.label = 2,row.names=FALSE)
    network(toxicity.spls,show.edge.labels = FALSE,cex.edge.label = 2,row.names=FALSE)
    
    # test show.color.key
    network(toxicity.spls,show.color.key = FALSE)
    
    # test symkey
    network(toxicity.spls,symkey = FALSE)
    
    # test symkey + color.edge
    network(toxicity.spls,symkey = FALSE, color.edge=c("red","blue","orange"),row.names=FALSE) # ??Color Key
    
    #test keysize
    network(toxicity.spls,keysize=c(1,2))
    network(toxicity.spls,keysize=c(2,1),symkey = FALSE)
    network(toxicity.spls,keysize=c(2,1),show.color.key=FALSE)
    
    # test interactive
    #network(toxicity.spls,interactive = TRUE)
    #network(toxicity.spls,threshold = 0.3, interactive = TRUE)
    #network(toxicity.spls,threshold = 0.6, interactive = TRUE)
    
    
    ####### RCC
    
    data(linnerud)
    
    X <- linnerud$exercise
    Y <- linnerud$physiological
    linn.res <- rcc(X, Y,ncomp=3)
    
    network(linn.res)
    
    # test comp
    network(linn.res,comp=1)
    network(linn.res,comp=c(1,2))
    network(linn.res,comp=c(3,2)) ####
    network(linn.res,comp=1:3)
    
    # test threshold
    network(linn.res, threshold = 0)
    network(linn.res, threshold = 0.8)
    network(linn.res, threshold = 0.4)
    
    # test row.names, col.names
    network(linn.res,row.names = FALSE)
    network(linn.res,col.names = FALSE)
    network(linn.res,row.names = FALSE, col.names = FALSE)
    
    # test color.node
    network(linn.res,color.node = c("orange","blue"))
    network(linn.res,color.node = c("orange","blue"),row.names = FALSE, col.names = FALSE)
    
    # test shape.node
    network(linn.res,shape.node = c('none','circle'))
    network(linn.res,shape.node = c('rectangle','none'))
    network(linn.res,shape.node = c('circle','rectangle'))
    
    # test shape.node + color.node
    network(linn.res,shape.node = c('circle','rectangle'),color.node = c("orange","blue"))
    network(linn.res,shape.node = c('rectangle','none'),color.node = c("orange","blue"))
    network(linn.res,shape.node = c('none','circle'),color.node = c("orange","blue"))
    
    # test cex.node.name
    network(linn.res,cex.node.name = 0.5)
    network(linn.res,cex.node.name = 1.5)
    
    # test color.edge
    network(linn.res,color.edge = c("blue","red"))
    network(linn.res,color.edge = c("red","orange","cyan","blue"))
    
    # test lty.edge, lwd.edge
    network(linn.res,lty.edge='dashed')
    network(linn.res,lty.edge='dotted',row.names = FALSE)
    network(linn.res,lwd.edge=3,row.names = FALSE)
    network(linn.res,lwd.edge=2,lty.edge = "twodash",row.names = FALSE)
    
    # test show.edge.label
    network(linn.res,show.edge.labels = TRUE,row.names=FALSE)
    # test cex.edge.label
    network(linn.res,show.edge.labels = TRUE,cex.edge.label = 2,row.names=FALSE)
    network(linn.res,show.edge.labels = FALSE,cex.edge.label = 2,row.names=FALSE)
    
    # test show.color.key
    network(linn.res,show.color.key = FALSE)
    
    # test symkey
    network(linn.res,symkey = FALSE)
    
    # test symkey + color.edge
    network(linn.res,symkey = FALSE, color.edge=c("red","blue","orange"),row.names=FALSE)
    
    #test keysize
    network(linn.res,keysize=c(1,2))
    network(linn.res,keysize=c(2,1),symkey = FALSE)
    network(linn.res,keysize=c(2,1),show.color.key=FALSE)
    
    # test interactive
    #network(linn.res,interactive = TRUE)
    #network(linn.res,threshold = 0.3, interactive = TRUE)
    #network(linn.res,threshold = 0.6, interactive = TRUE)
    
    
    ####### MULTILEVEL (SPLSDA)
    
    data(liver.toxicity)
    repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
    6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
    10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
    13, 14, 15, 16, 15, 16, 15, 16, 15, 16)
    
    design <- data.frame(sample = repeat.indiv)
    res.spls.1level <- spls(X = liver.toxicity$gene,
    Y=liver.toxicity$clinic,
    multilevel = design,
    ncomp = 3,
    keepX = c(50, 50, 50), keepY = c(5, 5, 5),
    mode = 'canonical')
    
    
    network(res.spls.1level)
    
    # test comp
    network(res.spls.1level,comp=1)
    network(res.spls.1level,comp=c(1,2))
    network(res.spls.1level,comp=c(3,2))
    network(res.spls.1level,comp=1:3)
    
    # test threshold
    network(res.spls.1level, threshold = 0)
    network(res.spls.1level, threshold = 0.8)
    network(res.spls.1level, threshold = 0.7)
    
    # test row.names, col.names
    network(res.spls.1level,row.names = FALSE)
    network(res.spls.1level,col.names = FALSE)
    network(res.spls.1level,row.names = FALSE, col.names = FALSE)
    
    # test color.node
    network(res.spls.1level,color.node = c("orange","blue"))
    network(res.spls.1level,color.node = c("orange","blue"),row.names = FALSE, col.names = FALSE)
    
    # test shape.node
    network(res.spls.1level,shape.node = c('none','circle'))
    network(res.spls.1level,shape.node = c('rectangle','none'))
    network(res.spls.1level,shape.node = c('circle','rectangle'))
    
    # test shape.node + color.node
    network(res.spls.1level,shape.node = c('circle','rectangle'),color.node = c("orange","blue"))
    network(res.spls.1level,shape.node = c('rectangle','none'),color.node = c("orange","blue"))
    network(res.spls.1level,shape.node = c('none','circle'),color.node = c("orange","blue"))
    
    # test cex.node.name
    network(res.spls.1level,cex.node.name = 0.5)
    network(res.spls.1level,cex.node.name = 1.5)
    
    # test color.edge
    network(res.spls.1level,color.edge = c("blue","red"))
    network(res.spls.1level,color.edge = c("red","orange","cyan","blue"))
    
    # test lty.edge, lwd.edge
    network(res.spls.1level,lty.edge='dashed')
    network(res.spls.1level,lty.edge='dotted',row.names = FALSE)
    network(res.spls.1level,lwd.edge=3,row.names = FALSE)
    network(res.spls.1level,lwd.edge=2,lty.edge = "twodash",row.names = FALSE)
    
    # test show.edge.label
    network(res.spls.1level,show.edge.labels = TRUE,row.names=FALSE)
    # test cex.edge.label
    network(res.spls.1level,show.edge.labels = TRUE,cex.edge.label = 2,row.names=FALSE)
    network(res.spls.1level,show.edge.labels = FALSE,cex.edge.label = 2,row.names=FALSE)
    
    # test show.color.key
    network(res.spls.1level,show.color.key = FALSE)
    
    # test symkey
    network(res.spls.1level,symkey = FALSE)
    
    # test symkey + color.edge
    network(res.spls.1level,symkey = FALSE, color.edge=c("red","orange","blue"),row.names=FALSE)
    
    #test keysize
    network(res.spls.1level,keysize=c(1,2))
    network(res.spls.1level,keysize=c(2,1),symkey = FALSE)
    network(res.spls.1level,keysize=c(2,1),show.color.key=FALSE)
    
    
    # test interactive
    #network(res.spls.1level,interactive = TRUE)
    #network(res.spls.1level,threshold = 0.3, interactive = TRUE)
    #network(res.spls.1level,threshold = 0.6, interactive = TRUE)
    
    
    ####### BLOCKS (SGCCA)
    
    ### 2 blocks
    data(nutrimouse)
    
    diet = unmap(nutrimouse$diet)
    blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE, dimnames = list(names(blocks), names(blocks)))
    
    nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = 4)
    
    network(nutri.sgcca)
    
    # test comp
    network(nutri.sgcca,comp=list(gene=1,lipid=1))
    network(nutri.sgcca,comp=list(gene=2,lipid=2))
    network(nutri.sgcca,comp=list(gene=c(1,2),lipid=c(2,3)))
    network(nutri.sgcca,comp=list(lipid=c(1,2),gene=c(2,3)))
    network(nutri.sgcca,comp=list(lipid=c(2,3),gene=c(1,2)))
    
    
    # test threshold
    network(nutri.sgcca, threshold = 0)
    network(nutri.sgcca, threshold = 0.8)
    network(nutri.sgcca, threshold = 0.7)
    
    # test block.var.names
    network(nutri.sgcca,block.var.names = FALSE)
    network(nutri.sgcca,block.var.names = c(FALSE,TRUE))
    
    
    # test color.node
    network(nutri.sgcca,color.node = c("orange","blue"))
    
    # test shape.node
    network(nutri.sgcca,shape.node = c(gene='none',lipid='circle'))
    network(nutri.sgcca,shape.node = c(gene='rectangle',lipid='none'))
    network(nutri.sgcca,shape.node = c(gene='circle',lipid='rectangle'))
    
    # test shape.node + color.node
    network(nutri.sgcca,shape.node = c('circle','rectangle'),color.node = c(lipid="orange",gene="blue")) #couleur invers?e
    network(nutri.sgcca,shape.node = c('rectangle','none'),color.node = c("orange","blue"))
    network(nutri.sgcca,shape.node = c('none','circle'),color.node = c("orange","blue"))
    
    # test cex.node.name
    network(nutri.sgcca,cex.node.name = 0.5)
    network(nutri.sgcca,cex.node.name = 1.5)
    
    if(FALSE)# we can only display 2 blocks
    {
        ### 3 blocks
        data(nutrimouse)
        
        diet = unmap(nutrimouse$diet)
        blocks = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, diet = diet)
        design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE, dimnames = list(names(blocks), names(blocks)))
        
        nutri.sgcca <- wrapper.sgcca(blocks,design=design, ncomp = 4)
        
        network(nutri.sgcca)
        
        # test comp
        network(nutri.sgcca,comp=list(gene=1,lipid=1,diet=c(1,2)),blocks=1:3)
        network(nutri.sgcca,comp=list(gene=c(1,2),lipid=c(2,3),diet=1:3),blocks=1:3)
        
        
        # test threshold
        network(nutri.sgcca, threshold = 0,blocks=1:3)
        network(nutri.sgcca, threshold = 0.8,blocks=1:3)
        network(nutri.sgcca, threshold = 0.7,blocks=1:3)
        
        # test block.var.names
        network(nutri.sgcca,block.var.names = FALSE,blocks=1:3)
        network(nutri.sgcca,block.var.names = c(FALSE,FALSE,TRUE),blocks=1:3)
        network(nutri.sgcca,block.var.names = list(rep("",120),rep("",21),c("d1","d2","d3","d4","d5")),blocks=1:3)
        
        # test color.node
        network(nutri.sgcca,color.node = c("orange","blue","cyan"),blocks=1:3)
        
        # test shape.node
        network(nutri.sgcca,shape.node = c(gene='none',lipid='circle',diet="rectangle"),blocks=1:3)
        network(nutri.sgcca,shape.node = c(gene='rectangle',lipid='none',diet="circle"),blocks=1:3)
    }
    
}

par(opar)
