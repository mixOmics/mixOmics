# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created:  18-02-2016
# last modified 19-02-2016
# Copyright (C) 2016
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


which.function.should.I.use=function()
{
    function.one.study.one.TypeOfData=c("pls","plsda","spls","splsda")
    function.one.study.multiple.TypeOfData=paste0("block.",function.one.study.one.TypeOfData) # block.pls
    function.multiple.study.one.TypeOfData=paste0("mint.",function.one.study.one.TypeOfData) # mint.pls
    function.multiple.study.multiple.TypeOfData=paste0("mint.",function.one.study.multiple.TypeOfData) # mint.block.pls
    
    
    all.functions=c(function.one.study.one.TypeOfData,function.one.study.multiple.TypeOfData,
    function.multiple.study.one.TypeOfData,function.multiple.study.multiple.TypeOfData)
    
    
    #------------
    # Question 1
    #------------
    message("Can I ask you a few questions to identify your needs?")
    Q1=readline("yes or no. \n ")
    while(is.na(charmatch(Q1,c("yes","no"))) | charmatch(Q1,c("yes","no"))==0)
    {
        Q1=readline("Please answer `yes' or `no'\n ")
    }
    
    if(!is.na(charmatch(Q1,"no")))
    {
        message("Understood. You can come back here anytime you want.")
    }else{
        #------------
        # Question 2: block analysis
        #------------
        message("Q2: Do you have multiple types of data that you want to analyse together, e.g. transcriptomics, proteomics, mintbolomics? If you are unsure, look up the example.")
        Q2=readline("yes, no or example. \n ")
        
        while(is.na(charmatch(Q2,c("yes","no","example"))) | charmatch(Q2,c("yes","no","example"))==0)
        {
            Q2=readline("Please answer `yes', `no' or `example' '\n ")
        }
        
        if(!is.na(charmatch(Q2,"example")))# Q2==example, plot an example
        {
            message("Example has been plotted, hope it helps.")
            #to run example, find the file to plot, import it and plot it!
            fpath <- system.file("extdata", "jellyfish.jpg", package="mixOmicsv6")
            library(jpeg)
            img=readJPEG(fpath)
            opar=par(no.readonly=TRUE)
            par(mar=c(1,1,3,1))
            plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "",xaxt="n",yaxt="n",main="Horizontal analysis",bty="n")
            rasterImage(img,100,300,250,450)
            par(opar)

            message("Q2: Do you have multiple types of data that you want to analyse together, e.g. transcriptomics, proteomics, mintbolomics? ")
            Q2=readline("yes or no. \n ")
            while(is.na(charmatch(Q2,c("yes","no"))) | charmatch(Q2,c("yes","no"))==0)
            {
                Q2=readline("Please answer `yes', or `no'' '\n ")
            }
        }
        
        if(!is.na(charmatch(Q2,"no")))# Q2==no, no block
        {
            final.function=c(function.one.study.one.TypeOfData,function.multiple.study.one.TypeOfData)#initialise the final answer
            cat("You are analysing a single type of data. You can use the following functions",final.function,"\n",sep="\n") # pls, mint.pls
        }else# Q2==yes, block
        {
            final.function=c(function.one.study.multiple.TypeOfData,function.multiple.study.multiple.TypeOfData)#initialise the final answer
            cat("You are doing a horizontal analysis. You can use the following functions\n",final.function,"\n",sep="\n") # mint.block.pls, block.pls,
        }
        
        
        
        
        #------------
        # Question 3: mint analysis
        #------------
        message("Q3: Do you have multiple experiments of the same type of data that you want to combine together, e.g. transcriptomics from several studies? If you are unsure, look up the example.")
        Q3=readline("yes, no or example. \n ")
        
        while(is.na(charmatch(Q3,c("yes","no","example"))) | charmatch(Q3,c("yes","no","example"))==0)
        {
            Q3=readline("Please answer `yes', `no' or `example' '\n ")
        }
        
        if(!is.na(charmatch(Q3,"example")))# Q2==example, plot an example
        {
            if(!is.na(charmatch(Q2,"no")))# no block
            {
                message("Example has been plotted, hope it helps.")
                fpath <- system.file("extdata", "cat.jpeg", package="mixOmicsv6")# mint without block
                library(jpeg)
                img=readJPEG(fpath)
                opar=par(no.readonly=TRUE)
                par(mar=c(1,1,3,1))
                plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "",xaxt="n",yaxt="n",main="Vertical analysis",bty="n")
                rasterImage(img,100,300,250,450)
                par(opar)

            }else{# block
                message("Example has been plotted, hope it helps.")
                fpath <- system.file("extdata", "dog.jpeg", package="mixOmicsv6")# mint and block
                library(jpeg)
                img=readJPEG(fpath)
                opar=par(no.readonly=TRUE)
                par(mar=c(1,1,3,1))
                plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "",xaxt="n",yaxt="n",main="Horizontal and vertical analysis",bty="n")
                rasterImage(img,100,300,250,450)
                par(opar)

            }
            message("Q3: Do you have multiple experiments of the same type of data that you want to combine together, e.g. transcriptomics from several studies? ")
            Q3=readline("yes or no. \n ")
            
            while(is.na(charmatch(Q3,c("yes","no"))) | charmatch(Q3,c("yes","no"))==0)
            {
                Q3=readline("Please answer `yes', or `no' '\n ")
            }
        }
        
        if(!is.na(charmatch(Q3,"no"))) # no mint
        {
            if(!is.na(charmatch(Q2,"no")))
            {
                final.function=function.one.study.one.TypeOfData
                cat("You are analysing a single type of data and a single study. You can use the following functions",final.function,"\n",sep="\n") # pls
            }else{# Q2==multiple
                final.function=function.one.study.multiple.TypeOfData
                cat("You are doing a horizontal analysis without combining different studies of the same type of data. You can use the following functions:",final.function,"\n",sep="\n")# block.pls
            }
        }else# Q2==yes, mint
        {
            if(!is.na(charmatch(Q2,"one")))
            {
                final.function=function.multiple.study.one.TypeOfData
                cat("You are analysing a single type of data and several studies. You can use the following functions",final.function,"\n",sep="\n") # mint.pls
            }else{# Q2==multiple
                final.function=function.multiple.study.multiple.TypeOfData
                cat("You are doing a horizontal and vertical analysis and combining different studies of the same type of data. You can use the following functions:",final.function,"\n",sep="\n")# mint.block.pls
            }
        }
        
        
        #------------
        # Question 4: Discriminant Analysis
        #------------
        message("Q4: Do you want to discrimate your samples? Is your outcome a categorical variable? If you are unsure, look up the example. ")
        Q4=readline("yes, no or example. \n ")
        
        
        while(is.na(charmatch(Q4,c("yes","no","example"))) | charmatch(Q4,c("yes","no","example"))==0)
        {
            Q4=readline("Please answer `yes', `no' or `example' '\n ")
        }
        
        if(!is.na(charmatch(Q4,"example")))# Q2==example, plot an example
        {
            message("Example has been plotted, hope it helps.")
            fpath <- system.file("extdata", "DA-analysis.jpg", package="mixOmicsv6")# mint without block
            library(jpeg)
            img=readJPEG(fpath)
            opar=par(no.readonly=TRUE)
            par(mar=c(1,1,3,1))
            plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "",xaxt="n",yaxt="n",main="Discriminant Analysis",bty="n")
            rasterImage(img,100,300,250,450)
            par(opar)

        }
        
        message("Q4: Do you want to discrimate your samples? Is your outcome a categorical variable? ")
        Q4=readline("yes or no. \n ")
        
        while(is.na(charmatch(Q4,c("yes","no"))) | charmatch(Q4,c("yes","no"))==0)
        {
            Q4=readline("Please answer `yes', or `no' '\n ")
        }
        
        if(!is.na(charmatch(Q4,"yes")))
        {
            final.function=final.function[grep("da",final.function)]
            cat("You want to perform a Discriminant Analysis. You can use the following functions",final.function,"\n",sep="\n") # mint.pls
            
            
        }else{
            final.function=final.function[-grep("da",final.function)]
            cat("You are not performing a Discriminant Analysis. You can use the following functions",final.function,"\n",sep="\n") # mint.pls
            
        }
        
        
        #------------
        # Question 5: sparse Analysis
        #------------
        message("Q5: Do you want to identify relevant features, as biomarkers/genes?")
        Q5=readline("yes or no. \n ")
        while(is.na(charmatch(Q5,c("yes","no"))) | charmatch(Q5,c("yes","no"))==0)
        {
            Q5=readline("Please answer `yes' or `no'\n")
        }
        
        if(!is.na(charmatch(Q5,"yes")))
        {
            final.function=final.function[grep("sp",final.function)]
            cat("You want to perform a sparse Analysis. You can use the following function",final.function,"\n",sep="\n") # mint.pls
            
            
        }else{
            final.function=final.function[-grep("sp",final.function)]
            cat("You are not performing a sparse Analysis. You can use the following function",final.function,"\n",sep="\n") # mint.pls
            
        }
    }
    
}








