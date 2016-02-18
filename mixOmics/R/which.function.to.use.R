# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created:  18-02-2016
# last modified 18-02-2016
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
    function.multiple.study.one.TypeOfData=paste0("meta.",function.one.study.one.TypeOfData) # meta.pls
    function.multiple.study.multiple.TypeOfData=paste0("meta.",function.one.study.multiple.TypeOfData) # meta.block.pls
    
    
    all.functions=c(function.one.study.one.TypeOfData,function.one.study.multiple.TypeOfData,
    function.multiple.study.one.TypeOfData,function.multiple.study.multiple.TypeOfData)
    
    #------------
    # Question 1
    #------------
    
    Q1=readline("We are going to ask you a few questions to identify your needs, ok?  yes or no. \n ")
    while(is.na(charmatch(Q1,c("yes","no"))))
    {
        Q1=readline("Please answer `yes' or `no'\n ")
    }
    
    if(!is.na(charmatch(Q1,"no")))
        stop("Understood. You can come back here anytime you want.")
    
    #------------
    # Question 2: block analysis
    #------------
    Q2=readline("Q2: Do you have multiple types of data that you want to analyse together, e.g. transcriptomics, proteomics, metabolomics? \n one or multiple. \n ")
    while(is.na(charmatch(Q2,c("one","multiple"))))
    {
        Q2=readline("Please answer `one' or `multiple'\n ")
    }
    
    if(!is.na(charmatch(Q2,"one")))
    {
        final.function=c(function.one.study.one.TypeOfData,function.multiple.study.one.TypeOfData)#initialise the final answer
        cat("You are analysing a single type of data. At this point, you can use the following functions",final.function,"\n",sep="\n") # pls, meta.pls
    }else{# Q2==multiple
        final.function=c(function.one.study.multiple.TypeOfData,function.multiple.study.multiple.TypeOfData)#initialise the final answer
        cat("You are doing a horizontal analysis. At this point, you can use the following functions\n",final.function,"\n",sep="\n") # meta.block.pls, block.pls,
    }
    
    
    #------------
    # Question 3: meta analysis
    #------------
    Q3=readline("Q3: Do you have multiple experiments of the same type of data that you want to combine together, e.g. transcriptomics from several studies? \n one or multiple. \n ")
    while(is.na(charmatch(Q2,c("one","multiple"))))
    {
        Q3=readline("Please answer `one' or `multiple'\n ")
    }
    
    if(!is.na(charmatch(Q3,"one"))) # no meta
    {
        if(!is.na(charmatch(Q2,"one")))
        {
            final.function=function.one.study.one.TypeOfData
            cat("You are analysing a single type of data and a single study. At this point, you can use the following functions",final.function,"\n",sep="\n") # pls
        }else{# Q2==multiple
            final.function=function.one.study.multiple.TypeOfData
            cat("You are doing a horizontal analysis without combining different studies of the same type of data. At this point, you can use the following functions:",final.function,"\n",sep="\n")# block.pls

        }

        
    }else{# Q3==multiple, meta
        if(!is.na(charmatch(Q2,"one")))
        {
            final.function=function.multiple.study.one.TypeOfData
            cat("You are analysing a single type of data and several studies. At this point, you can use the following functions",final.function,"\n",sep="\n") # meta.pls
        }else{# Q2==multiple
            final.function=function.multiple.study.multiple.TypeOfData
            cat("You are doing a horizontal analysis and combining different studies of the same type of data. At this point, you can use the following functions:",final.function,"\n",sep="\n")# meta.block.pls
        }

        
    }
    
    
    #------------
    # Question 4: Discriminant Analysis
    #------------
    Q4=readline("Q4: Is your outcome Y a categorical variable, i.e. do you want to classify your samples? \n yes or no. \n ")
    while(is.na(charmatch(Q4,c("yes","no"))))
    {
        Q4=readline("Please answer `yes' or `no'\n")
    }
    
    if(!is.na(charmatch(Q4,"yes")))
    {
        final.function=final.function[grep("da",final.function)]
        cat("You want to perform a Discriminant Analysis. At this point, you can use the following functions",final.function,"\n",sep="\n") # meta.pls


    }else{
        final.function=final.function[-grep("da",final.function)]
        cat("You are not performing a Discriminant Analysis. At this point, you can use the following functions",final.function,"\n",sep="\n") # meta.pls

    }
    
    
    #------------
    # Question 5: sparse Analysis
    #------------
    Q5=readline("Q5: Do you want to identify relevant features, as biomarkers/genes? \n yes or no. \n ")
    while(is.na(charmatch(Q5,c("yes","no"))))
    {
        Q5=readline("Please answer `yes' or `no'\n")
    }
    
    if(!is.na(charmatch(Q5,"yes")))
    {
        final.function=final.function[grep("sp",final.function)]
        cat("You want to perform a sparse Analysis. At this point, you can use the following function",final.function,"\n",sep="\n") # meta.pls
        
        
    }else{
        final.function=final.function[-grep("sp",final.function)]
        cat("You are not performing a sparse Analysis. At this point, you can use the following function",final.function,"\n",sep="\n") # meta.pls
        
    }
    
    
    
}








