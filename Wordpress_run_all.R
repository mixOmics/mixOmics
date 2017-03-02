# created on 25-01-2017
# last modified: 25-01-2017
# Author: Florian Rohart
# purpose: knit/sweave all the rmd file from the mixOmics.org website

#this file needs to be sourced in R

scripts.files=list.files("wordpress_mixomics/rmd/")
scripts.files=scripts.files[grep(".Rmd",scripts.files)]#[5]
scripts.files

out=lapply(1:length(scripts.files),function(x){print(scripts.files[x]);try(suppressMessages(rmarkdown::render(paste0("wordpress_mixomics/rmd/",scripts.files[x]), output_dir = "wordpress_mixomics/html/")),silent=TRUE)})
a=lapply(out,class)
ind.error=grep("try-error",a)


for(i in 1:length(scripts.files))
{
    cat("====================\n")
    print(scripts.files[i])
    if(i%in%ind.error)
    {
        cat(out[[i]],"\n")
    }else{
        cat("all good\n\n")
    }
    cat("====================\n")
    
}

if(length(ind.error)>0)
{
    cat("\n\nerror found in the following wordpress files:\n\n")
    for(i in ind.error)
    {
        cat(scripts.files[i],"\n")
        
    }
    cat("\n\n the best way forward is to rmarkdown::render each faulty script and traceback() to get the line that broke the code\n")
}else{
    
    cat("\n\n NO ERROR found! \n\n")
    
}