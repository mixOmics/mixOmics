# ------------------------
# date: 28/10/2015
# author: Florian Rohart
# run this script before releasing an update of the package; run all test scripts and output the errors
# ---------------------

# need to install the latest version to be tested
library(mixOmics)

scripts.files=list.files("running_scripts/")
scripts.files=scripts.files[grep(".R",scripts.files)]


out=lapply(1:length(scripts.files),function(x){try(suppressMessages(source(paste0("running_scripts/",scripts.files[x]))),silent=TRUE)})
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

cat("\n\nerror found in the following test files:\n\n")
for(i in ind.error)
{
    cat(scripts.files[i],"\n")
    
}
