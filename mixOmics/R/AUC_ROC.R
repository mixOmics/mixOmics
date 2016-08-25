#############################################################################################################
# Authors:
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Kim-Anh Le Cao, University of Queensland Diamantina Institute, Brisbane, Australia
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 23-08-2016
# last modified:  23-08-2016
#
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
#############################################################################################################

auc = function(object, ...)
UseMethod("auc")

auc.plsda = auc.splsda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
plot = TRUE,
ncomp = 1,
...)
{
    if(dim(newdata)[[1]]!=length(outcome.test))
    stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]], " elements.",call. = FALSE)
    
    data=list()
    statauc=list()
    data$outcome=as.factor(outcome.test)
    
    res.predict  =  predict(object, newdata = newdata, method = "max.dist")$predict
    
    for (i in 1:object$ncomp)
    {
        data$data=res.predict[,,i]
        title=paste("ROC Curve Comp",i)
        statauc[[paste("Comp", i, sep = "")]]=statauc(data, plot = ifelse(i%in%ncomp,plot,FALSE), title = title)
    }
    return(statauc)
}

auc.mint.plsda = auc.mint.splsda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
study.test = object$study,
plot = TRUE,
ncomp = 1,
...)
{

    if(dim(newdata)[[1]]!=length(outcome.test))
    stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]], " elements.",call. = FALSE)
    
    if(dim(newdata)[[1]]!=length(study.test))
    stop("Factor study.test must be a factor with ",dim(newdata)[[1]], " elements.",call. = FALSE)
    study.test=as.factor(study.test)
    
    data=list()
    statauc=list()
    data$outcome=as.factor(outcome.test)
    
    res.predict  =  predict(object, newdata = newdata, method = "max.dist",  study.test = study.test)$predict

    for (i in 1:object$ncomp)
    {
        data$data=res.predict[,,i]
        title=paste("ROC Curve Comp",i)
        statauc[[paste("Comp", i, sep = "")]]=statauc(data, plot = ifelse(i%in%ncomp,plot,FALSE), title = title)
    }
    return(statauc)
}

auc.sgccda = function(
object,
newdata = object$X,
outcome.test = as.factor(object$Y),
plot = TRUE,
block = 1,
ncomp = 1,
...)
{
    
    data=list()
    auc.mean=list()
    data$outcome=as.factor(outcome.test)

    res.predict  =  predict(object, newdata = newdata, method = "max.dist")$predict
    block.all = names(res.predict)
    block.temp = names(res.predict[block])
    
    for(j in 1:length(res.predict))
    {
        for (i in 1:object$ncomp[j])
        {
            data$data=res.predict[[j]][,,i]
            title=paste("ROC Curve\nBlock: ", names(res.predict)[j], ", comp: ",i, sep="")
            
            plot.temp = ifelse(i%in%ncomp && names(res.predict)[j]%in%block.temp, plot, FALSE)
            auc.mean[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = statauc(data, plot = plot.temp, title = title)
            
        }
    }
    return(auc.mean)
}

