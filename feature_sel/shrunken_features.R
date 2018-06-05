library(pamr)
library(pROC)
library(mccr)
library(foreach)
library(doParallel)
source('feature_sel/pamr.listgenes.R')
source('get_results.R')
source('helper_func.R')
#' Extracts genes from the shrunken object
#' 
#' @param shrunken.genes.df.list list containing data frames of output pamr.listgenes() w.r.t each group
#' @return List containing gene names for each group
get.genes.shrunken <- function(shrunken.genes.df.list)
{
  genes = list()
  for(i in seq_along(shrunken.genes.df.list))
  {
    genes[[i]] = shrunken.genes.df.list[[i]][['genes.list']][,2]
  }
  return(genes)
}

#' Gives the shrunken gene object for given data w.r.t each group
#' 
#' @param data normalised data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds
#' @return List containing gene object for all groups and AFs
get.shrunken.object <- function(data, train.ind.list, stages, cores, type = 1, min.genes = 1, min.range = 0.02)
{
  registerDoParallel(cores = cores)
  pamr.genes.list <- foreach(i = 1:length(train.ind.list)) %dopar%
  {
    if(type == 1)
      train.ind <- sort(unlist(train.ind.list[-i]))
    else
      train.ind <- sort(unlist(train.ind.list[[i]]))
    train.model <- pamr.train(list(x = as.matrix(t(data[train.ind,])), 
                                            y = stages[train.ind]))
    cv.model <- pamr.cv(train.model, 
                                 data = list(x = t(as.matrix(data[train.ind,])),
                                             y = stages[train.ind]),nfold = 10)
    type <- as.factor(stages.levels.comb[train.ind])[1]
    mccs <- sapply(seq_along(cv.model$threshold), function(x)
    {
      mccr(get.order(stages[train.ind], type), get.order(cv.model$yhat[,x], type))
    })
    thr.ind = sort(which(mccs == max(mccs)), decreasing = T)[1]
    
    if(min.genes != 1)
      thr.ind <- max(which(mccs > (mccs[thr.ind] - min.range)))  
    
    
    genes.list <- pamr.listgene(train.model, 
                                          data = list(x=as.matrix(t(data[train.ind,])),
                                                      y=stages[train.ind]), 
                                          threshold = cv.model$threshold[thr.ind],
                                          fitcv = cv.model, genenames = T)
    return(list(genes.list = genes.list, aucs = aucs))
  }
  return(pamr.genes.list)
}
