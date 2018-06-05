library(varSelRF)
library(foreach)
library(doParallel)
#' Gets the varself object w.r.t each group
#' 
#' @param data normalised data containing samples as rows and columns as genes
#' @param train.ind.list List containing indexes of the folds
#' @param stages Stage of every sample in data
#' @param cores Number of CPU cores to be used
#' @param type Indicates whether there are 1 or multiple(type = 1) folds in train.ind.list
#' @return List containing varselrf objects for all groups and AFs
get.varselRF.object <- function(data, train.ind.list, stages, cores = 1, type = 1)
{
  #writeLines(c(""), "log.txt")
  registerDoParallel(cores)
  #sink("log.txt", append = T)
  models.vars <- foreach(i=1:length(train.ind.list)) %dopar% 
  {
    if(type == 1)
    {
      train.indexes = sort(unlist(train.ind.list[-i]))
      #print(train.indexes)
    }
    else
      train.indexes = sort(unlist(train.ind.list[[i]]))
    model <- varSelRF(xdata = as.matrix(data[train.indexes,]), 
                       Class = as.factor(stages[train.indexes]))
    
    return(model)
  }
  return(models.vars)
}

#' Gets the genes having minimum OOB error 
#' 
#' @param varselRf.ob.list List containing varselrf object w.r.t each group
#' @return List containing genes w.r.t each group
get.min.oob.varselRf <- function(varselRf.ob.list)
{
  sel.genes.list <- sapply(varselRf.ob.list, function(obj)
  {
    obj = obj$selec.history
    min.oob.error <- min(obj[,3])
    best.pos <-
      which(obj[,3] == min.oob.error)[which.min(obj[,1][which(obj[,3] == min.oob.error)])]
    genes <- obj[best.pos, 2]
    genes <- as.character(genes)
    genes <- strsplit(genes, split = ' + ', fixed = T)
  })
  return(sel.genes.list)
}
