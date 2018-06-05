library(randomForest)
library(foreach)
library(doParallel)

#' Gets the cross valiation stages for Random Forest for a given data
#' 
#' @param tr.data normalised training data containing the gene set
#' @param tr.model training model
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages after 10 fold cross validation
cv.rf <- function(tr.data, tr.model, folds, stages.train, cores = 1, sampsize = if (replace) nrow(data) else ceiling(.632*nrow(data)))
{
  #writeLines(c(""), "log.txt")
  #sink('log.txt', append = T)
  total.samp <- nrow(tr.data)
  gr <- build.groups(total.samp, folds)
  registerDoParallel(cores = cores)
  predicted.prob <- foreach(i = 1:folds, .combine = c) %dopar%
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[[i]]))
    rf.model <- randomForest(x = tr.data[train.index, ], y = stages.train[train.index])
    
    return(get.pr(tr.model = rf.model, test.data = tr.data[test.index,]))
  }
  predicted.prob[unlist(gr)] <- predicted.prob
  registerDoSEQ()
  return(list(pr.prob = predicted.prob, gr = gr, cla.pr = tr.model[['cla.pr']], 
              cla.rem = tr.model[['cla.rem']]))
}

#' Gets the cross valiation models for RFs for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param train.model.list List of trained models w.r.t feature selection algorithm
#' @param folds no of folds for k-fold CV
#' @param genes.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages on CV SVM models for AFs
cv.rf.list <- function(tr.data, train.model.list, folds, genes.list, stages.train, cores = 1,
                       sampsize = if (replace) nrow(data)
                       else ceiling(.632*nrow(data)))
{
  rf.pred <- mclapply(seq_along(genes.list), function(i)
  {
    if(is.null(genes.list[[i]]) | (length(genes.list[[i]]) == 1))
      return(NULL)
    else
    return(cv.rf(tr.data = tr.data[, genes.list[[i]]], tr.model = train.model.list[[i]],
                 folds =  folds, 
                 stages.train = stages.train, cores = cores))
  }, mc.cores = cores)
  names(rf.pred) <- names(genes.list)
  return(rf.pred)
}
