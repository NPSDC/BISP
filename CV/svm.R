library(e1071)
library(foreach)
library(doParallel)
source('helper_func.R')
source('predict_models.R')
#' Gets the cross valiation stages for SVM for a given data
#' 
#' @param tr.data normalised training data containing the gene set
#' @param tr.model training model
#' @param folds no of folds for k-fold CV
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The predicted stages after 10 fold cross validation
cv.svm <- function(tr.data, tr.model, folds, stages.train, cores = 1, gamma = 0, kernel = 'linear', cost =1,
                   class.weights = sapply(levels(as.factor(stages.train)), function(x){
                     c(x=1)
                   } ))
{
  #sink("log.txt", append=TRUE)
  names(class.weights) = levels(as.factor(stages.train))
  total.samp <- nrow(tr.data)
  gr <- build.groups(total.samp, folds)
  
  if(length(gr) != folds)
    return('not correct')
  registerDoParallel(cores = cores)
  predicted.prob <- foreach(i = 1:folds, .combine = c) %dopar%
  {
    train.index = sort(unlist(gr[-i]))
    test.index = sort(unlist(gr[[i]]))
    svm.model <- svm(x = tr.data[train.index, ], y = stages.train[train.index], kernel = kernel, gamma = gamma)
    return(get.pr(tr.model = svm.model, test.data = tr.data[test.index,]))
  }
  predicted.prob[unlist(gr)] <- predicted.prob
  registerDoSEQ()
  return(list(pr.prob = predicted.prob, gr = gr, cla.pr = tr.model[['cla.pr']],
              cla.rem = tr.model[['cla.rem']]))
}

#' Gets the predicted stages for SVM for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param train.model.list List of training models
#' @param folds no of folds for k-fold CV
#' @param genes.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return List of predicted stages on CV SVM models for AFs
cv.svm.list <- function(tr.data, train.model.list, folds, genes.list, stages.train,
                        gamma = 0, kernel = 'linear', cores = 1, cost =1,
                        class.weights =if(length(levels(stages.train)) == 4) 
                          c('stage i' = 1, 'stage ii' =1, 'stage iii' = 1, 'stage iv' =1)
                        else c('stage i' = 1, 'stage iv' = 1))
{
  #writeLines(c(""), "log.txt")
  
  svm.pred <- mclapply(seq_along(genes.list), function(i)
  {
    if(is.null(genes.list[[i]]) | (length(genes.list[[i]]) == 1))
      return(NULL)
    else
     return(cv.svm(tr.data = tr.data[,genes.list[[i]]], folds = folds,
                   stages.train = stages.train, tr.model = train.model.list[[i]], cores = cores))
  }, mc.cores = cores)
  names(svm.pred) <- names(genes.list)
  return(svm.pred)
}