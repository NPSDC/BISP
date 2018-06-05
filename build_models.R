library(pamr)
library(randomForest)
library(e1071)
library(foreach)
library(doParallel)

#' Gets the trained Shrunken Centroid classifiers for the AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param genes.list List of AFs of a particular feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.shrunken.classifier <- function(tr.data, genes.list, stages.train, cores = 1)
{
  cla.pr <- levels(as.factor(stages.train))[1]
  cla.rem <- levels(as.factor(stages.train))[2]
  shrunken.train.list <- mclapply(genes.list, function(genes)
  {
    mod <- pamr.train(data = list(x = as.matrix(t(tr.data[,genes])), y = stages.train))
    mod$cla.pr <- cla.pr
    mod$cla.rem <- cla.rem
    mod
  }, mc.cores = cores)
  names(shrunken.train.list) <- names(genes.list)
  return(shrunken.train.list)
}

#' Returns the trained Random Forest classifiers for AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param genes.list List of AFs of a feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.rf.classifier <- function(tr.data, genes.list, stages.train, cores = 1)
{
  cla.pr <- levels(as.factor(stages.train))[1]
  cla.rem <- levels(as.factor(stages.train))[2]
  rf.train.list <- mclapply(genes.list, function(genes)
  {
    rf.mod <- randomForest(tr.data[,genes], stages.train)
    rf.mod$cla.pr <- cla.pr
    rf.mod$cla.rem <- cla.rem
    rf.mod
    
  }, mc.cores = cores)
  names(rf.train.list) = names(genes.list)
  return(rf.train.list)
}

#' Returns the trained SVM classifiers for AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param genes.list List of AFs of a feature selection algorithm
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.svm.classifier <- function(tr.data, genes.list, stages.train, cores = 1,
                                 gamma = 0, kernel = 'linear', cost =1,
                                 class.weights =if(length(levels(stages.levels)) == 4) 
                                   c('stage i' = 1, 'stage ii' =1, 'stage iii' = 1, 'stage iv' =1)
                                 else c('stage i' = 1, 'stage iv' = 1))
{
  svm.train.list <- mclapply(genes.list, function(genes)
  {
    svm.model <- svm(x = tr.data[, genes], y = stages.train, kernel = kernel, gamma = gamma, probability = T)
    if(svm.model$decision.values[1,1] > 0  & stages.train[1] !=
             levels(as.factor(stages.train))[1])
      svm.model$decision.values[,1] <- svm.model$decision.values[,1]*-1
      
      
    svm.model$cla.pr <- levels(as.factor(stages.train))[1]
    svm.model$cla.rem <- levels(as.factor(stages.train))[2]
    svm.model
  }, mc.cores = cores)
  names(svm.train.list) <- names(genes.list)
  return(svm.train.list)
}

#' Returns the trained naive bayes classifiers for AFs of a feature selection algorithm
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param genes.list List of AFs of a feature selection algorithm
#' @param cores Number of CPU cores to be used
#' @return The list containing the trained models for the AFs
build.nb.classifier <- function(tr.data, genes.list, stages.train, cores = 1)
{
  cla.pr <- levels(as.factor(stages.train))[1]
  cla.rem <- levels(as.factor(stages.train))[2]
  nb.train.list <- mclapply(seq_along(genes.list), function(i){
  nb.mod <- naiveBayes(tr.data[, genes.list[[i]]], stages.train)
  nb.mod$cla.pr <- cla.pr
  nb.mod$cla.rem <- cla.rem
  nb.mod
  }, mc.cores = cores)
  names(nb.train.list) <- names(genes.list)
  return(nb.train.list)
}

#' Returns the different training classifiers for the AFs of different feature selection algorithms
#' 
#' @param tr.data normalised training data containing samples as rows and columns as genes
#' @param fea.list List of list of AFs yielded by different feature selection algorithms 
#' @param stages.train Stage of every sample in training data
#' @param cores Number of CPU cores to be used
#' @return The list containing the different trained models for every feature in fea.list
get.train.model <- function(tr.data, fea.list, stages.train, cores = 1)
{
  registerDoParallel(cores = cores)
  train.model <- foreach(i = 1:length(fea.list)) %dopar%
  {
    fea.name <- names(fea.list)[i]
    model <- list()
    model[['shrunken']] <- build.shrunken.classifier(tr.data = tr.data, genes.list = fea.list[[fea.name]], 
                                                     stages.train = stages.train, cores = cores)
    model[['rf']] <- build.rf.classifier(tr.data = tr.data, genes.list = fea.list[[fea.name]], 
                                         stages.train = stages.train, cores = cores)
    model[['nb']] <- build.nb.classifier(tr.data = tr.data, genes.list = fea.list[[fea.name]], 
                                         stages.train = stages.train, cores = cores)
    model[['svm']] <- build.svm.classifier(tr.data = tr.data, genes.list = fea.list[[fea.name]], 
                                           stages.train = stages.train, cores = cores)
    return(model)
  }
  names(train.model) = names(fea.list)
  return(train.model)
}
