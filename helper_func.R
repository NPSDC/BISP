#' Creates random partitions based on total number required and total samples
#' 
#' @param total.samples The total number of samples
#' @param  num.group The number of partitions that need to be created
#' @return a list containing the indexes within with each partition
build.groups <- function(total.samples, num.group)
{
  gr=NULL
  num.sel=floor(total.samples/num.group)
  num.add=total.samples%%num.group
  range=1:total.samples
  rr=total.samples
  for(i in 1:num.group)
  {
    vrem=sample(1:rr,size=num.sel)
    sel=range[vrem]
    gr=c(gr,list(sel))
    range=range[-vrem]
    rr=rr-num.sel
  }
  
  if(num.add>0)
  {
    vrem=sample(1:num.group,num.add)
    for(i in 1:num.add)
    {
      gr[[vrem[i]]]=c(gr[[vrem[i]]],range[i])
    }
  }
  gr <- lapply(gr, sort)
  return(gr)
}

#' Gives the stage distribution of a partition/s
#' @param  gr List containing the partitions
#' @param stages Vector containing the overall stages
#' @return A list containing the stage distribution of each sample within the partition 
get.stage.distribution <- function(gr, stages)
{
  stage.dist <- lapply(gr, function(x)
  {
    table(stages[x])
  })
  return(stage.dist)
}

#' Gives the genes present in atleast max no of models from the given genes list
#' 
#' @param genes.list List containing genes across the different groups
#' @param max.no.of.models No of groups should the genes atleast be present in
#' @return final genes that are atleast present in the asked number of groups
get.genes.common <- function(genes.list, max.no.of.models)
{
  total.genes <- Reduce(union, genes.list)
  genes.req <- c()
  for(i in total.genes)
  {
    count = 0
    for(j in seq_along(genes.list))
    {
      if(i %in%  genes.list[[j]])
        count = count + 1
      if(count == max.no.of.models)
      {
        genes.req <- c( genes.req, i)
        break
      }
    }
  }
  return(genes.req)
}

#'Creates the deseq2 res object
#'
#' @param dds_comp deseq object that would be used for differential analysis
#' @param contrast column in colData of dds object that would be used for comparison
#' @param stage.cont  base stage that would be used for comparision
#' @param stages.tum vector of stages that would be compared  against stage.cont
#' @return The list of resut matrix for comparision of stages.tum against stage.cont
comp.res <- function(dds.comp, contrast, stage.cont, stages.tum)
{
  stages.comp.cont <- lapply(stages.tum, function(x)
  {
    res = results(dds.comp, contrast = c(contrast, x, stage.cont))
    indexes = is.na(res[,5]) | is.na(res[,6])
    res = res[!indexes,]
  })
  names(stages.comp.cont) = stages.tum
  return(stages.comp.cont)
}

#' Gets the stages in 1 or 0 form with 1 being the dominant one
#' 
#' @param stages Vector that needs to be converted in the above form
#' @param type Level of stages that would be the dominant category
#' @return Vector with stages converted into 1 or 0
get.order <- function(stages, type)
{
  return(ifelse(stages == type, 1, 0))
}
#' Gets the genes from result object passing the threshold of log2FC and pvalue
#' 
#' @param res The result object created by deseq2
#' @param logfc log2FC threshold
#' @param adj.pval adjusted pvalue cutoff
#' @return the genes passing the threshold criteria
#' @param pval pvalue cutoff  
get.genes <- function(res, logfc, adj.pval, pval)
{
  return(rownames(res)[abs(res[,2]) > logfc & res[,6] < adj.pval & res[,5] < pval])
}

#' Gets the Capitalised full Classifier Name
#' 
#' @param cla classifier name to be converted
#' @param return updated classifier name
get.class.name <- function(cla)
{
  if(cla == 'knn')
    return(c('KNN'))
  else if(cla == 'nb')
    return(c('NB'))
  else if(cla == 'rf')
    return(c('RF'))
  else if(cla == 'shrunken')
    return(c('SC'))
  else
    return(c('SVM'))
}

#' Creates a consolidated data frame containing the value of an evaluation metric w.r.t all AFs for every feature selection w.r.t every classifier
#' 
#' @param test.ac List returned by get.cv.res or get.test.res
#' @return data.frame that can be used by ggplot for drawing plots w.r.t each metric
create.net.df <- function(test.ac)
{
  req.df <- list()
  req.df$roc_auc <- c()
  req.df$accuracy <- c()
  req.df$sens <- c()
  req.df$spec <- c()
  req.df$mcc <- c()
  req.df$pr_auc
  req.df$f_val <- c()
  req.df$classifier <- c()
  req.df$feature_sel <- c()
  for(j in seq_along(test.ac))
  {
    test.ob <- test.ac[[j]]
    for(i in seq_along(test.ob))
    {
      req.df$roc_auc <- c(req.df$roc_auc, unlist(get.aucs(test.ob[[i]])))
      req.df$accuracy <- c(req.df$accuracy, unlist(get.accuracy(test.ob[[i]])))
      req.df$sens <- c(req.df$sens, unlist(get.sens(test.ob[[i]])))
      req.df$spec <- c(req.df$spec, unlist(get.spec(test.ob[[i]])))
      req.df$f_val <- c(req.df$f_val, unlist(get.f(test.ob[[i]])))
      req.df$mcc <- c(req.df$mcc, unlist(get.mcc(test.ob[[i]])))
      req.df$pr_auc <- c(req.df$pr_auc, unlist(get.pr_aucs(test.ob[[i]])))
      req.df$classifier <- c(req.df$classifier,rep(get.class.name(names(test.ob)[i]), 
                                                   length(unlist(get.f(test.ob[[i]])))))
      req.df$feature_sel <- c(req.df$feature_sel,rep(names(test.ac)[j],
                                                     length(unlist(get.f(test.ob[[i]])))))
    }
  }
  df <- data.frame(ROC_AUC <- req.df$roc_auc,  Accuracy <- req.df$accuracy, 
                   Sensitivity <- req.df$sens, Specificity <- req.df$spec,
                   F_Value <- req.df$f_val, MCC <- req.df$mcc, PR_AUC <- req.df$pr_auc,
                   Classifier <- req.df$classifier,
                   Feature_Selection <- req.df$feature_sel)
  colnames(df) <- c('ROC_AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'F1 Score', 'MCC', 'PR-AUC',
                    'Classifier', 'Feature_Selection')
  return(df)
}

create.net.df.one <- function(res.ob, ind = 1, df.name)
{
  req.df <- list()
  req.df$auc <- c()
  req.df$accuracy <- c()
  req.df$sens <- c()
  req.df$spec <- c()
  req.df$f_val <- c()
  req.df$mcc <- c()
  req.df$classifier <- c()
  req.df$feature_sel <- c()
  req.df$df_name <- c()
  for(j in seq_along(res.ob))
  {
      class.name <- sapply(names(get.aucs(res.ob[[j]])), get.class.name)
      fea.name <- names(res.ob)[j]
      temp.ind = ind
      if(fea.name != 'varSelRF')
        temp.ind = 1
      req.df$auc <- c(req.df$auc, get.aucs(res.ob[[j]], type = 2, ind = temp.ind))
      req.df$accuracy <- c(req.df$accuracy, unlist(get.accuracy(res.ob[[j]], 
                                                                type = 2, ind = temp.ind)))
      req.df$sens <- c(req.df$sens, unlist(get.sens(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$spec <- c(req.df$spec, unlist(get.spec(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$f_val <- c(req.df$f_val, unlist(get.f(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$mcc <- c(req.df$mcc, unlist(get.mcc(res.ob[[j]], type = 2, ind = temp.ind)))
      req.df$classifier <- c(req.df$classifier, class.name)
      req.df$feature_sel <- c(req.df$feature_sel, rep(fea.name, length(class.name)))
      req.df$df_name <- c(req.df$df_name, rep(df.name, length(class.name)))
  }
  df <- data.frame(AUC <- req.df$auc,  Accuracy <- req.df$accuracy, 
                   Sensitivity <- req.df$sens, Specificity <- req.df$spec,
                   F_Value <- req.df$f_val, MCC <- req.df$mcc, Classifier <- req.df$classifier,
                   Feature_Selection <- req.df$feature_sel, Representation <- req.df$df_name)
  colnames(df) <- c('AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'F_val', 'MCC',
                   'Classifier', 'Feature_Selection', 'Representation')
  return(df)
}