source('helper_func.R')
source('get_eval_mat.R')
library(ggpubr)

#' Runs the repeated commands across the various ggplots
#' 
#' @param p.g Plot corresponding to a particular measure
#' @param df dataframe containing the requried balues
#' @param size text.size
#' @param labs labels for annotation
#' @param y.max Vector containing the maximum height to be taken in respective graph
#' @return modified p.g with common repeated modifications
create.rep <- function(p.g, df, size, labs, y.max)
{
  scaleFUN <- function(x) sprintf("%.2f", x)
  p.g <- p.g + geom_boxplot() +
    scale_fill_hue(labels = c(bquote(paste('DESeq2 | ', log[2], 'FC | '>=' 1')),
                                                 bquote(paste('DESeq2 | ', log[2], 'FC | '>=' 1.5')), 
                                                 bquote(paste('DESeq2 | ', log[2], 'FC | '>=' 2')),
                                                 bquote(paste('SAMseq | ', log[2], 'FC | '>=' 1')), 
                                                 bquote(paste('SAMseq | ', log[2], 'FC | '>=' 1.5')),
                                                 bquote(paste('SAMseq | ', log[2], 'FC | '>=' 2')), 'Shrunken Centroids', 'varSelRF')) +
    theme(legend.text=element_text(size = 12), legend.title=element_text(size = 12, face = 'bold')) +
    theme(panel.background = element_rect(fill = 'white', color = 'black')) +
    theme(axis.text.x = element_text(size = size, color = 'black'), axis.title.x = element_text(size = size)) +
    theme(axis.text.y = element_text(size = size, color = 'black'), axis.title.y = element_text(size = size)) +
    scale_y_continuous(labels=scaleFUN) #+ scale_x_discrete(labels = c('Naive Bayes' = ))
    geom_text(aes_string(label=labs, y = y.max), position = position_dodge(0.75), size = 2.3)
  return(p.g)
}

#' Gets the plots for the various evaluation metrics
#' 
#' @param results List returned by get.test.results or get.cv.results
#' @param size Font size to be used for labeling
#' @return List of plots for the individual evaluation metric
create.gridplot <- function(results, size = 12)
{
   
  net.df <- create.net.df(results)
  colnames(net.df)[c(5,9)] = c('F1 Score', 'Feature Selection')
  # y.pos <- rep(0,160)
  # for(i in seq(1,160,4))
  #   y.pos[i:(i+3)] = max(net.df$AUC[i:(i+3)]) + 0.004
    
  net.df[,9] = factor(net.df[,9], levels=levels(net.df[,9])[c(2,1,3,5,4,6:8)])
  #net.df$labs <- rep(c('G','H', 'A', 'B', 'C', 'D', 'E', 'F'), each = 20)
  #net.df <- net.df[order(net.df[,7], net.df[,6]),]
  
  # x.start = 0.672
  # x.delta = 0.093
  # line.size = 0.1
  # line.type = 'dotted'
  
  #y.max <- sapply(seq(1:5), function(i) max(net.df[,i]) + sd(net.df[,i])/1.2)
  p1 <- ggplot(net.df, aes(x=Classifier, y=`ROC_AUC`, fill=`Feature Selection`))  
  p1 <- create.rep(p1, net.df, size, 'labs', y.max[1])
 
  p2 <- ggplot(net.df, aes(x=Classifier, y=Accuracy, fill=`Feature Selection`)) 
  p2 <- create.rep(p2, net.df, size, 'labs', y.max[2])  

  p3 <- ggplot(net.df, aes(x=Classifier, y=Sensitivity, fill=`Feature Selection`)) 
  p3 <- create.rep(p3, net.df, size, 'labs', y.max[3])  
  
  p4 <- ggplot(net.df, aes(x=Classifier, y=Specificity, fill=`Feature Selection`)) 
  p4 <- create.rep(p4, net.df, size, 'labs', y.max[4])
  
  p5 <- ggplot(net.df, aes(x=Classifier, y=`F1 Score`, fill=`Feature Selection`)) 
  p5 <- create.rep(p5, net.df, size, 'labs', y.max[5]) 
  
  p6 <- ggplot(net.df, aes(x=Classifier, y=`MCC`, fill=`Feature Selection`)) 
  p6 <- create.rep(p6, net.df, size, 'labs', y.max[6]) 
  
  p7 <- ggplot(net.df, aes(x=Classifier, y=`PR-AUC`, fill=`Feature Selection`)) 
  p7 <- create.rep(p7, net.df, size, 'labs', y.max[7]) 
  
  p <- list(p1,p2,p3,p4,p5,p6,p7)
   # for(j in seq(5))
   # {
   #   for(i in seq(8))
   #   {
   #      p <- lapply(seq(p), function(ind)
   #      {
   #        p[[ind]] <- p[[ind]] + geom_segment(aes_string(x = x.start + (j-1)*1.001 + (i-1)*0.093, 
   #                                         y = max(net.df[,ind][(i-1)*20+ (j-1)*4 + 1:4]) + 0.002,
   #                                         xend =x.start + (j-1)*1.001 + (i-1)*0.093, 
   #                                         yend = y.max[ind] - 0.005),
   #                                         linetype = line.type, size = line.size)  
   #      })
   #   }
   # }
  
  return(p)
}


# load('../papillary/environment/accuracy_feature/updated/cv_mirc.RData')
# load('../papillary/environment/accuracy_feature/updated/new_data/test_results.RData')
# load('../papillary/environment/accuracy_feature/updated/new_data/cv_results.RData')
# labels = c('DESeq2 |logFC > 1|', 'DESeq2 |logFC > 1.5|', 'DESeq2 |logFC > 2|',
#            'SAMseq |logFC > 1|', 'SAMseq |logFC > 1.5|', 'SAMseq |logFC > 2|',
#            'Shrunken Centroids', 'varSelRF')
# p <- create.gridplot(cv.micr.res,  size = 12)
# p.micr <- ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], ncol = 2, nrow = 3, 
#                       common.legend = T, legend = 'bottom')
# ggsave('micr_results.eps', plot = p.micr, units = c("in"), dpi = 300, compression = 'lzw')
# 
# p <- create.gridplot(test.trial.results, size = 12)
# p.test <- ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], ncol = 2, nrow = 3, 
#                      common.legend = T, legend = 'bottom')
# ggsave('test_results.eps', plot = p.test, units = c("in"), dpi = 300, compression = 'lzw')

# p <- create.gridplot(cv.trial.results, size = 12)
# p.cv <- ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], ncol = 2, nrow = 3, 
#                      common.legend = T, legend = 'bottom')
# ggsave('cv_results.eps', plot = p.cv, units = c("in"), dpi = 300)
# tiff("cv_results.tif", res=1200, compression = "lzw")
# p.im
# dev.off()
