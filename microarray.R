source('helper_func.R')
source('build_models.R')
source('CV/cv_models.R')
source('get_results.R')
source('plots.R')
load('environment/prcc/micr_df.RData')
load('environment/prcc/colData_micr.RData')
load('environment/prcc/fea.RData')

micr_genes <- colnames(micr_df)
fea.micr <- list()
for(i in seq_along(fea.trial.list))
  fea.micr[[names(fea.trial.list)[i]]] <- lapply(fea.trial.list[[names(fea.trial.list)[i]]], function(g) intersect(g, micr_genes))

train.micr <- get.train.model(tr.data = micr_df, fea.list = fea.micr, stages.train = as.factor(sample_micro_info$stage), cores = 1)
cv.micr <- get.cv.model(tr.data = micr_df, fea.list = fea.micr, folds = 10,
                        stages.train = as.factor(sample_micro_info$stage), 
                        tr.model =  train.micr, cores = 1)
cv.micr.res <- get.cv.res(stages.train = as.factor(sample_micro_info$stage), 
                          cv.model = cv.micr, cores = 1)
p <- create.gridplot(cv.micr.res, size = 12)
p.micr <- ggarrange(p[[7]], p[[2]], p[[3]], p[[4]], p[[5]],p[[6]], ncol = 2, nrow = 3, 
          common.legend = T, legend = 'bottom')
ggsave('cv_micr.tiff', plot = p.micr, units = c("in"), dpi = 300, compression='lzw')
