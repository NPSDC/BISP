source('helper_func.R')
source('feature_sel/feature_sel_list.R')
source('build_models.R')
source('CV/cv_models.R')
source('predict_models.R')
source('get_results.R')
source('plots.R')

load('environment/prcc/prcc_norm.RData')
load('environment/prcc/counts.RData')
load('environment/prcc/stages_prcc.RData')


gr <- build.groups(nrow(prcc.norm.data), 5)
get.stage.distribution(gr, stages.prcc)
test.ind <- sort(unlist(gr[[5]]))
test.data <- prcc.norm.data[test.ind,]
stages.prcc.test <- stages.prcc[test.ind]

train.ind <- sort(unlist(gr[c(1:4)]))
train.data <- prcc.norm.data[train.ind,]
stages.prcc.train <- stages.prcc[train.ind]
gr.train <- build.groups(nrow(train.data), 4)

cores = detectCores() - 1

pap.fea.ob <- get.features.object(train.count.data = counts.data[train.ind,]),
                                  train.norm.data = train.data, train.stages = stages.prcc.train,
                                  train.ind.list = gr.train, cores.list = c(cores, cores, cores, 1),
                                  type = 1)
save(pap.fea.ob, file = 'environment/prcc/pap_fea_ob.RData')
pap.fea.list <- get.class.fea(pap.fea.ob)
save(pap.fea.list, file = 'environment/prcc/pap_fea_list.RData')
pap.train.model <- get.train.model(tr.data = train.data, fea.list = pap.fea.list, stages.train = stages.prcc.train, 
                                 cores = 1)
save(pap.train.model, file = 'environment/prcc/pap_train_model.RData')
##2 cores used since using more than this was not optimal with the given machine
pap.cv.model <- get.cv.model(tr.data = train.data, fea.list = pap.fea.list, folds = 10, 
                       stages.train = stages.prcc.train, tr.model = pap.train.model, cores = 1)
save(pap.cv.model, file = 'environment/prcc/pap_cv_model.RData')
pap.test.pred <- get.test.pred(tr.data = train.data, te.data = test.data, fea.list = pap.fea.list, 
                               stages.tr = stages.prcc.train, tr.model = pap.train.model, cv.model = pap.cv.model, cores = 1
                                 )
save(pap.test.pred, file = 'environment/prcc/pap_test_pred.RData')
pap.cv.res <- get.cv.res(stages.train = stages.prcc.train, cv.model = pap.cv.model, cores = 1)
save(pap.cv.res, file = 'environment/prcc/pap_cv_res.RData')
pap.test.res <- get.test.res(stages.test = stages.prcc.test, test.pred = pap.test.pred, cores = 1)
save(pap.test.res, file = 'environment/prcc/pap_test_res.RData')
