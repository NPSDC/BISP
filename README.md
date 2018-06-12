# BISP-PRCC - Biomarker Identification and Stage Prediction for Papillary Renal Cell Carcinoma

This repository contains the code for analysing the PRCC tumor samples according to the pipeline described in Figure 2 of the manuscript. After creating the folds in training data, groups are formed and on each group feature selection is performed. The features from each group are aggregated into feature sets AF1, AF2, AF3 and AF4 based on the number of groups they appear in. Different classifier models are built on the entire training data for all the AFs w.r.t each feature selection algorithm. These are then evaluated using 10 fold cross-validation and on an independent test dataset. A 10 fold cross validation is also done on a microarray dataset.

Following is the code organisation and its description in this repository:

**CV** - Contains scripts for the k-fold cross validation w.r.t different classifier models and a unified script performing all the above cross validations.

**feature_sel** - Contains scripts for applying the different feature selection algorithms and a unified script that applies all the feature selection algorithms.

**build_models.R** - Script for building different training models for all AFs w.r.t each feature selection algorithm.

**get_eval_mat.R** - Script for extracting the various evaluation metrics computed for all AFs w.r.t each feature selection algorithm.

**get_results.R** - Script for computing the evaluation metrics for both the predicted test samples and cross validated models.

**helper_func.R** - Script providing helper functions used in running of the other scripts.

**papillary.R** - Script analysing PRCC dataset using the scripts in the repository,

**plots.R** - Script for visualsing the results w.r.t various metrics for the different models built as a part of the pipeline.

**predict_models.R** - Script for predicitng the test  outputs from the trained models.

**microarray.R** - Script for performing 10-fold CV on microarray dataset.

