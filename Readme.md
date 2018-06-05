File Description
Feature_Extract.R - contains all the necessary function for extracting the raw expression data and storing them into a matrix of standard form. Also map information  patient and file ids.

extract_stage.R - Returning a data frame that combines stage and type wise information

finding_deg.R

get_features.R - Uses Feature_Extract to extract features for counts and fpqm 

gen_ext.R

onlytumors.R - creates data frames/dds objects for only tumor samples of exp_prof,exp_fpqm/log
Also separates them into only those for which we stage data.
Uses function.R and it is in this files we call functions and save images for visualisation
sample_info.R

Visualisation.R

##Environment Variables and proj variables
##Working for papillary
exp_prof - Expression matrix corresponding to raw counts
exp_fpqm - Expression matrix corresponding to fpqm-ql values
ens.ids.all - All ens ids in file matrix
ens.ids - ens ids corresponding to protein coding regions
genes.entrez - mapping between ens.ids and genes.entrez
sample.file.map - mapping between count files and patient ids
pap.comb - data frame that tells the stage a patient is in along with its type
int.type.stage.indexes - contains the indexes according pap.comb$name. It is a list which contains intersection of stage and type
df.stages - Coldata for dds_tumor, similar to pap.comb but does not contain type
dds_tumor - DESeq object for exp_prof_tumor
dds_tumor_reported - DESeq object for exp_prof_tumor_reported
