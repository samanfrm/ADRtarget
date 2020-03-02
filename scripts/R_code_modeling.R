#################################################
### Seda Arat ###################################
### Manuscript/Modeling Code ####################
### Feb 29, 2020 ################################
#################################################

### Citation:
#Robert Ietswaart<sup>\*,#</sup>, Seda Arat<sup>\*,#</sup>, Amanda X. Chen<sup>\*</sup>, 
#Saman Farahmand<sup>\*</sup>, Bumjun Kim, William DuMouchel, 
#Duncan Armstrong, Alexander Fekete, Jeffrey J. Sutherland<sup>#</sup>, Laszlo Urban<sup>#</sup>  
#*Machine learning guided association of adverse drug reactions with in vitro target-based 
#pharmacology* (2019), [BioRxiv; 750950](https://www.biorxiv.org/content/10.1101/750950v2).


# Please make sure that SOC/ and HLGT/ folders are in the same folder with this R code (R_code_modeling.R)

rm (list = ls ())

# Please set your own working directory (where this R code, SOC/ and HLGT/folders are) here
setwd ("/Users/arats/Documents/PCS_sa/")


##### Generating base matrix data and adr data for both training and test cases #####

get.bm.and.adr = function (path4bm, path4adr, path4output) {
  
  ### ADR
  ADR = read.csv (path4adr, header = T, row.names = 1, 
                  stringsAsFactors = F, na.strings = "")
  ADR[is.na (ADR)] = 0
  
  # use N_results for determining training set and testing set
  # if N_results = -1, then testing set, otherwise training set
  indeces4testing = which (ADR$N_reports < 1)
  compounds4testing = rownames (ADR) [indeces4testing]
  ADR.testing = ADR[indeces4testing, ]
  
  ADR.training = ADR[-indeces4testing, ]
  compounds4training = rownames (ADR.training)
  
  # baseMatrix
  baseMatrix = readRDS (path4bm)
  bm.training = baseMatrix[rownames (baseMatrix) %in% compounds4training, ]
  bm.testing = baseMatrix[rownames (baseMatrix) %in% compounds4testing, ]
  
  # saving files
  saveRDS (bm.training, paste0 (path4output, "bm_tr.Rdata"))
  #saveRDS (bm.testing, paste0 (path4output, "bm_te.Rdata"))
  saveRDS (ADR.training[, -1], paste0 (path4output, "adr_tr.Rdata"))
  #saveRDS (ADR.testing[, -1], paste0 (path4output, "adr_te.Rdata"))
}

get.bm.and.adr (path4bm = "./SOC/baseMatrix.Rdata", 
                path4adr = "./SOC/v1_compounds_FDA_model_format_SOC_ocr_bool.csv", 
                path4output = "./SOC/")

get.bm.and.adr (path4bm = "./HLGT/baseMatrix.Rdata", 
                path4adr = "./HLGT/v1_compounds_FDA_model_format_HLGT_ocr_bool.csv", 
                path4output = "./HLGT/")

get.bm.and.adr (path4bm = "./HLGT/BSEP_limits_30-300/baseMatrix.Rdata", 
                path4adr = "./HLGT/v1_compounds_FDA_model_format_HLGT_ocr_bool.csv", 
                path4output = "./HLGT/BSEP_limits_30-300/")

get.bm.and.adr (path4bm = "./HLGT/DiffSeed_49/baseMatrix.Rdata", 
                path4adr = "./HLGT/v1_compounds_FDA_model_format_HLGT_ocr_bool.csv", 
                path4output = "./HLGT/DiffSeed_49/")

#################################################

##### Generating random forest models #####

library (randomForest)
library (mldr)
library (utiml)
options (utiml.empty.prediction = T)

library (mltools)
#library (pROC)

#getwd ()

rf.model = function (mldr.obj, filename4model) {
  model.obj = br (mldr.obj, "RF", importance = T, seed = 7332) # please change the seed to 49 for HLGT/DiffSeed_49 model
  saveRDS (model.obj, filename4model)
  return (model.obj)
}

rf.model.prediction = function (mldr.obj, model.obj, filename4pred) {
  pred = predict (model.obj, mldr.obj, probability = T)
  saveRDS (pred, paste0 (filename4pred, ".Rdata"))
  write.csv (pred, paste0 (filename4pred, ".csv"), quote = F)
  return (pred)
}

model.prediction.metrics = function (mldr.obj, pred.obj, filename4confMtrx, ADRlabelIndeces) {
  predicted = as.bipartition (pred.obj)
  conf.matrix = multilabel_confusion_matrix (mldr.obj, predicted)
  saveRDS (conf.matrix, filename4confMtrx)
  
  actual = mldr.obj$dataset[, ADRlabelIndeces]
  metrics = t (data.frame (accuracy = round (accuracy (actual, predicted), 3), 
                           macro_precision = round (macro_precision (actual, predicted), 3), 
                           macro_recall = round (macro_recall (actual, predicted), 3), 
                           mcc = round (mcc (TP = length (which (conf.matrix$TP == TRUE)), 
                                             FP = length (which (conf.matrix$FP == TRUE)), 
                                             TN = length (which (conf.matrix$TN == TRUE)), 
                                             FN = length (which (conf.matrix$FN == TRUE))), 3), 
                           hamming_loss = round (hamming_loss (actual, predicted), 3), 
                           f1 = round (fmeasure (actual, predicted), 3)
                           ))
  
  return (metrics)
}

sequencial_split = function (mdata, r) {
  S = list()
  
  amount = trunc (r * mdata$measures$num.instances)
  indexes = c (0, cumsum (amount))
  indexes[length (r)+1] = mdata$measures$num.instances
  
  S = lapply (seq (length (r)), function (i) {
    seq (indexes[i]+1, indexes[i+1])
  })
  
  return (S)
}

crossVal = function (k.fold, mldr.obj, ADRlabelIndeces) {
  dir.create ("CrossVal_models")
  
  k5 = create_kfold_partition (mldr.obj, k = k.fold, "sequencial_split")
  saveRDS (k5, "./CrossVal_models/the_k5.Rdata")
  
  metrics = c ()
  for (k in 1:k.fold) {
    
    print (c (k, date ()))
    
    toy = partition_fold (k5, k) 
    saveRDS (toy$train, paste0 ("./CrossVal_models/",  k, "-train.Rdata"))
    saveRDS (toy$test, paste0 ("./CrossVal_models/",  k, "-test.Rdata"))
    
    mdl = rf.model (toy$train, filename4model = paste0 ("./CrossVal_models/",  k, "-the_model.Rdata"))
    
    pred.train.prob = rf.model.prediction (toy$train, mdl, filename4pred = paste0 ("./CrossVal_models/",  k, "-train_pred_prob"))
    tr.metric = model.prediction.metrics (toy$train, pred.train.prob, paste0 ("./CrossVal_models/",  k, "-train_confMatrix.Rdata"), ADRlabelIndeces)
    
    pred.test.prob = rf.model.prediction (toy$test, mdl, filename4pred = paste0 ("./CrossVal_models/",  k, "-test_pred_prob"))
    te.metric = model.prediction.metrics (toy$test, pred.test.prob, paste0 ("./CrossVal_models/",  k, "-test_confMatrix.Rdata"), ADRlabelIndeces)
    
    metrics = cbind (metrics, cbind (tr.metric, te.metric))
  }
  
  saveRDS (metrics, "./CrossVal_models/metrics.Rdata")
  return (metrics)
}

features.ADRs.predictions = function (bm.training, ADR.tr.disc, model.obj) {
  features = colnames (bm.training)
  n.features = ncol (bm.training)
  
  all.features.matrix = diag (nrow = n.features, ncol = n.features)
  rownames (all.features.matrix) = features
  colnames (all.features.matrix) = features
  
  adr = matrix (0, nrow = n.features, ncol = ncol (ADR.tr.disc))
  
  input = mldr_from_dataframe (data.frame (cbind (all.features.matrix, adr)), 
                               labelIndices = (nrow (all.features.matrix)+1):(nrow (all.features.matrix)+ncol (adr)))
  
  prediction = predict (model.obj, input)
  
  write.csv (prediction, "./Features_ADRs_predictions.csv", quote = F)
}

get.var.importance = function (model.obj, ADRlabelIndeces) {
  dir.create ("Importance_Gini")
  significant.imp.vars = c ()
  n.labels = length (ADRlabelIndeces)
  
  for (i in 1:n.labels) {
    imp = model.obj$models[[i]]$importance
    if (is.null (imp)) { # some models are empty for some reason
      print (c (i, "empty importance"))
      next
    }
    write.csv (imp, paste0 ("./Importance_Gini/imp_for_ADR-", i, ".csv"), quote = F)
    
    imp.gini = imp[, 4]
    signf.indeces = which (imp.gini > quantile (imp.gini[imp.gini > 0], probs = 0.95))
    write.csv (imp[signf.indeces, ], 
               paste0 ("./Importance_Gini/imp_for_ADR-", i, "_signf.csv"), quote = F)
    
    significant.imp.vars = c (significant.imp.vars, rownames (imp)[signf.indeces])
  }
}

threshold = function (x) {
  if (x < 0.5) {
    return (0)
  }
  else {
    return (1)
  }
}

get.summ.metrics.ADR = function () {
  actual.all = readRDS ("./adr_tr.Rdata")
  n = nrow (actual.all)
  n.models = ncol (actual.all)
  
  predicted.all = read.csv ("./the_model_pred_prob.csv", 
                            stringsAsFactors = F, row.names = 1)
  predicted.all.discr = apply (predicted.all, c (1, 2), threshold)
  
  res.table = c ()
  
  for (i in 1:n.models) {
    actual = actual.all[, i]
    predicted = predicted.all.discr[, i]
    
    tn = length (which (actual == 0 & predicted == 0))
    fp = length (which (actual == 0 & predicted == 1))
    fn = length (which (actual == 1 & predicted == 0))
    tp = length (which (actual == 1 & predicted == 1))
    metrics = data.frame (ADR = i, 
                          TN = tn, 
                          FP = fp, 
                          FN = fn, 
                          TP = tp, 
                          accuracy = round ((tp + tn) / n, 3), 
                          mcc = round (mcc (TN = tn, FN = fn, FP = fp, TP = tp), 3), 
                          precision = round (tp / (tp + fp), 3), 
                          recall = round (tp / (tp + fn), 3),
                          hamming_loss = round (hamming_loss (actual, predicted), 3), 
                          f1 = round (2 * tp / (2 * tp + fp + fn), 3)
    )
    
    res.table = rbind (res.table, metrics)
  }
  
  write.csv (res.table, "./summ_metrics_ADRs.csv", row.names = F, quote = F)
  
  return (res.table)
}

get.rf.model = function (path4res) {
  print (date ())
  setwd (path4res)
  
  bm.training = readRDS ("./bm_tr.Rdata")
  ADR.tr.disc = readRDS ("./adr_tr.Rdata")
  ADR.label.indeces = (ncol (bm.training)+1):(ncol (bm.training)+ncol (ADR.tr.disc))
  
  my.mldr = mldr_from_dataframe (as.data.frame (cbind (bm.training, ADR.tr.disc)), 
                                 labelIndices = ADR.label.indeces)
  saveRDS (my.mldr, "./the_mldr.Rdata")
  
  ### Cross Validation
  summ.metrics = crossVal (k.fold = 5, my.mldr, ADR.label.indeces)
  #print (c (date (), "5-fold CV done"))
  
  ### Modeling and Predictions
  my.rf.model = rf.model (my.mldr, filename4model = "./the_model.Rdata")
  my.pred = rf.model.prediction (my.mldr, my.rf.model, filename4pred = "./the_model_pred_prob")
  
  #print (c (date (), "Model & Predictions done."))
  
  summ.metrics = cbind (summ.metrics, model.prediction.metrics (my.mldr, my.pred, 
                                                                filename4confMtrx = "./the_confMatrix.Rdata",
                                                                ADR.label.indeces))
  cnames = c (paste0 (1:5, "-train"), paste0 (1:5, "-test"), "Model (tr)")
  colnames (summ.metrics) = cnames[c (1, 6, 2, 7, 3, 8, 4, 9, 5, 10, 11)]
  write.csv (summ.metrics, "./summ_metrics.csv", quote = F)
  
  
  ### Individual Features
  features.ADRs.predictions (bm.training, ADR.tr.disc, my.rf.model)
  
  #print (c (date (), "Features & ADRs predictions done."))
  
  ### Variable Importance
  # There are (# of ADRs) models, each of which consists of 500 trees (RF)
  # Each model is for each label (ADR)
  get.var.importance (my.rf.model, ADR.label.indeces)
  
  #print (c (date (), "Variable Importance done."))
  
  ### Summary metrics for each ADR model
  summ.metrics.ADR = get.summ.metrics.ADR ()
}

get.rf.model (path4res = "/Users/arats/Documents/PCS_sa/SOC/") # please set your working directory for SOC here

get.rf.model (path4res = "/Users/arats/Documents/PCS_sa/HLGT/") #  please set your working directory for HLGT here

get.rf.model (path4res = "/Users/arats/Documents/PCS_sa/HLGT/BSEP_limits_30-300/") #  please set your working directory for HLGT/BSEP_limits_30-300 here

# please change the seed to 49 in rf.model function for HLGT/DiffSeed_49 model below
get.rf.model (path4res = "/Users/arats/Documents/PCS_sa/HLGT/DiffSeed_49/") #  please set your working directory for HLGT/DiffSeed_49 here

#################################################

##### Re-calculating the summary metrics for random forest models #####
# in case there is discrepancey in results due to different versions of utiml package

get.summ.metrics = function (path4res) {
  
  setwd (path4res)
  
  output = c ()
  
  all.cv.files = list.files ("./CrossVal_models", pattern = "_confMatrix.Rdata")
  n.files = length (all.cv.files)
  
  for (i in 1:n.files) {
    cm = readRDS (paste0 ("./CrossVal_models/", all.cv.files[i]))
    tp = length (which (cm$TP == TRUE))
    fp = length (which (cm$FP == TRUE))
    tn = length (which (cm$TN == TRUE))
    fn = length (which (cm$FN == TRUE))
    
    metrics = data.frame (file = all.cv.files[i], 
                          TN = tn, 
                          FP = fp, 
                          FN = fn, 
                          TP = tp, 
                          accuracy = round ((tp + tn) / (tp + fp + tn + fn), 3), 
                          mcc = round (mcc (TN = tn, FN = fn, FP = fp, TP = tp), 3), 
                          precision = round (tp / (tp + fp), 3), 
                          recall = round (tp / (tp + fn), 3),
                          hamming_loss = round ((fp + fn) / (tp + fp + tn + fn), 3), 
                          f1 = round (2 * tp / (2 * tp + fp + fn), 3)
    )
    
    output = rbind (output, metrics)
  }
  
  cm = readRDS ("./the_confMatrix.Rdata")
  tp = length (which (cm$TP == TRUE))
  fp = length (which (cm$FP == TRUE))
  tn = length (which (cm$TN == TRUE))
  fn = length (which (cm$FN == TRUE))
  
  metrics = data.frame (file = "the_confMatrix.Rdata", 
                        TN = tn, 
                        FP = fp, 
                        FN = fn, 
                        TP = tp, 
                        accuracy = round ((tp + tn) / (tp + fp + tn + fn), 3), 
                        mcc = round (mcc (TN = tn, FN = fn, FP = fp, TP = tp), 3), 
                        precision = round (tp / (tp + fp), 3), 
                        recall = round (tp / (tp + fn), 3),
                        hamming_loss = round ((fp + fn) / (tp + fp + tn + fn), 3), 
                        f1 = round (2 * tp / (2 * tp + fp + fn), 3)
  )
  
  output = rbind (output, metrics)
  write.csv (output, "./summ_metrics_V2.csv", row.names = F, quote = F)
}

get.summ.metrics (path4res = "/Users/arats/Documents/PCS_sa/SOC/") # please set your working directory for SOC here

get.summ.metrics (path4res = "/Users/arats/Documents/PCS_sa/HLGT/") #  please set your working directory for HLGT here

get.summ.metrics (path4res = "/Users/arats/Documents/PCS_sa/HLGT/BSEP_limits_30-300/") #  please set your working directory for HLGT/BSEP_limits_30-300 here

get.summ.metrics (path4res = "/Users/arats/Documents/PCS_sa/HLGT/DiffSeed_49/") #  please set your working directory for HLGT/DiffSeed_49 here

#################################################

##### Chronological Validation #####

library (utiml)
setwd ("/Users/arats/Documents/PCS_sa/3classes/")

actual.2014 = readRDS ("./2014_v1_HLGT/adr_tr.rds")
pred.2014 = as.data.frame (as.bipartition (readRDS ("./2014_v1_HLGT/the_model_pred_prob.rds")))
actual.2019 = readRDS ("./2019_v1_HLGT/adr_tr.rds")

all (colnames (actual.2014) == colnames (pred.2014))
all (colnames (actual.2014) == colnames (actual.2019))

all (rownames (actual.2014) == rownames (pred.2014))
all (rownames (actual.2014) == rownames (actual.2019))

actual.2019 = actual.2019[rownames (actual.2014), colnames (actual.2014)]

get.indeces = function (row.index, col.index) {
  if (actual.2014[row.index, col.index] == 1 && pred.2014[row.index, col.index] == 0 && actual.2019[row.index, col.index] == 0) {
    return (c (row.index, col.index))
  }
  else {
    return (NA)
  }
}

cnt = 0
for (i in 1:1351) {
  #i = 1
  for (j in 1:321) {
    #j = 1
    indeces = get.indeces (i, j)
    if (!is.na (indeces)) {
      print (indeces)
      cnt = cnt + 1
    }
  }
}

# There are 13 FPs becoming TPs in 2019
#[1] 203  57
#[1] 203 291
#[1] 269 291
#[1] 314  57
#[1] 430 138
#[1] 808  57
#[1] 912 167
#[1] 944  73
#[1] 1140  272
#[1] 1140  291
#[1] 1153  303
#[1] 1209   57
#[1] 1280  193

drugs.2019 = rownames (actual.2019)[c (203, 269, 314, 430, 808, 912, 944, 1140, 1153, 1209, 1280)]

###

get.indeces.v2 = function (row.index, col.index) {
  if (actual.2014[row.index, col.index] == 1 && actual.2019[row.index, col.index] == 0) {
    return (c (row.index, col.index))
  }
  else {
    return (NA)
  }
}

count = 0
for (i in 1:1351) {
  #i = 1
  for (j in 1:321) {
    #j = 1
    indeces = get.indeces.v2 (i, j)
    if (!is.na (indeces)) {
      print (indeces)
      count = count + 1
    }
  }
}


actual.2018 = readRDS ("./v1_HLGT/adr_tr.Rdata")
common.drugs = intersect (rownames (actual.2014), rownames (actual.2018))
actual.2014 = actual.2014[common.drugs, ]
pred.2014 = pred.2014[common.drugs, ]
actual.2018 = actual.2018[common.drugs, colnames (actual.2014)]

all (colnames (actual.2014) == colnames (pred.2014))
all (colnames (actual.2014) == colnames (actual.2018))

all (rownames (actual.2014) == rownames (pred.2014))
all (rownames (actual.2014) == rownames (actual.2018))

get.indeces.2018 = function (row.index, col.index) {
  if (actual.2014[row.index, col.index] == 0 & pred.2014[row.index, col.index] == 1 & actual.2018[row.index, col.index] == 1) {
    return (c (row.index, col.index))
  }
  else {
    return (NA)
  }
}

for (i in 1:1319) {
  #i = 1
  for (j in 1:321) {
    #j = 1
    indeces = get.indeces.2018 (i, j)
    if (!is.na (indeces)) {
      print (indeces)
    }
  }
}

# There are 14 FPs becoming TPs in 2018
#[1] 195  57
#[1] 195 291
#[1] 235 261
#[1] 258 291
#[1] 303  57
#[1] 418 138
#[1] 640 291
#[1] 797  57
#[1] 1115  272
#[1] 1115  291
#[1] 1127  303
#[1] 1132  272
#[1] 1132  291
#[1] 1182   57

drugs.2018 = rownames (actual.2018)[c (195, 235, 258, 303, 418, 640, 797, 1115, 1127, 1132, 1182)]

intersect (drugs.2018, drugs.2019)

#################################################

##### Models with different Seeds for HLGT Q4_2018, different BSEP limits, SOC and HLGT Q4_2014 and HLGT Q2_2019  #####

main = function (path4res) {
  
  setwd (path4res)
  
  output = c ()
  
  all.cv.files = list.files ("./CrossVal_models", pattern = "_confMatrix.rds")
  n.files = length (all.cv.files)
  
  for (i in 1:n.files) {
    cm = readRDS (paste0 ("./CrossVal_models/", all.cv.files[i]))
    tp = length (which (cm$TP == TRUE))
    fp = length (which (cm$FP == TRUE))
    tn = length (which (cm$TN == TRUE))
    fn = length (which (cm$FN == TRUE))
    
    metrics = data.frame (file = all.cv.files[i], 
                          TN = tn, 
                          FP = fp, 
                          FN = fn, 
                          TP = tp, 
                          accuracy = round ((tp + tn) / (tp + fp + tn + fn), 3), 
                          mcc = round (mcc (TN = tn, FN = fn, FP = fp, TP = tp), 3), 
                          precision = round (tp / (tp + fp), 3), 
                          recall = round (tp / (tp + fn), 3),
                          hamming_loss = round ((fp + fn) / (tp + fp + tn + fn), 3), 
                          f1 = round (2 * tp / (2 * tp + fp + fn), 3), 
                          specificity = round (tn / (tn + fp), 3), 
                          sensitivity = round (tp / (tp + fn), 3)
    )
    
    output = rbind (output, metrics)
  }
  
  cm = readRDS ("./the_confMatrix.rds")
  tp = length (which (cm$TP == TRUE))
  fp = length (which (cm$FP == TRUE))
  tn = length (which (cm$TN == TRUE))
  fn = length (which (cm$FN == TRUE))
  
  metrics = data.frame (file = "the_confMatrix.rds", 
                        TN = tn, 
                        FP = fp, 
                        FN = fn, 
                        TP = tp, 
                        accuracy = round ((tp + tn) / (tp + fp + tn + fn), 3), 
                        mcc = round (mcc (TN = tn, FN = fn, FP = fp, TP = tp), 3), 
                        precision = round (tp / (tp + fp), 3), 
                        recall = round (tp / (tp + fn), 3),
                        hamming_loss = round ((fp + fn) / (tp + fp + tn + fn), 3), 
                        f1 = round (2 * tp / (2 * tp + fp + fn), 3), 
                        specificity = round (tn / (tn + fp), 3), 
                        sensitivity = round (tp / (tp + fn), 3)
  )
  
  output = rbind (output, metrics)
  write.csv (output, "./summ_metrics.csv", row.names = F, quote = F)
  
  return (output)
}

# please remember setting your working directory here, for all code below
main (path4res = "/Users/arats/Documents/PCS_sa/3classes/v1_HLGT_diffSeed_112358/") # please set your working directory here
main (path4res = "/Users/arats/Documents/PCS_sa/3classes/v1_HLGT_diffSeed_2662851/")
main (path4res = "/Users/arats/Documents/PCS_sa/3classes/v1_HLGT_diffSeed_49/")
main (path4res = "/Users/arats/Documents/PCS_sa/3classes/v1_HLGT_diffSeed_5332728/")

main (path4res = "/Users/arats/Documents/PCS_sa/3classes/3classes_v1HLGT_BSEPlimits_100-300/")
main (path4res = "/Users/arats/Documents/PCS_sa/3classes/3classes_v1HLGT_BSEPlimits_30-100/")
main (path4res = "/Users/arats/Documents/PCS_sa/3classes/3classes_v1HLGT_BSEPlimits_30-300/")

soc2014 = main (path4res = "/Users/arats/Documents/PCS_sa/3classes/2014_v1_SOC/")
hlgt2014 = main (path4res = "/Users/arats/Documents/PCS_sa/3classes/2014_v1_HLGT/")
hlgt2019 = main (path4res = "/Users/arats/Documents/PCS_sa/3classes/2019_v1_HLGT/")

#################################################
