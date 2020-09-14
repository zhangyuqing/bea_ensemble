####  Wrappers
ensemble_wrapper <- function(pred_test_lst, cs_zmat, reg_ssl_res, 
                             batches_ind, y_sgbatch){
  pred_mat <- do.call(cbind, pred_test_lst)
  
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  # n-Avg: sample-size-weighted average
  pred_N_avg <- pred_mat %*% (as.matrix(sapply(batches_ind, length)) / sum(sapply(batches_ind, length)))
  # CS-Avg: replicability weights
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                             n_seq=sapply(batches_ind, length))
  pred_reg_a <- pred_mat %*% reg_a_beta
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label <- do.call(c, lapply(y_sgbatch, as.character))
  reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
  reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  return(list(pred_avg=pred_avg, pred_N_avg=pred_N_avg, pred_cs_avg=pred_cs_avg,
              pred_reg_a=pred_reg_a, pred_reg_s=pred_reg_s))
}


ensemble_crossmod_wrapper <- function(preds_cm, cs_zmat_lst, cm_reg_ssl_res,
                                      batches_ind, y_sgbatch, learners){
  N_batch <- length(batches_ind)
  
  pred_mat_crossmod <- do.call(cbind, preds_cm)
  pred_mat_crossmod_reg <- lapply(1:N_batch, function(i){
    do.call(cbind, lapply(preds_cm, function(pred)pred[,i]))
  })
  pred_mat_crossmod_reg <- do.call(cbind, pred_mat_crossmod_reg)
  # Avg
  cm_avg <- rowMeans(pred_mat_crossmod)
  # n-Avg
  navg_weights <- rep((as.matrix(sapply(batches_ind, length)) / sum(sapply(batches_ind, length))), 
                      length(learners))
  cm_N_avg <- pred_mat_crossmod %*% (navg_weights / sum(navg_weights))
  # CS-Avg
  cm_cs_weights_seq <- CS_weight_crossmod(cs_zmat_lst)
  cm_cs_avg <- pred_mat_crossmod %*% cm_cs_weights_seq
  # Reg-a
  cm_reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, cm_reg_ssl_res$coef), 
                                n_seq=sapply(batches_ind, length))
  cm_reg_a <- pred_mat_crossmod_reg %*% cm_reg_a_beta
  # Reg-s
  cm_stacked_pred <- do.call(rbind, cm_reg_ssl_res$pred)
  cm_stacked_label <- do.call(c, lapply(y_sgbatch, as.character))
  cm_reg_s_beta <- nnls(A=cm_stacked_pred, b=as.numeric(cm_stacked_label))$x
  cm_reg_s <- pred_mat_crossmod_reg %*% (cm_reg_s_beta / sum(cm_reg_s_beta))
  
  return(list(cm_avg=cm_avg, cm_N_avg=cm_N_avg, cm_cs_avg=cm_cs_avg, 
              cm_reg_a=cm_reg_a, cm_reg_s=cm_reg_s))
}


perf_wrapper <- function(perf_names, tst_scores, ytest){
  perf_df <- lapply(perf_names, function(perf_name){
    as.data.frame(t(sapply(tst_scores, function(preds){
      if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  
      rocr_pred <- prediction(preds, as.numeric(as.character(ytest)))
      if(perf_name %in% c("acc", "f", "err")){
        curr_perf <- performance(rocr_pred, perf_name)  
        return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
      }else{
        curr_perf <- performance(rocr_pred, perf_name)
        return(as.numeric(curr_perf@y.values))
      }
    })))
  })
  names(perf_df) <- perf_names
  perf_df <- do.call(rbind, perf_df)
  return(perf_df)
}


ensemble_wrapper_realdata <- function(sgbatch_mods, l_type, tst_norm, 
                                      navg_weights, cs_weights_seq, reg_a_beta, reg_s_beta){
  
  pred_test_lst <- lapply(sgbatch_mods[[l_type]], function(curr_mod){
    predWrapper(curr_mod, tst_norm, l_type)
  })
  ##  Aggregate with different weights
  pred_mat <- do.call(cbind, pred_test_lst)
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  # n-Avg: sample-size-weighted average
  pred_N_avg <- pred_mat %*% navg_weights
  # CS-Avg: replicability weights
  pred_cs_avg <- pred_mat %*% cs_weights_seq[[l_type]]
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  pred_reg_a <- pred_mat %*% reg_a_beta[[l_type]]
  # Reg-s:
  pred_reg_s <- pred_mat %*% reg_s_beta[[l_type]]
  
  return(list(pred_avg=pred_avg, pred_N_avg=pred_N_avg, pred_cs_avg=pred_cs_avg,
              pred_reg_a=pred_reg_a, pred_reg_s=pred_reg_s, pred_test_lst=pred_test_lst))
}


ensemble_crossmod_wrapperNew <- function(preds_cm, N_batch, 
                                         cm_navg_weights, cm_cs_weights_seq, cm_reg_a_beta, cm_reg_s_beta){
  pred_mat_crossmod <- do.call(cbind, preds_cm)
  pred_mat_crossmod_reg <- lapply(1:N_batch, function(i){
    do.call(cbind, lapply(preds_cm, function(pred)pred[,i]))
  })
  pred_mat_crossmod_reg <- do.call(cbind, pred_mat_crossmod_reg)
  # Avg
  cm_avg <- rowMeans(pred_mat_crossmod)
  # n-Avg
  cm_N_avg <- pred_mat_crossmod %*% cm_navg_weights
  # CS-Avg
  cm_cs_avg <- pred_mat_crossmod %*% cm_cs_weights_seq
  # Reg-a
  cm_reg_a <- pred_mat_crossmod_reg %*% cm_reg_a_beta
  # Reg-s
  cm_reg_s <- pred_mat_crossmod_reg %*% cm_reg_s_beta
  
  return(list(cm_avg=cm_avg, cm_N_avg=cm_N_avg, cm_cs_avg=cm_cs_avg, 
              cm_reg_a=cm_reg_a, cm_reg_s=cm_reg_s))
}


####  Helpers related to inverse gamma distribution
mv2ab <- function(m, v){
  a <- 2 + m^2/v
  b <- m * (a-1)
  return(list(alpha=a, beta=b))
}

ab2mv <- function(a, b){
  m <- b / (a-1)
  v <- b^2 / ((a-1)^2*(a-2))
  return(list(mean=m, var=v))
}


####  Take subset of dataset
reduceSize <- function(dat, y, N){
  reduced_ctrl <- sample(which(y==0)) 
  reduced_case <- sample(which(y==1))
  #identical(sort(reduced_ctrl), which(y==0))
  
  reduced_ctrl <- reduced_ctrl[1:N] 
  reduced_case <- reduced_case[1:N]
  
  reduced_indices <- c(reduced_ctrl, reduced_case)
  dat <- dat[, reduced_indices]
  y <- y[reduced_indices]
  return(list(dat=dat, y=y))
}


####  Split a dataset into man-made batches
splitBatch <- function(condition, N_batch){
  # split samples into case / control groups
  case_ind <- which(condition==1)
  ctrl_ind <- which(condition==0)
  
  # split each condition group into N_batch batches
  batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))
  batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))
  #print((sum(sapply(batches_ind_case,length))==length(case_ind)) & (sum(sapply(batches_ind_ctrl,length))==length(ctrl_ind)))
  
  # combine case / control samples in each batch
  batches_ind <- list()
  for(i in 1:N_batch){
    batches_ind[[i]] <- sort(c(batches_ind_case[[i]], batches_ind_ctrl[[i]]))
  }
  return(batches_ind)
}


####  Take subset from dataset
subsetBatch <- function(condition, N_sample_size, N_batch){
  # split samples into case / control groups
  case_ind <- which(condition==1)
  ctrl_ind <- which(condition==0)
  
  # number of controls and cases to take
  N_ctrl <- N_sample_size / 2
  N_case <- N_sample_size / 2
  if(N_ctrl*N_batch > length(ctrl_ind) | N_case*N_batch > length(case_ind)){stop("Not enough samples to subset!")}
  
  # split each condition group into N_batch batches
  batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))
  while(any(sapply(batches_ind_case, length) < N_case)){batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))}
  batches_ind_case <- lapply(batches_ind_case, function(x){x[1:N_case]})
  
  batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))
  while(any(sapply(batches_ind_ctrl, length) < N_ctrl)){batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))}
  batches_ind_ctrl <- lapply(batches_ind_ctrl, function(x){x[1:N_ctrl]})
  
  # combine case / control samples in each batch
  batches_ind <- list()
  for(i in 1:N_batch){
    batches_ind[[i]] <- sort(c(batches_ind_case[[i]], batches_ind_ctrl[[i]]))
  }
  return(batches_ind)
}


####  Simulate batch effect based on ComBat assumption
simBatch <- function(dat, condition, batches_ind, batch, hyper_pars){
  n_batches <- sapply(batches_ind, length) # number of samples in each batch
  n_genes <- nrow(dat)
  
  ## Organize hyper batch parameters
  batch_par <- list()
  for(i in 1:length(n_batches)){
    batch_par[[i]] <- sapply(hyper_pars, function(item){item[i]}) # mean, sd of gaussian; alpha, beta for InvGamma
  }
    
  ## Simulate batch parameters from hyper-pars
  gamma <- delta2 <- list()
  for(i in 1:length(n_batches)){
    gamma[[i]] <- rnorm(n_genes, mean=batch_par[[i]]["hyper_mu"], sd=batch_par[[i]]["hyper_sd"])
    delta2[[i]] <- rinvgamma(n_genes, shape=batch_par[[i]]["hyper_alpha"], scale=batch_par[[i]]["hyper_beta"])
  }
    
  ## Simulate batch effect
  # fit linear model to data with no batch parameters, calculate residual variance
  X <- model.matrix(~Condition, data=data.frame(Condition=condition))
  beta <- solve(t(X) %*% X) %*% t(X) %*% t(dat)
  resid <- dat - t(X %*% beta)
  #range(apply(resid,1,mean)); range(apply(resid,1,var))
  
  # spike-in batch variance: multiply by condition adjusted data with delta
  resid_varbatch <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
  for(j in 1:length(n_batches)){
    curr_resid <- resid[, batches_ind[[j]]]
    spikein_var <- lapply(1:n_batches[j], function(col_ind){curr_resid[, col_ind] * sqrt(delta2[[j]])})
    resid_varbatch[, batches_ind[[j]]] <- do.call(cbind, spikein_var) 
  }
  #sapply(1:5, function(k){mean(apply(resid[, batches_ind[[k]]],1,var))})
  #sapply(1:5, function(k){mean(apply(resid_varbatch[, batches_ind[[k]]],1,var))})
  
  # construct mean batch parameter design matrix using gamma
  X_batch <- model.matrix(~-1+Batch, data=data.frame(Batch=factor(batch)))
  gamma_vec <- do.call(rbind, gamma)  #apply(gamma_vec,1,mean)
  
  # new data with added batch effect
  new_dat <- t(cbind(X, X_batch) %*% rbind(beta, gamma_vec)) + resid_varbatch
  if(!identical(rownames(new_dat), rownames(dat))){stop("BUG in simBatch function!")
  }else{colnames(new_dat) <- colnames(dat)}
  
  res <- list(new_dat=new_dat, batch_par=batch_par)
  return(res)
}


####  Gene-wise normalize datasets (z-score scaling)
normalizeData <- function(dat){
  dat_norm <- t(apply(dat, 1, scale, center=TRUE, scale=TRUE))
  dimnames(dat_norm) <- dimnames(dat)
  return(dat_norm)
}


#### Train pipeline: ref combat test with train, and fit learner
#train_set=train_expr_norm; train_label=y_train; test_set=test_expr_norm; lfit=learner_fit
trainPipe <- function(train_set, test_set, train_label, lfit=learner_fit, use_ref_combat=TRUE){
  if(use_ref_combat){  # use ref combat to adjust test set to match train set
    cmb_dat <- cbind(test_set, train_set)
    tmp_batch <- c(rep(2,ncol(test_set)), rep(1,ncol(train_set)))
    test_refadj <- ComBat(cmb_dat, batch=tmp_batch, mod=NULL, ref.batch=1)[, 1:ncol(test_set)]
    if(!identical(dim(test_refadj), dim(test_set))){stop("Error in trainPipe!")}
  }else{
    test_refadj <- test_set
  }
  pred_res <- lfit(trn_set=train_set, y_trn=train_label, tst_set=test_refadj)
  return(pred_res)
}


####  Get functions corresponding to learner type
getPredFunctions <- function(learner_type){
  if(learner_type=="lasso"){return(predLasso_pp)  #return(predLasso)
  }else if(learner_type=="elnet"){return(predElnet)
  }else if(learner_type=="naivebayes"){return(predNB)
  }else if(learner_type=="svm"){return(predSVM)
  }else if(learner_type=="knn"){return(predKNN) 
  }else if(learner_type=="rf"){return(predRF_pp)  #return(predRF)
  }else if(learner_type=="rf_fs"){return(predRF_fs_pp) 
  }else if(learner_type=="nnet"){return(predNnet_pp)  #return(predNnet)
  }else if(learner_type=="plusminus"){return(predMas)
  }else{stop("Method not supported!")}
}


####  Some performance metrics 
LogLossBinary <- function(actual, predicted, eps=1e-15) {
  if(class(actual)=="factor"){actual <- as.numeric(as.character(actual))}
  predicted = pmin(pmax(predicted, eps), 1-eps)
  - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
}

AccuracyBinary <- function(actual, predicted) {
  if(class(actual)=="factor"){actual <- as.numeric(as.character(actual))}
  predicted_bi <- as.numeric(predicted >= 0.5)
  sum(predicted_bi==actual) / length(actual)
}


####  Prediction functions
# lasso
predLasso <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set,
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(glmnet)
  obj <- cv.glmnet(x=t(trn_set), y=factor(y_trn), family="binomial", alpha=1,
                   lambda=10^(seq(from=-4, to=4, by=0.1)),
                   type.measure="mse", nfolds=10)#, intercept=FALSE)
  best_lambda <- obj$lambda.1se
  mod_logit <- glmnet(x=t(trn_set), y=factor(y_trn), family="binomial", alpha=1,
                      lambda=best_lambda)#, intercept=FALSE)

  # predictions
  pred_train_prob <- predict(mod_logit, t(trn_set), type="response")[,1]
  pred_test_prob <-  predict(mod_logit, t(tst_set), type="response")[,1]
  #type "response" gives the fitted probabilities for "binomial",
  #type "class" produces the class label corresponding to the maximum probability

  res <- list(beta=c(mod_logit$a0, mod_logit$beta[,1]),
              pred_trn_prob=pred_train_prob, pred_tst_prob=pred_test_prob)#,
              #pred_trn_class=pred_train_class, pred_tst_class=pred_test_class)
  return(res)
}

# elastic net
# predElnet <- function(
#   trn_set,
#   # gene-by-sample expression matrix for training
#   tst_set, 
#   # gene-by-sample expression matrix for test
#   y_trn 
#   # response of training set, binary & numeric
# ){
#   library(caret)
#   parGrid <- expand.grid(lambda=exp(seq(from=-10, to=10, by=1)),
#                          alpha=seq(from=0, to=1, by=0.1))
#   ctrl <- trainControl(method = "cv", number=4)
#   mod_elnet <- train(x=t(trn_set), y=as.factor(y_trn), family = "binomial",
#                      method="glmnet",
#                      trControl=ctrl,
#                      tuneGrid=parGrid)
#   pred_train_elnet <- predict(mod_elnet, t(trn_set), type="prob")[,"1"] 
#   pred_test_elnet <- predict(mod_elnet, t(tst_set), type="prob")[,"1"] 
#   # either "raw" or "prob", for the number/class predictions or class probabilities
#   res <- list(pred_trn=pred_train_elnet, pred_tst=pred_test_elnet)
#   return(res)
# }

# SVM
predSVM <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set = NULL, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(e1071)
  tune_ctrl <- tune.control(sampling="cross", cross=4)#, error.fun=LogLossBinary)
  obj <- tune(svm, train.x=t(trn_set), train.y=as.factor(y_trn),
              tunecontrol=tune_ctrl,
              ranges=list(type="C-classification",
                          kernel="linear",
                          cost=exp(seq(from=-10, to=10, by=2))))
  best_cost <- obj$best.parameters[,"cost"]
  mod_svm <- svm(x=t(trn_set), y=as.factor(y_trn),
                 type="C-classification", kernel="linear",
                 cost=best_cost, probability=TRUE)
  pred_train_svm <- predict(mod_svm, t(trn_set), probability=TRUE)
  pred_train_svm <- attr(pred_train_svm, "probabilities")[,"1"]
  if(!is.null(tst_set)){
    pred_test_svm <- predict(mod_svm, t(tst_set), probability=TRUE)
    pred_test_svm <- attr(pred_test_svm, "probabilities")[,"1"]
  }else{
    pred_test_svm <- NULL
  }
  res <- list(mod=mod_svm, pred_trn_prob=pred_train_svm, pred_tst_prob=pred_test_svm)
  return(res)
}

# random forest
predRF <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(caret)
  
  training_df <- data.frame(t(trn_set), as.factor(y_trn))
  colnames(training_df) <- c(paste("gene", 1:nrow(trn_set), sep=""), "response")
  rownames(training_df) <- 1:ncol(trn_set)
  
  test_df <- data.frame(t(tst_set))
  colnames(test_df) <- paste("gene", 1:nrow(tst_set), sep="")
  rownames(test_df) <- 1:ncol(tst_set)
  
  ctrl <- trainControl(method="cv", number=10)
  f <- as.formula(paste("response ~ ", paste(colnames(training_df)[-ncol(training_df)], collapse= "+",sep="")))
  mod_rf <- train(form=f, data=training_df, 
                  method="rf", metric="Accuracy", 
                  tuneLength=5, trControl=ctrl)
  
  pred_train_rf <- predict(mod_rf, training_df, type="prob")[,"1"]
  pred_test_rf <- predict(mod_rf, test_df, type="prob")[,"1"]
  
  res <- list(pred_trn_prob=pred_train_rf, pred_tst_prob=pred_test_rf)
  return(res)
}

# neural net 
predNnet <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(caret)
  parGrid <- expand.grid(size=seq(from=2, to=3, by=1),
                         decay=10^seq(from=-4, to=-1, by=0.5))
  ctrl <- trainControl(method = "cv", number=10)
  mod_nnet <- train(x=t(trn_set), y=as.factor(y_trn), maxit=1000, MaxNWts=50000,
                    method="nnet", trace=FALSE,
                    trControl=ctrl, softmat=TRUE,
                    tuneGrid=parGrid)
  pred_train_nnet <- predict(mod_nnet, t(trn_set), type="prob")[,"1"]
  pred_test_nnet <- predict(mod_nnet, t(tst_set), type="prob")[,"1"] 
  
  res <- list(pred_trn_prob=pred_train_nnet, pred_tst_prob=pred_test_nnet)
  return(res)
}

# mas-o-menos 
predMas <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  trn_set_norm <- t(scale(t(trn_set), center=TRUE, scale=TRUE))
  tst_set_norm <- t(scale(t(tst_set), center=TRUE, scale=TRUE))
  
  training_df_norm <- data.frame(t(trn_set_norm), as.factor(y_trn))
  colnames(training_df_norm) <- c(paste("gene", 1:nrow(trn_set_norm), sep=""), "response")
  rownames(training_df_norm) <- 1:ncol(trn_set_norm)
  
  alpha <- rep(0, nrow(trn_set_norm))
  for(j in 1:nrow(trn_set_norm)){
    f <- as.formula(paste("response ~ 0 +", paste(colnames(training_df_norm)[j], 
                                                  collapse= "+",sep="")))
    ctr <- glm.control(maxit=1000)
    mod_tmp <- glm(f, data=training_df_norm, family=binomial, control=ctr)
    alpha[j] <- coef(mod_tmp)
  }
  v <- (2*(alpha>0)-1)/sqrt(nrow(trn_set))
  
  pred_train_plusminus <- as.numeric(t(trn_set_norm) %*% as.matrix(v))
  pred_train_plusminus <- 1/(1+exp(- pred_train_plusminus))
  pred_test_plusminus <- as.numeric(t(tst_set_norm) %*% as.matrix(v))
  pred_test_plusminus <- 1/(1+exp(- pred_test_plusminus))
  
  res <- list(pred_trn_prob=pred_train_plusminus, pred_tst_prob=pred_test_plusminus)
  return(res)
}


####  Ensemble with different weighting methods
## CS-Avg
# train in each batch, predict on all other batches
CS_zmatrix <- function(study_lst, label_lst, lfit, perf_name, use_ref_combat){
  n_batch <- length(study_lst)
  zmat <- matrix(0, nrow=n_batch, ncol=n_batch)
  for(i in 1:n_batch){
    for(j in 1:n_batch){
      if(i!=j){
        tmp_pred <- trainPipe(train_set=study_lst[[i]], train_label=label_lst[[i]], 
                              test_set=study_lst[[j]], 
                              lfit=lfit, use_ref_combat=use_ref_combat)$pred_tst_prob
        if(perf_name=="mxe"){tmp_pred <- pmax(pmin(tmp_pred, 1 - 1e-15), 1e-15)}
        # avoid Inf in computing cross-entropy loss
        rocr_pred <- prediction(tmp_pred, as.numeric(as.character(label_lst[[j]])))
        perf_tmp <- performance(rocr_pred, perf_name)  # mean cross entropy
        zmat[i,j] <- as.numeric(perf_tmp@y.values)
      }
    }
  }
  return(zmat)
}

# calculate weights for CS-Avg
CS_weight <- function(cs_zmat){
  z_seq <- sqrt(rowSums(cs_zmat) / (nrow(cs_zmat)-1))
  weights_seq <- abs(z_seq - max(z_seq))
  if(sum(weights_seq)==0){weights_seq <- rep(1, length(z_seq))}
  weights_seq <- weights_seq / sum(weights_seq)
  return(weights_seq)
}

CS_weight_crossmod <- function(cs_zmat_lst){
  z_seq_lst <- lapply(cs_zmat_lst, function(cs_zmat){sqrt(rowSums(cs_zmat) / (nrow(cs_zmat)-1))})
  z_seq <- do.call(c, z_seq_lst)
  weights_seq <- abs(z_seq - max(z_seq))
  if(sum(weights_seq)==0){weights_seq <- rep(1, length(z_seq))}
  weights_seq <- weights_seq / sum(weights_seq)
  return(weights_seq)
}

## Regression weights
# for each batch, train on each batch, prediction on it
Reg_SSL_pred <- function(study_lst, label_lst, lfit, use_ref_combat){
  SSL_pred_lst <- SSL_coef_lst <- list()
  n_batch <- length(study_lst)
  for(k in 1:n_batch){
    tmp <- lapply(1:n_batch, function(batch_id){
      return(trainPipe(train_set=study_lst[[batch_id]], train_label=label_lst[[batch_id]], 
                       test_set=study_lst[[k]], 
                       lfit=lfit, use_ref_combat=use_ref_combat)$pred_tst_prob)
    })
    names(tmp) <- paste0("Batch", 1:n_batch)
    SSL_pred_lst[[k]] <- do.call(cbind, tmp)
    coef_k <- nnls(A=SSL_pred_lst[[k]], b=as.numeric(as.character(label_lst[[k]])))$x  
    SSL_coef_lst[[k]] <- coef_k
  }
  return(list(pred=SSL_pred_lst, coef=SSL_coef_lst))
}

Reg_a_weight <- function(coef_mat, n_seq){
  weights_seq <- rep(0, length(n_seq))
  for(i in 1:length(n_seq)){
    weights_seq <- weights_seq + coef_mat[i, ] * n_seq[i]
  }
  weights_seq <- weights_seq / sum(weights_seq)
  return(weights_seq)
}

#study_lst=train_lst; label_lst=y_sgbatch_train; learner_lst=learner_types
crossmod_Reg_SSL_pred <- function(study_lst, label_lst, learner_lst, use_ref_combat){
  SSL_pred_lst <- SSL_coef_lst <- list()
  n_batch <- length(study_lst)
  for(k in 1:n_batch){
    tmp <- lapply(1:n_batch, function(batch_id){
      cm_tmp <- list()
      for(lfit_name in learner_lst){
        lfit <- getPredFunctions(lfit_name)
        cm_tmp[[lfit_name]] <- trainPipe(train_set=study_lst[[batch_id]], train_label=label_lst[[batch_id]], 
                                         test_set=study_lst[[k]], lfit=lfit, use_ref_combat=use_ref_combat)$pred_tst_prob
      }
      return(do.call(cbind, cm_tmp))
    })
    names(tmp) <- paste0("Batch", 1:n_batch)
    SSL_pred_lst[[k]] <- do.call(cbind, tmp)
    coef_k <- nnls(A=SSL_pred_lst[[k]], b=as.numeric(as.character(label_lst[[k]])))$x  
    SSL_coef_lst[[k]] <- coef_k
  }
  return(list(pred=SSL_pred_lst, coef=SSL_coef_lst))
}


####  Prediction functions slightly organized from Prasad Patil's github
predLasso_pp <- function(trn_set, tst_set=NULL, y_trn, ...){
  mod <- glmnet::cv.glmnet(x=t(trn_set), y=as.numeric(as.character(y_trn)), family="binomial", ...)
  pred_trn_prob <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="response"))
  pred_trn_class <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="class"))
  if(!is.null(tst_set)){
    pred_tst_prob <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
    pred_tst_class <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="class"))
  }else{
    pred_tst_prob <- NULL; pred_tst_class <- NULL
  }
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

#trn_set=train_set; y_trn=train_label; tst_set=test_refadj
predRF_pp <- function(trn_set, tst_set=NULL, y_trn){
  data <- data.frame(y=as.factor(y_trn), t(trn_set))
  mod <- ranger::ranger(y ~ ., data = data, write.forest=TRUE)
  mod_prob <- ranger::ranger(y ~ ., data = data, write.forest=TRUE, probability=TRUE)
  pred_trn_prob <- predict(mod_prob, data = data.frame(t(trn_set)))$predictions[, "1"]
  pred_trn_class <- predict(mod, data = data.frame(t(trn_set)))$predictions
  if(!is.null(tst_set)){
    newdata <- data.frame(t(tst_set))
    pred_tst_prob <- predict(mod_prob, data = newdata)$predictions[, "1"]
    pred_tst_class <- predict(mod, data = newdata)$predictions
  }else{
    pred_tst_prob <- NULL; pred_tst_class <- NULL
  }
  return(list(mod = mod_prob, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

predRF_fs_pp <- function(trn_set, tst_set, y_trn){
  data <- data.frame(y=as.factor(y_trn), t(trn_set))
  newdata <- data.frame(t(tst_set))
  tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
  data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
  mod <- ranger(y ~ ., data = data, write.forest = TRUE, probability = T)
  predict(mod, data=newdata)$predictions[,2]
}

predNnet_pp <- function(trn_set, tst_set, y_trn){
  data <- data.frame(y=as.factor(y_trn), t(trn_set))
  newdata <- data.frame(t(tst_set))
  mod <- nnet::nnet(y ~ ., data = data, size = 10, MaxNWts = 10000, linout = F, trace = F)
  pred_trn_prob <- as.vector(predict(mod, newdata = data[,-1]))
  pred_trn_class <- as.vector(predict(mod, newdata = data[,-1], type="class"))
  pred_tst_prob <- as.vector(predict(mod, newdata = newdata))
  pred_tst_class <- as.vector(predict(mod, newdata = newdata, type="class"))
  return(list(pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}


predWrapper <- function(mod, tst_set, function_name){
  if(function_name=='lasso'){
    res <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
  }else if(function_name=='rf'){
    newdata <- data.frame(t(tst_set))
    res <- predict(mod, data = newdata)$predictions[, "1"]
  }else if(function_name=='svm'){
    res <- predict(mod, t(tst_set), probability=TRUE)
    res <- attr(res, "probabilities")[,"1"]
  }
  return(res)
}
