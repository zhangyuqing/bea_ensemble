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


####  Gene-wise normalize datasets
normalizeData <- function(dat){
  dat_norm <- t(apply(dat, 1, scale, center=TRUE, scale=TRUE))
  dimnames(dat_norm) <- dimnames(dat)
  return(dat_norm)
}


#### Train pipeline: ref combat test with train, and fit learner
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
  if(learner_type=="lasso"){return(predLasso)
  }else if(learner_type=="elnet"){return(predElnet)
  }else if(learner_type=="naivebayes"){return(predNB)
  }else if(learner_type=="svm"){return(predSVM)
  }else if(learner_type=="knn"){return(predKNN) 
  }else if(learner_type=="rf"){return(predRF)
  }else if(learner_type=="nnet"){return(predNnet)
  }else if(learner_type=="plusminus"){return(predMas)
  #}else if(learner_type=="cart"){return(predCART)
  }else{stop("Method not supported!")}
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
                   lambda=seq(0,1,0.001),#10^(seq(from=-3, to=3, by=1)), 
                   type.measure="class", nfolds=10)#, intercept=FALSE)
  best_lambda <- obj$lambda.min
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

# naive bayes
# predNB <- function(
#   trn_set,
#   # gene-by-sample expression matrix for training
#   tst_set, 
#   # gene-by-sample expression matrix for test
#   y_trn
#   # response of training set, binary & numeric
# ){
#   library(e1071)
#   mod_nb <- naiveBayes(x=t(trn_set), y=as.factor(y_trn))
#   pred_train_nb <- predict(mod_nb, t(trn_set), type="raw")[,"1"]
#   pred_test_nb <- predict(mod_nb, t(tst_set), type="raw")[,"1"]
#   # If "raw", the conditional a-posterior probabilities for each class are returned
#   res <- list(pred_trn=pred_train_nb, pred_tst=pred_test_nb)
#   return(res)
# }

# SVM
predSVM <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(e1071)
  tune_ctrl <- tune.control(sampling="cross", cross=10)
  obj <- tune(svm, train.x=t(trn_set), train.y=as.factor(y_trn),
              tunecontrol=tune_ctrl,
              ranges=list(type="C-classification",
                          kernel="linear",
                          cost=exp(seq(from=-10, to=10, by=1))))
  best_cost <- obj$best.parameters[,"cost"]
  mod_svm <- svm(x=t(trn_set), y=as.factor(y_trn),
                 type="C-classification", kernel="linear",
                 cost=best_cost, probability=TRUE)
  pred_train_svm <- predict(mod_svm, t(trn_set), probability=TRUE)
  pred_train_svm <- attr(pred_train_svm, "probabilities")[,"1"]
  pred_test_svm <- predict(mod_svm, t(tst_set), probability=TRUE)
  pred_test_svm <- attr(pred_test_svm, "probabilities")[,"1"]
  
  res <- list(pred_trn_prob=pred_train_svm, pred_tst_prob=pred_test_svm)
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
                         decay=10^seq(from=-2, to=-1, by=1))
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
CS_zmatrix <- function(study_lst, label_lst, lfit, perf_name){
  n_batch <- length(study_lst)
  zmat <- matrix(0, nrow=n_batch, ncol=n_batch)
  for(i in 1:n_batch){
    for(j in 1:n_batch){
      if(i!=j){
        tmp_pred <- trainPipe(train_set=study_lst[[i]], train_label=label_lst[[i]], 
                              test_set=study_lst[[j]], lfit=lfit)$pred_tst_prob
        if(perf_name=="mxe"){tmp_pred <- pmax(pmin(tmp_pred, 1 - 1e-15), 1e-15)} # avoid Inf in computing cross-entropy loss
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
  weights_seq <- weights_seq / sum(weights_seq)
  return(weights_seq)
}


## Regression weights
# for each batch, train on each batch, prediction on it
Reg_SSL_pred <- function(study_lst, label_lst, lfit){
  SSL_pred_lst <- SSL_coef_lst <- list()
  n_batch <- length(study_lst)
  for(k in 1:n_batch){
    tmp <- lapply(1:n_batch, function(batch_id){
      return(trainPipe(train_set=study_lst[[batch_id]], train_label=label_lst[[batch_id]], 
                       test_set=study_lst[[k]], lfit=lfit)$pred_tst_prob)
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

