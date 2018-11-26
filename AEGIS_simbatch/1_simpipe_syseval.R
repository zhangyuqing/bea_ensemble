rm(list=ls())
#setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/AEGIS_simbatch/")
setwd("~/Dropbox/Work/MultiStudy_BatchEffect/AEGIS_simbatch/")
sapply(c("sva", "MCMCpack", "BatchQC", "ROCR", "ggplot2", "limma", "nnls", 
         "glmnet", "rpart", "genefilter", "nnet", "e1071", "RcppArmadillo", "foreach", 
         "parallel", "doParallel", "MLmetrics"), require, character.only=TRUE)
#library("ranger", lib.loc="/restricted/projectnb/combat/apps/")
source("../TB/helper_new.R")
#load("AEGIS_data.RData")
load("AEGIS_data_new.RData")
#set.seed(123)


####  Parameters 
#command_args <- commandArgs(trailingOnly=TRUE)  
command_args <- c("sub", "lasso", "3", "2", "2", "TRUE", "50", "demo")
if(length(command_args)!=8){stop("Not enough input parameters!")}

##  Use data with selected genes or not
#_sel represents data with 100 top DE genes, _whole represents data without genes filtering
if(command_args[1]=="whole"){
  train_expr <- train_expr_whole
  test_expr <- test_expr_whole
  # train_expr <- test_expr_whole
  # test_expr <- train_expr_whole
  # tmp <- y_test; y_test <- y_train; y_train=tmp; rm(tmp)
}else if(command_args[1]=="sub"){
  train_expr <- train_expr_sel
  test_expr <- test_expr_sel
  # train_expr <- test_expr_sel #train_expr_sel
  # test_expr <- train_expr_sel #test_expr_sel
  # tmp <- y_test; y_test <- y_train; y_train=tmp; rm(tmp)
}else{stop("Wrong command argument 1.")}
print(dim(train_expr))

# select 1000 high variance features
# var_trn <- rowVars(train_expr)
# #var_tst <- rowVars(test_expr)
# genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:1000]]
# train_expr <- train_expr[genes_sel, ]; test_expr <- test_expr[genes_sel, ]


## Prediction model
learner_type <- command_args[2] # "lasso" 
learner_fit <- getPredFunctions(learner_type)
perf_measure_name <- "mxe" 

## Degree of batch effect (strength of signal)
N_batch <- as.numeric(command_args[3])  # 5
max_batch_mean <- as.numeric(command_args[4]) 
# median 8, range(4, 12), recommended values 0-3 (additive)
max_batch_var <- as.numeric(command_args[5]) 
# median 0.1, range(0.01, 1), recommended values 1-10 (multiplicative)
hyper_pars <- list(hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  #rep(0,5),
                   hyper_sd=sqrt(rep(0.01, N_batch)),
                   hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                     v=rep(0.01, N_batch))$alpha,  #c(3.28, 2.02, 2.845, 2.32, 2.125),
                   hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                    v=rep(0.01, N_batch))$beta)  #c(0.1824, 0.0102, 0.12, 0.0528, 0.028))
# sanity checks
if(!identical(ab2mv(a=hyper_pars$hyper_alpha, b=hyper_pars$hyper_beta)$var, rep(0.01, N_batch))){
  stop("Error in generating hyper pars for invgamma!!")
}
cat("\nBatch changes\n");
print(hyper_pars$hyper_mu);
print(hyper_pars$hyper_beta/(hyper_pars$hyper_alpha-1))

## Sample size
reduce_size <- as.logical(command_args[6])  # whether to take random subset of the original training set 
N_reduce_size <- as.numeric(command_args[7])   # max 74 
if(reduce_size){
  trainRed <- reduceSize(dat=train_expr, y=y_train, N=N_reduce_size)
  train_expr <- trainRed$dat; y_train <- trainRed$y
  
  testRed <- reduceSize(dat=test_expr, y=y_test, N=N_reduce_size)
  test_expr <- testRed$dat; y_test <- testRed$y
}
print(table(y_train))
print(table(y_test))

## Pipeline
iterations <- 3
plot_sim_ind <- NULL 
norm_data <- FALSE #as.logical(command_args[7])  # whether to normalize datasets by features
use_ref_combat <- FALSE #as.logical(command_args[8])  # whether to use ref combat to adjust test set against training set
test_item <- command_args[8]  # factor being tested, c("sig", "size", "nbatch")


####  Run pipeline
ID <- 0
pred_mat_lst <- hyper_pars_lst <- perfstats_lst <- list()
exp_name <- sprintf('%s_%s_%s_batchN%s_m%s_v%s_size%s', #_useref%s', 
                    test_item, command_args[1], learner_type, N_batch, 
                    gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T),
                    ifelse(reduce_size, N_reduce_size, 'Org'))  #, ifelse(use_ref_combat, 'T', 'F'))

while(ID < iterations){
  ID <- ID + 1
  print(paste("Simulation:", ID))
  
  ####  Spike in batch effect
  ## Split training set (AEGIS-2) in batches
  batches_ind <- splitBatch(condition=y_train, N_batch=N_batch)
  #tmp=do.call(rbind,lapply(1:N_batch,function(i){table(y_train[batches_ind[[i]]])}));print(tmp);print(apply(tmp,2,sum))
  y_sgbatch_train <- lapply(1:N_batch, function(k){y_train[batches_ind[[k]]]})
  if(any(sapply(y_sgbatch_train,table)<9)){ID <- ID - 1; next}
  batch <- rep(0, ncol(train_expr))
  for(i in 1:N_batch){batch[batches_ind[[i]]] <- i}
  
  ## Randomly shuffle batch parameters to different combinations of mean-var batch effect
  shuffled_order <- sample(1:N_batch, N_batch, replace=F)
  hyper_pars_shuffled <- c(hyper_pars[1:2], lapply(hyper_pars[3:4], function(x){x[shuffled_order]}))
  if(!identical(ab2mv(a=hyper_pars_shuffled$hyper_alpha, b=hyper_pars_shuffled$hyper_beta)$var, 
                rep(0.01, N_batch))){stop("Error in shuffling hyper pars!!")}
  
  ## Simulate batch effect 
  sim_batch_res <- simBatch(dat=train_expr, condition=y_train, batches_ind=batches_ind, 
                            batch=batch, hyper_pars=hyper_pars_shuffled)
  if(ID==1){print(round(do.call(cbind, sim_batch_res$batch_par), 3))}
  train_expr_batch <- sim_batch_res$new_dat
  #sapply(1:N_batch, function(i){mean(apply(train_expr_batch[, batches_ind[[i]]], 1, mean))})
  
  
  ####  Normalize datasets before training
  if(norm_data){
    if(ID==1){print("Normalizing data.")}
    train_expr_norm <- normalizeData(train_expr)
    test_expr_norm <- normalizeData(test_expr)
    # for train set with batch effect: normalize as a whole
    train_expr_batch_whole_norm <- normalizeData(train_expr_batch)
    # normalize within each batch
    train_expr_batch_norm <- matrix(NA, nrow=nrow(train_expr_batch), ncol=ncol(train_expr_batch), 
                                    dimnames=dimnames(train_expr_batch))
    for(k in 1:N_batch){ 
      train_expr_batch_norm[, batches_ind[[k]]] <- normalizeData(train_expr_batch[, batches_ind[[k]]])
    }
  }else{
    if(ID==1){print("Datasets are NOT normalized.")}
    train_expr_norm <- train_expr
    test_expr_norm <- test_expr
    train_expr_batch_whole_norm <- train_expr_batch_norm <- train_expr_batch
  }
  
  
  ####  Training
  ## Prediction from original training to test, without batch effect
  pred_base_res <- try(trainPipe(train_set=train_expr_norm, train_label=y_train, 
                                 test_set=test_expr_norm, lfit=learner_fit, 
                                 use_ref_combat=use_ref_combat))
  if(class(pred_base_res)=="try-error"){ID <- ID - 1; next}
  
  # ## Prediction from training WITH batch effect to test
  # pred_batch_res <- try(trainPipe(train_set=train_expr_batch_whole_norm, train_label=y_train,
  #                                 test_set=test_expr_norm, lfit=learner_fit, 
  #                                 use_ref_combat=use_ref_combat))
  # if(class(pred_batch_res)=="try-error"){ID <- ID - 1; next}
  # 
  # ##  Prediction from training after batch adjustment (Merged)  
  # train_expr_combat <- try(ComBat(train_expr_batch, batch=batch, mod=model.matrix(~y_train)))
  # if(class(train_expr_combat)=="try-error"){ID <- ID - 1; next}
  # train_expr_combat_norm <- normalizeData(train_expr_combat)
  # pred_combat_res <- try(trainPipe(train_set=train_expr_combat_norm, train_label=y_train, 
  #                                  test_set=test_expr_norm, lfit=learner_fit, 
  #                                  use_ref_combat=use_ref_combat))
  # if(class(pred_combat_res)=="try-error"){ID <- ID - 1; next}
  # 
  # ## Obtain predictions from learner trained within each batch 
  # pred_sgbatch_res <- try(lapply(1:N_batch, function(batch_id){
  #   trainPipe(train_set=train_expr_batch_norm[, batches_ind[[batch_id]]], train_label=y_sgbatch_train[[batch_id]], 
  #             test_set=test_expr_norm, lfit=learner_fit, use_ref_combat=use_ref_combat)
  # }))
  # if(class(pred_sgbatch_res)=="try-error"){ID <- ID - 1; next}
  # names(pred_sgbatch_res) <- paste0("Batch", 1:N_batch)
  # 
  # ##  Aggregate with different weights
  # pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp$pred_tst_prob)})
  # pred_mat <- do.call(cbind, pred_test_lst)
  # 
  # # Avg: simple average
  # pred_avg <- rowMeans(pred_mat)
  # 
  # # n-Avg: sample-size-weighted average 
  # pred_N_avg <- pred_mat %*% (as.matrix(sapply(batches_ind, length)) / length(y_train))
  # 
  # # CS-Avg: replicability weights
  # train_lst <- lapply(batches_ind, function(ind){train_expr_batch_norm[, ind]})
  # cs_zmat <- try(CS_zmatrix(study_lst=train_lst, label_lst=y_sgbatch_train, 
  #                           lfit=learner_fit, perf_name=perf_measure_name, 
  #                           use_ref_combat=use_ref_combat))
  # if(class(cs_zmat)=="try-error"){ID <- ID - 1; next}
  # cs_weights_seq <- CS_weight(cs_zmat)
  # pred_cs_avg <- pred_mat %*% cs_weights_seq
  # 
  # # Reg-a: use each function to predict on one study, bind predictions and do regression
  # reg_ssl_res <- try(Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train, 
  #                                 lfit=learner_fit, use_ref_combat=use_ref_combat))
  # if(class(reg_ssl_res)=="try-error"){ID <- ID - 1; next}
  # reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), n_seq=sapply(batches_ind, length))
  # pred_reg_a <- pred_mat %*% reg_a_beta
  # 
  # # Reg-s:
  # stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  # stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
  # reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x  
  # reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  # pred_reg_s <- pred_mat %*% reg_s_beta
  # 
  # 
  # ####  Evaluate performance 
  # tst_scores <- c(list(NoBatch=pred_base_res$pred_tst_prob, Batch=pred_batch_res$pred_tst_prob),
  #                 pred_test_lst,   
  #                 list(ComBat=pred_combat_res$pred_tst_prob, 
  #                      Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, 
  #                      Reg_a=pred_reg_a, Reg_s=pred_reg_s))
  # 
  # # calculate performance
  # perf_df <- as.data.frame(t(sapply(tst_scores, function(preds){
  #   if(perf_measure_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)} 
  #   # avoid Inf in computing cross-entropy loss
  #   rocr_pred <- prediction(preds, as.numeric(as.character(y_test)))
  #   perf_measure_tmp <- performance(rocr_pred, perf_measure_name)  # mean cross entropy
  #   return(as.numeric(perf_measure_tmp@y.values))
  # })))
  # 
  # 
  # ####  Output results
  # first_file <- !file.exists(sprintf('SysEval_new/%s.csv', exp_name))
  # # write.table(perf_df, sprintf('SysEval_new/%s.csv', exp_name),
  # #             append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
  # pred_mat_lst[[ID]] <- pred_mat
  # hyper_pars_lst[[ID]] <- hyper_pars_shuffled
  perfstats_lst[[ID]] <- c(logloss_trn=LogLossBinary(y_train, pred_base_res$pred_trn_prob),
                           logloss_tst=LogLossBinary(y_test, pred_base_res$pred_tst_prob),
                           accuracy_trn=AccuracyBinary(y_train, pred_base_res$pred_trn_prob),
                           accuracy_tst=AccuracyBinary(y_test, pred_base_res$pred_tst_prob),
                           #auc_trn=MLmetrics::AUC(pred_base_res$pred_trn_prob, y_train),
                           #auc_tst=MLmetrics::AUC(pred_base_res$pred_tst_prob, y_test),
                           f1_trn=MLmetrics::F1_Score(y_train, pred_base_res$pred_trn_class),
                           f1_tst=MLmetrics::F1_Score(y_test, pred_base_res$pred_tst_class))
}

perfstats_lst <- do.call(rbind, perfstats_lst)
print(perfstats_lst)
print(colMeans(perfstats_lst, na.rm=T))

#save(pred_mat_lst, hyper_pars_lst, perfstats_lst, file=sprintf("SysEval_new/predScores_%s.RData", exp_name))
