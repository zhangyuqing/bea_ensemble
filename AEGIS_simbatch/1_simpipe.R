rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/AEGIS_simbatch/")
#setwd("~/Dropbox/Work/MultiStudy_BatchEffect/AEGIS_simbatch/")
sapply(c("sva", "MCMCpack", "BatchQC", "ROCR", "ggplot2", "limma", "nnls"), require, character.only=TRUE)
source("helper.R")
load("AEGIS_data.RData")


####  Parameters 
command_args <- commandArgs(trailingOnly=TRUE)  
# command_args <- c("nnet", "TRUE", "TRUE")
# prediction
learner_type <- command_args[1] # "lasso" 
learner_fit <- getPredFunctions(learner_type)
perf_measure_name <- "mxe" #command_args[3]

# batch
N_batch <- 5 #as.numeric(command_args[1])
hyper_pars <- list(hyper_mu=c(-1, -0.5, 0, 0.5, 1),  #rep(0,5),
                   hyper_sd=sqrt(rep(0.05, 5)),
                   hyper_alpha=c(2.2, 2.8, 7, 22, 82),  #c(3.28, 2.02, 2.845, 2.32, 2.125),
                   hyper_beta=c(0.12, 0.36, 3, 21, 162))  #c(0.1824, 0.0102, 0.12, 0.0528, 0.028))
if(any(sapply(hyper_pars,length)<N_batch)){stop("Not enough hyper parameters for batch effect!")}
reduce_size <- as.logical(command_args[2])  # "FALSE"  # whether to take subsets to reduce batch sizes
N_reduce_size <- 70

# pipeline
iterations <- 5
plot_sim_ind <- NULL #seq(10,100,20)
use_ref_combat <- as.logical(command_args[3])  # whether to use ref combat to adjust test set against training set
set.seed(1)


####  Make training data smaller (reduce sample size)
if(reduce_size){
  reduced_ctrl <- sample(which(y_train==0)); reduced_case <- sample(which(y_train==1))#identical(sort(reduced_ctrl), which(y_train==0))
  reduced_ctrl <- reduced_ctrl[1:N_reduce_size]; reduced_case <- reduced_case[1:N_reduce_size]
  
  reduced_indices <- c(reduced_ctrl, reduced_case)
  train_expr <- train_expr[, reduced_indices]
  y_train <- y_train[reduced_indices]
}
print(table(y_train))


####  Run pipeline
ID <- 0
pred_mat_lst <- list()
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
    
  ## Simulate batch effect 
  sim_batch_res <- simBatch(dat=train_expr, condition=y_train, batches_ind=batches_ind, 
                            batch=batch, hyper_pars=hyper_pars)
  if(ID==1){print(round(do.call(cbind, sim_batch_res$batch_par), 3))}
  train_expr_batch <- sim_batch_res$new_dat
  #sapply(1:N_batch, function(i){mean(apply(train_expr_batch[, batches_ind[[i]]], 1, mean))})
  
  ## Visualize: PCA of training set with batch
  if(ID %in% plot_sim_ind){
    pca_res <- prcomp(t(train_expr_batch), center=TRUE, scale.=TRUE)
    pca_plt_obj <- data.frame(PC1=pca_res$x[, 1], PC2=pca_res$x[, 2],
                              Condition=y_train, Batch=as.factor(batch))
    png(sprintf("PCAplot_simBatch_ID%s.png", ID), width=6, height=5, units="in", res=300)
    ggplot(pca_plt_obj, aes(x=PC1, y=PC2)) +
      geom_point(aes(color=Batch, shape=Condition)) +
      scale_shape_manual(values=c(16, 17))+
      ggtitle("AEGIS-2 with simulated batch effect")
    dev.off()
    rm(pca_res, pca_plt_obj)
  }
  
  
  ####  Normalize datasets before training
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
  
  
  ####  Training
  ## Prediction from original training to test, without batch effect
  pred_base_res <- try(trainPipe(train_set=train_expr_norm, train_label=y_train, 
                                 test_set=test_expr_norm, lfit=learner_fit, use_ref_combat=use_ref_combat))
  if(class(pred_base_res)=="try-error"){ID <- ID - 1; next}
  
  ## Prediction from training WITH batch effect to test
  pred_batch_res <- try(trainPipe(train_set=train_expr_batch_whole_norm, train_label=y_train,
                                  test_set=test_expr_norm, lfit=learner_fit, use_ref_combat=use_ref_combat))
  if(class(pred_batch_res)=="try-error"){ID <- ID - 1; next}
  
  ##  Prediction from training after batch adjustment (Merged)  
  train_expr_combat <- try(ComBat(train_expr_batch, batch=batch, mod=model.matrix(~y_train)))
  if(class(train_expr_combat)=="try-error"){ID <- ID - 1; next}
  train_expr_combat_norm <- normalizeData(train_expr_combat)
  pred_combat_res <- try(trainPipe(train_set=train_expr_combat_norm, train_label=y_train, 
                                   test_set=test_expr_norm, lfit=learner_fit, use_ref_combat=use_ref_combat))
  if(class(pred_combat_res)=="try-error"){ID <- ID - 1; next}
  
  ## Obtain predictions from learner trained within each batch 
  pred_sgbatch_res <- try(lapply(1:N_batch, function(batch_id){
    trainPipe(train_set=train_expr_batch_norm[, batches_ind[[batch_id]]], train_label=y_sgbatch_train[[batch_id]], 
              test_set=test_expr_norm, lfit=learner_fit, use_ref_combat=use_ref_combat)
  }))
  if(class(pred_sgbatch_res)=="try-error"){ID <- ID - 1; next}
  names(pred_sgbatch_res) <- paste0("Batch", 1:N_batch)
  
  ##  Aggregate with different weights
  pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp$pred_tst_prob)})
  pred_mat <- do.call(cbind, pred_test_lst)
  
  # Avg: simple average
  pred_avg <- rowMeans(pred_mat)
  
  # n-Avg: sample-size-weighted average 
  pred_N_avg <- pred_mat %*% (as.matrix(sapply(batches_ind, length)) / length(y_train))
  
  # CS-Avg: replicability weights
  train_lst <- lapply(batches_ind, function(ind){train_expr_batch_norm[, ind]})
  cs_zmat <- try(CS_zmatrix(study_lst=train_lst, label_lst=y_sgbatch_train, 
                            lfit=learner_fit, perf_name=perf_measure_name))
  if(class(cs_zmat)=="try-error"){ID <- ID - 1; next}
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq
  
  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- try(Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train, lfit=learner_fit))
  if(class(reg_ssl_res)=="try-error"){ID <- ID - 1; next}
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), n_seq=sapply(batches_ind, length))
  pred_reg_a <- pred_mat %*% reg_a_beta
  
  # Reg-s:
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
  reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x  
  reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  pred_reg_s <- pred_mat %*% reg_s_beta
  
  
  ####  Evaluate performance 
  tst_scores <- c(list(NoBatch=pred_base_res$pred_tst_prob, Batch=pred_batch_res$pred_tst_prob),
                  pred_test_lst,   
                  list(ComBat=pred_combat_res$pred_tst_prob, 
                       Avg=pred_avg, n_Avg=pred_N_avg, CS_Avg=pred_cs_avg, 
                       Reg_a=pred_reg_a, Reg_s=pred_reg_s))
  
  # calculate performance
  perf_df <- as.data.frame(t(sapply(tst_scores, function(preds){
    if(perf_measure_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)} # avoid Inf in computing cross-entropy loss
    rocr_pred <- prediction(preds, as.numeric(as.character(y_test)))
    perf_measure_tmp <- performance(rocr_pred, perf_measure_name)  # mean cross entropy
    return(as.numeric(perf_measure_tmp@y.values))
  })))
  
  
  ####  Output results
  first_file <- !file.exists(sprintf('batchCSL_AEGIS_%s_%s_reduce%s_useref%s.csv', 
                                     perf_measure_name, learner_type, 
                                     ifelse(reduce_size, 'T', 'F'), ifelse(use_ref_combat, 'T', 'F')))
  write.table(perf_df, sprintf('batchCSL_AEGIS_%s_%s_reduce%s_useref%s.csv', 
                               perf_measure_name, learner_type, 
                               ifelse(reduce_size, 'T', 'F'), ifelse(use_ref_combat, 'T', 'F')), 
              append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
  pred_mat_lst[[ID]] <- pred_mat
}

save(pred_mat_lst, file=sprintf("testPredScores_%s_%s_Reduce%s_useRef%s.RData", 
                                perf_measure_name, learner_type, 
                                ifelse(reduce_size, 'T', 'F'), ifelse(use_ref_combat, 'T', 'F')))
