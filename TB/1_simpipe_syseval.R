rm(list=ls())
demo <- FALSE
if(demo){
  setwd("~/Dropbox/Work/MultiStudy_BatchEffect/TB/")
  rds_obj <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB/combined.rds")
}else{
  setwd("~/yuqingz/multistudy_batcheffect/TB")
  rds_obj <- readRDS("combined.rds")
}
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "limma", "nnls", 
             "glmnet", "rpart", "genefilter", "nnet", "e1071", "RcppArmadillo", "foreach", "parallel", 
             "doParallel", "MLmetrics", "ranger"), require, character.only=TRUE))
source("helper_new.R")
set.seed(123)


####  Load data
data_name <- "logfpkm"   #"logtpm"

tb_africa <- rds_obj[, rds_obj$SequencingBatch=="Africa"]
tb_africa_sub <- tb_africa[, tb_africa$Label %in% c("Non-progressor", "Progressor")]
tb_g6 <- rds_obj[, rds_obj$SequencingBatch=="G6"]

train_expr <- assays(tb_africa_sub)[[data_name]]
y_train <- revalue(as.factor(as.character(tb_africa_sub$Label)), c("Non-progressor"="0", "Progressor"="1"))
test_expr <- assays(tb_g6)[[data_name]]
y_test <- revalue(as.factor(as.character(tb_g6$Label)), c("Non-progressor"="0", "Progressor"="1"))

rm(tb_africa, tb_africa_sub, tb_g6)

# select top G highly variable genes
# G_seq <- c(seq(50, 450, 50), seq(500, 5000, 500)) #c(300, seq(500, 5000, 500))  
# var_trn <- rowVars(train_expr)
# logloss_fs <- sapply(G_seq, function(G){
#   genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:G]]
#   curr_train_expr <- train_expr[genes_sel, ]; curr_test_expr <- test_expr[genes_sel, ]
#   pred_res <- predLasso_pp(curr_train_expr, curr_test_expr, y_train)
#   MLmetrics::LogLoss(pred_res$pred_tst_prob, as.numeric(as.character(y_test)))
# })
# G_sel <- G_seq[which.min(logloss_fs)]
# genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:G_sel]]
# train_expr <- train_expr[genes_sel, ]
# test_expr <- test_expr[genes_sel, ]

var_trn <- rowVars(train_expr)
genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:300]]
train_expr <- train_expr[genes_sel, ]
test_expr <- test_expr[genes_sel, ]


####  Parameters 
command_args <- commandArgs(trailingOnly=TRUE)  
#command_args <- c("rf", "5", "4", "4", "sig")
if(length(command_args)!=5){stop("Not enough input parameters!")}

## Prediction model
learner_type <- command_args[1]   # "lasso" 
learner_fit <- getPredFunctions(learner_type)

## Degree of batch effect (strength of signal)
N_batch <- as.numeric(command_args[2])  # 5
max_batch_mean <- as.numeric(command_args[3]) 
# median 8, range(4, 12), recommended values 0-3 (additive)
max_batch_var <- as.numeric(command_args[4]) 
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

## Pipeline
iterations <- 100
norm_data <- TRUE #as.logical(command_args[7])  # whether to normalize datasets by features
use_ref_combat <- TRUE #as.logical(command_args[8])  # whether to use ref combat to adjust test set against training set
test_item <- command_args[5]  # factor being tested, c("sig", "size", "nbatch")
exp_name <- sprintf('%s_%s_batchN%s_m%s_v%s', #_size%s', #_useref%s', 
                    test_item, learner_type, N_batch, 
                    gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T))#,
                    #ifelse(reduce_size, N_reduce_size, 'Org'))  #, ifelse(use_ref_combat, 'T', 'F'))


####  Run pipeline
ID <- 0
pred_mat_lst <- hyper_pars_lst <- perfstats_lst <- list()
while(ID < iterations){
  ID <- ID + 1
  print(paste("Simulation:", ID))
  
  ####  Spike in batch effect
  ## Split training set in batches
  batches_ind <- splitBatch(condition=y_train, N_batch=N_batch)
  #tmp=do.call(rbind,lapply(1:N_batch,function(i){table(y_train[batches_ind[[i]]])}));print(tmp);print(apply(tmp,2,sum))
  y_sgbatch_train <- lapply(1:N_batch, function(k){y_train[batches_ind[[k]]]})
  if(any(sapply(y_sgbatch_train,table)<9)){ID <- ID - 1; next}
  batch <- rep(0, ncol(train_expr))
  for(i in 1:N_batch){batch[batches_ind[[i]]] <- i}
  
  ## Randomly shuffle batch parameters to different combinations of mean-var batch effect
  # shuffled_order <- sample(1:N_batch, N_batch, replace=F)
  # hyper_pars_shuffled <- c(hyper_pars[1:2], lapply(hyper_pars[3:4], function(x){x[shuffled_order]}))
  # if(!identical(ab2mv(a=hyper_pars_shuffled$hyper_alpha, b=hyper_pars_shuffled$hyper_beta)$var, 
  #               rep(0.01, N_batch))){stop("Error in shuffling hyper pars!!")}
  hyper_pars_shuffled <- hyper_pars
    
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
  
  ## Prediction from training WITH batch effect to test
  pred_batch_res <- try(trainPipe(train_set=train_expr_batch_whole_norm, train_label=y_train,
                                  test_set=test_expr_norm, lfit=learner_fit,
                                  use_ref_combat=use_ref_combat))
  if(class(pred_batch_res)=="try-error"){ID <- ID - 1; next}

  ##  Prediction from training after batch adjustment (Merged)
  train_expr_combat <- try(ComBat(train_expr_batch, batch=batch, mod=model.matrix(~y_train)))
  if(class(train_expr_combat)=="try-error"){ID <- ID - 1; next}
  if(norm_data){
    train_expr_combat_norm <- normalizeData(train_expr_combat)
  }else{
    train_expr_combat_norm <- train_expr_combat
  }
  pred_combat_res <- try(trainPipe(train_set=train_expr_combat_norm, train_label=y_train,
                                   test_set=test_expr_norm, lfit=learner_fit,
                                   use_ref_combat=use_ref_combat))
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
                            lfit=learner_fit, perf_name="mxe",
                            use_ref_combat=use_ref_combat))
  if(class(cs_zmat)=="try-error"){ID <- ID - 1; next}
  cs_weights_seq <- CS_weight(cs_zmat)
  pred_cs_avg <- pred_mat %*% cs_weights_seq

  # Reg-a: use each function to predict on one study, bind predictions and do regression
  reg_ssl_res <- try(Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train,
                                  lfit=learner_fit, use_ref_combat=use_ref_combat))
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
  perf_df <- lapply(c("mxe", "auc", "acc", "f"), function(perf_name){
    as.data.frame(t(sapply(tst_scores, function(preds){
      if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  # avoid Inf in computing cross-entropy loss
      rocr_pred <- prediction(preds, as.numeric(as.character(y_test)))
      if(perf_name %in% c("acc", "f")){
        curr_perf<- performance(rocr_pred, perf_name)  # mean cross entropy
        return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
      }else{
        curr_perf<- performance(rocr_pred, perf_name)
        return(as.numeric(curr_perf@y.values))
      }
    })))
  })
  names(perf_df) <- c("mxe", "auc", "acc", "f")
  perf_df <- do.call(rbind, perf_df)

  ####  Output results
  first_file <- !file.exists(sprintf('results/mxe_%s.csv', exp_name))
  wo <- sapply(c("mxe", "auc", "acc", "f"), function(perf_name){
    write.table(perf_df[perf_name, ], sprintf('results/%s_%s.csv', perf_name, exp_name),
                append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
  })
}

#perfstats_lst <- do.call(rbind, perfstats_lst)
#print(colMeans(perfstats_lst, na.rm=T))
#save(pred_mat_lst, hyper_pars_lst, perfstats_lst, file=sprintf("results/predScores_%s.RData", exp_name))
