rm(list=ls()); demo <- FALSE
if(demo){
  setwd("~/Dropbox/Work/MultiStudy_BatchEffect/TB/")
  rds_obj <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB/combined.rds")
}else{
  setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/TB_crossmod_reproduce")
  rds_obj <- readRDS("../TB_realdata_bootstrap/combined.rds")
  #setwd("~/yuqingz/multistudy_batcheffect/TB_crossmod_reproduce")
  #rds_obj <- readRDS("../TB/combined.rds")
}
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "limma", "nnls",  "glmnet", "rpart", 
             "genefilter", "nnet", "e1071", "RcppArmadillo", "foreach", "parallel", "doParallel",  
             "ranger", "scales"), require, character.only=TRUE))
source("helper_crossmod.R")
#set.seed(123)


####  Load data
data_name <- "logfpkm"   #"logtpm"

tb_africa <- rds_obj[, rds_obj$SequencingBatch=="Africa"]
tb_africa_sub <- tb_africa[, tb_africa$Label %in% c("Non-progressor", "Progressor")]
tb_g6 <- rds_obj[, rds_obj$SequencingBatch=="G6"]

train_expr <- assays(tb_africa_sub)[[data_name]]
y_train <- revalue(as.factor(as.character(tb_africa_sub$Label)), 
                   c("Non-progressor"="0", "Progressor"="1"))
test_expr <- assays(tb_g6)[[data_name]]
y_test <- revalue(as.factor(as.character(tb_g6$Label)), 
                  c("Non-progressor"="0", "Progressor"="1"))
identical(rownames(train_expr), rownames(test_expr))
rm(tb_africa, tb_africa_sub, tb_g6)

#select 1000 genes with largest variance in training set (Africa)
var_trn <- rowVars(train_expr)
genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:1000]]
train_expr <- train_expr[genes_sel, ]
test_expr <- test_expr[genes_sel, ]


####  Parameters 
command_args <- commandArgs(trailingOnly=TRUE)  
#command_args <- c("20", "3", "4")
if(length(command_args)!=3){stop("Not enough input parameters!")}

## Prediction model
learner_types <- c("lasso", "rf", "svm")   #c("lasso", "rf", "nnet", "svm") 

## Degree of batch effect (strength of signal)
N_batch <- 3
N_sample_size <- as.numeric(command_args[1])   # number of samples per batch (case + control)
max_batch_mean <- as.numeric(command_args[2]) 
# median 8, range(4, 12), recommended values 0-3 (additive)
max_batch_var <- as.numeric(command_args[3]) 
# median 0.1, range(0.01, 1), recommended values 1-10 (multiplicative)
hyper_pars <- list(hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  #rep(0,5),
                   hyper_sd=sqrt(rep(0.01, N_batch)),
                   hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                     v=rep(0.01, N_batch))$alpha,  #c(3.28, 2.02, 2.845, 2.32, 2.125),
                   hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                    v=rep(0.01, N_batch))$beta)  #c(0.1824, 0.0102, 0.12, 0.0528, 0.028))
cat("\nBatch changes\n");
print(hyper_pars$hyper_mu);
print(hyper_pars$hyper_beta/(hyper_pars$hyper_alpha-1))
#print((hyper_pars$hyper_beta^2)/((hyper_pars$hyper_alpha-1)^2*(hyper_pars$hyper_alpha-2)))

## Pipeline
iterations <- 100
norm_data <- TRUE #as.logical(command_args[7])  # whether to normalize datasets by features
use_ref_combat <- FALSE #as.logical(command_args[8])  
# whether to use ref combat to adjust test set against training set
exp_name <- sprintf('batchN%s_m%s_v%s', N_sample_size, 
                    gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T))  
perf_measures <- c("mxe", "auc")    #c("mxe", "auc", "acc", "f")


####  Run pipeline
#ID = 1; l_type = learner_types[1]
#start_time <- Sys.time()
for(ID in 1:iterations){
  tst_scores_modlst <- cs_zmat_lst <- list()  # list of results for each model
  
  ## Subset training set in batches
  batches_ind <- subsetBatch(condition=y_train, N_sample_size=N_sample_size, N_batch=N_batch)
  y_sgbatch_train <- lapply(1:N_batch, function(k){y_train[batches_ind[[k]]]})
  batch <- rep(0, ncol(train_expr)); for(i in 1:N_batch){batch[batches_ind[[i]]] <- i}
  curr_train_expr <- train_expr[, do.call(c, batches_ind)]
  curr_y_train <- y_train[do.call(c, batches_ind)]
  batch <- batch[do.call(c, batches_ind)]  
  batches_ind <- lapply(1:N_batch, function(i){which(batch==i)})
  
  ## Remove genes with only 0 values in any batch in current training set
  g_keep <- lapply(1:N_batch, function(j){which(apply(curr_train_expr[, batch==j], 1, function(x){!all(x==0)}))})
  g_keep <- Reduce(intersect, g_keep)  
  curr_train_expr <- curr_train_expr[g_keep, ]
  curr_test_expr <- test_expr[g_keep, ]
  
  ## Simulate batch effect 
  sim_batch_res <- simBatch(dat=curr_train_expr, condition=curr_y_train, batches_ind=batches_ind, 
                            batch=batch, hyper_pars=hyper_pars)
  #if(ID==1){print(round(do.call(cbind, sim_batch_res$batch_par), 3))}
  train_expr_batch <- sim_batch_res$new_dat
  #sapply(1:N_batch, function(i){mean(apply(train_expr_batch[, batches_ind[[i]]], 1, mean))})
  
  ## Normalize datasets before training
  if(norm_data){
    if(ID==1){print("Normalizing data.")}
    train_expr_norm <- normalizeData(curr_train_expr)
    test_expr_norm <- normalizeData(curr_test_expr)
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
    train_expr_norm <- curr_train_expr;  test_expr_norm <- curr_test_expr
    train_expr_batch_whole_norm <- train_expr_batch_norm <- train_expr_batch
  }
  train_lst <- lapply(batches_ind, function(ind){train_expr_batch_norm[, ind]})
  
  
  ####  Training with single model
  for(l_type in learner_types){
    learner_fit <- getPredFunctions(l_type)
    print(sprintf("Simulation: %s, Model: %s", ID, l_type))
    
    ## Prediction from original training to test, without batch effect
    pred_base_res <- try(trainPipe(train_set=train_expr_norm, train_label=curr_y_train, 
                                   test_set=test_expr_norm, 
                                   lfit=learner_fit, use_ref_combat=use_ref_combat))
    if(class(pred_base_res)=="try-error"){ID <- ID - 1; break; next}
    
    ## Prediction from training WITH batch effect to test
    pred_batch_res <- try(trainPipe(train_set=train_expr_batch_whole_norm, train_label=curr_y_train, 
                                    test_set=test_expr_norm, 
                                    lfit=learner_fit, use_ref_combat=use_ref_combat))
    if(class(pred_batch_res)=="try-error"){ID <- ID - 1; break; next}
    
    ##  Prediction from training after batch adjustment (Merged)
    train_expr_combat <- try(ComBat(train_expr_batch, batch=batch, mod=model.matrix(~curr_y_train)))
    if(class(train_expr_combat)=="try-error"){ID <- ID - 1; break; next}
    if(norm_data){
      train_expr_combat_norm <- normalizeData(train_expr_combat)
    }else{
      train_expr_combat_norm <- train_expr_combat
    }
    pred_combat_res <- try(trainPipe(train_set=train_expr_combat_norm, train_label=curr_y_train, 
                                     test_set=test_expr_norm, 
                                     lfit=learner_fit, use_ref_combat=use_ref_combat))
    if(class(pred_combat_res)=="try-error"){ID <- ID - 1; break; next}
    
    ## Obtain predictions from learner trained within each batch
    pred_sgbatch_res <- try(lapply(1:N_batch, function(batch_id){
      trainPipe(train_set=train_expr_batch_norm[, batches_ind[[batch_id]]], 
                train_label=y_sgbatch_train[[batch_id]],
                test_set=test_expr_norm, lfit=learner_fit, use_ref_combat=use_ref_combat)
    }))
    if(class(pred_sgbatch_res)=="try-error"){ID <- ID - 1; break; next}
    names(pred_sgbatch_res) <- paste0("Batch", 1:N_batch)
    
    ##  Aggregate with different weights
    pred_test_lst <- lapply(pred_sgbatch_res, function(tmp){return(tmp$pred_tst_prob)})
    pred_mat <- do.call(cbind, pred_test_lst)
    
    # Avg: simple average
    pred_avg <- rowMeans(pred_mat)
    # n-Avg: sample-size-weighted average
    pred_N_avg <- pred_mat %*% (as.matrix(sapply(batches_ind, length)) / length(curr_y_train))
    
    # CS-Avg: replicability weights
    cs_zmat <- try(CS_zmatrix(study_lst=train_lst, label_lst=y_sgbatch_train, lfit=learner_fit, 
                              perf_name="mxe", use_ref_combat=use_ref_combat))
    if(class(cs_zmat)=="try-error"){ID <- ID - 1; break; next}
    cs_weights_seq <- CS_weight(cs_zmat)
    pred_cs_avg <- pred_mat %*% cs_weights_seq
    
    # Reg-a: use each function to predict on one study, bind predictions and do regression
    reg_ssl_res <- try(Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train, 
                                    lfit=learner_fit, use_ref_combat=use_ref_combat))
    if(class(reg_ssl_res)=="try-error"){ID <- ID - 1; break; next}
    reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef), 
                               n_seq=sapply(batches_ind, length))
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
                    list(ComBat=pred_combat_res$pred_tst_prob, Avg=pred_avg, n_Avg=pred_N_avg, 
                         CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s))
    
    perf_df <- lapply(perf_measures, function(perf_name){
      as.data.frame(t(sapply(tst_scores, function(preds){
        if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  
        # avoid Inf in computing cross-entropy loss
        rocr_pred <- prediction(preds, as.numeric(as.character(y_test)))
        if(perf_name %in% c("acc", "f")){
          curr_perf <- performance(rocr_pred, perf_name)  
          return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
        }else{
          curr_perf <- performance(rocr_pred, perf_name)
          return(as.numeric(curr_perf@y.values))
        }
      })))
    })
    names(perf_df) <- perf_measures
    perf_df <- do.call(rbind, perf_df)
    
    
    ####  Output results
    first_file <- !file.exists(sprintf('results/%s_mxe_%s.csv', l_type, exp_name))
    wo <- sapply(perf_measures, function(perf_name){
      write.table(perf_df[perf_name, ], sprintf('results/%s_%s_%s.csv', l_type, perf_name, exp_name),
                  append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
    })
    tst_scores_modlst[[l_type]] <- tst_scores
    cs_zmat_lst[[l_type]] <- cs_zmat
  }
  
  
  ####  Ensemble across models
  ## Aggregate with different weights
  preds_crossmod <- lapply(tst_scores_modlst, function(x){do.call(cbind, x[paste0("Batch", 1:N_batch)])})
  pred_mat_crossmod <- do.call(cbind, preds_crossmod)
  pred_mat_crossmod_reg <- lapply(1:N_batch, function(i){
    do.call(cbind, lapply(preds_crossmod, function(pred)pred[,i]))
  })
  pred_mat_crossmod_reg <- do.call(cbind, pred_mat_crossmod_reg)
  # Avg
  cm_avg <- rowMeans(pred_mat_crossmod)
  # n-Avg
  navg_weights <- rep((as.matrix(sapply(batches_ind, length)) / length(curr_y_train)), length(learner_types))
  cm_N_avg <- pred_mat_crossmod %*% (navg_weights / sum(navg_weights))
  # CS-Avg
  cm_cs_weights_seq <- CS_weight_crossmod(cs_zmat_lst)
  cm_cs_avg <- pred_mat_crossmod %*% cm_cs_weights_seq
  # Reg-a
  cm_reg_ssl_res <- try(crossmod_Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train,
                                              learner_lst=learner_types, use_ref_combat=use_ref_combat))
  if(class(cm_reg_ssl_res)=="try-error"){ID <- ID - 1; next}
  cm_reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, cm_reg_ssl_res$coef), n_seq=sapply(batches_ind, length))
  cm_reg_a <- pred_mat_crossmod_reg %*% cm_reg_a_beta
  # Reg-s
  cm_stacked_pred <- do.call(rbind, cm_reg_ssl_res$pred)
  cm_stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
  cm_reg_s_beta <- nnls(A=cm_stacked_pred, b=as.numeric(cm_stacked_label))$x
  cm_reg_s <- pred_mat_crossmod_reg %*% (cm_reg_s_beta / sum(cm_reg_s_beta))
  
  
  ## Calculate performance
  tst_cm_scores <- list(Avg=cm_avg, n_Avg=cm_N_avg, CS_Avg=cm_cs_avg, Reg_a=cm_reg_a, Reg_s=cm_reg_s)
  perf_crossmod_df <- lapply(perf_measures, function(perf_name){
    as.data.frame(t(sapply(tst_cm_scores, function(preds){
      if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  
      # avoid Inf in computing cross-entropy loss
      rocr_pred <- prediction(preds, as.numeric(as.character(y_test)))
      if(perf_name %in% c("acc", "f")){
        curr_perf <- performance(rocr_pred, perf_name)  
        return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
      }else{
        curr_perf <- performance(rocr_pred, perf_name)
        return(as.numeric(curr_perf@y.values))
      }
    })))
  })
  names(perf_crossmod_df) <- perf_measures
  perf_crossmod_df <- do.call(rbind, perf_crossmod_df)

  ## output ensemble results across models
  first_file <- !file.exists(sprintf('results/crossmod_mxe_%s.csv', exp_name))
  wocm <- sapply(perf_measures, function(perf_name){
    write.table(perf_crossmod_df[perf_name, ], sprintf('results/crossmod_%s_%s.csv', perf_name, exp_name),
                append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
  })
}
#end_time <- Sys.time()
#print(end_time - start_time)
