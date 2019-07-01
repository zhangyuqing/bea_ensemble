rm(list=ls()); demo <- TRUE
if(demo){
  setwd("~/Documents/MSBE/TB_realdata_Testswap/")
  script_dir <- "~/Dropbox/Work/MultiStudy_BatchEffect/TB"
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"
}else{
  #setwd("~/yuqingz/multistudy_batcheffect/TB_realdata_bootstrapTest_swap")
  #data_dir <- "../TB"
  setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/TB_realdata_bootstrapTest")
  data_dir <- "../TB_realdata_bootstrap"
  script_dir <- "."
}
sapply(c("glmnet", "SummarizedExperiment", "sva", "DESeq2", "ROCR", "ggplot2", "gridExtra", 
         "reshape2", "dplyr", "nnls"), require, character.only=TRUE)
rds_obj <- readRDS(file.path(data_dir, "combined.rds"))
source(file.path(script_dir, "helper_crossmod.R"))
set.seed(12345)



####  Parameters  ####
#command_args <- commandArgs(trailingOnly=TRUE)
norm_data <- TRUE #as.logical(command_args[1])  # whether to normalize datasets by features using z-score scaling
use_ref_combat <- FALSE #as.logical(command_args[8])  # whether to use ref combat to adjust test set against training set
match_preval <- TRUE #as.logical(command_args[1])

n_highvar_genes <- 1000  # number of highly variable genes to use in feature reduction
B <- 100  # bootstrap samples
data_name <- "logfpkm"
study_names <- c("Brazil_1", "India", "Africa")
learner_types <- c("lasso", "rf", "svm")  #, "nnet") 
perf_measures <- c("mxe", "auc")  #, "acc", "f")
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  #, "Accuracy", "F1 score")
names(perf_measures_names) <- perf_measures

# take brazil batch 1, india, africa studies
batch_filter <- rds_obj$SequencingBatch %in% study_names
group_filter <- rds_obj$Label %in% c("Non-progressor", "Active")
rds_obj <- rds_obj[, batch_filter & group_filter]

# take a subset of Africa to make the 3 data balanced
if(match_preval){
  africa_id <- which(rds_obj$SequencingBatch=="Africa")
  rm_id <- africa_id[1:round(length(africa_id) - table(rds_obj$Label[africa_id])['Active']/0.588)]
  rds_obj <- rds_obj[, -rm_id]
  new_africa_group <- rds_obj$Label[rds_obj$SequencingBatch=="Africa"]
  print(sprintf("Prevalence in Africa: %s", 
                round(table(new_africa_group)['Active'] / sum(table(new_africa_group)), 2)))
}

# remove genes with too many 0s or 0 variance in each study
matA <- assay(rds_obj[, rds_obj$SequencingBatch=="Africa"], data_name)
keepA <- rowVars(matA) > 0 & rowSums(matA!=0) > 2 
matB <- assay(rds_obj[, rds_obj$SequencingBatch=="Brazil_1"], data_name)
keepB <-  rowVars(matB) > 0 & rowSums(matB!=0) > 2 
matI <- assay(rds_obj[, rds_obj$SequencingBatch=="India"], data_name)
keepI <- rowVars(matI) > 0 & rowSums(matI!=0) > 2
rds_obj <- rds_obj[keepA & keepB & keepI, ]
rm(matA, matB, matI, keepA, keepB, keepI)


####  Start pipeline  ####
#s = "India"
start_time <- Sys.time()
for(s in study_names){
  ## Get into training & test set
  test_name <- s
  train_name <- setdiff(study_names, test_name)
  
  # training set
  rds_obj_train <- rds_obj[, (rds_obj$SequencingBatch %in% train_name)]
  dat <- assay(rds_obj_train, data_name)
  batch_org <- as.character(colData(rds_obj_train)$SequencingBatch)
  batch <- rep(1, length(batch_org)); batch[batch_org==train_name[2]] <- 2
  batch_names <- levels(factor(batch))
  group_org <- as.character(colData(rds_obj_train)$Label)
  group <- rep(0, length(group_org)); group[group_org=="Active"] <- 1
  covar <- rds_obj_train$Sex
  rm(rds_obj_train)
  
  # test
  rds_obj_test <- rds_obj[, (rds_obj$SequencingBatch==test_name) & (rds_obj$Label %in% c("Non-progressor", "Active"))]
  dat_test <- assay(rds_obj_test, data_name)
  group_test_org <- as.character(colData(rds_obj_test)$Label)
  group_test <- rep(0, length(group_test_org)); group_test[group_test_org=="Active"] <- 1
  covar_test <- rds_obj_test$Sex
  rm(rds_obj_test)

  
  ####  Preprocess
  ## use ComBat to adjust training data - estimate training set without batch effect
  dat_fs <- ComBat(dat, batch=batch, mod=model.matrix(~group+covar))
  
  ## feature reduction - select highly variable genes in training data
  genes_sel_names <- names(sort(rowVars(dat_fs), decreasing=TRUE))[1:n_highvar_genes]
  dat <- dat[genes_sel_names, ]
  dat_test <- dat_test[genes_sel_names, ]
  
  ## batch correction again after feature
  dat_combat <- ComBat(dat, batch=batch, mod=model.matrix(~group+covar))
  
  ## normalize features
  if(norm_data){
    print("Normalizing data.")
    dat_whole_norm <- normalizeData(dat)  # norm training set as a whole
    dat_batch_norm <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))  
    # norm training set within each batch
    for(k in batch_names){dat_batch_norm[, batch==k] <- normalizeData(dat[, batch==k])}
    dat_combat_norm <- normalizeData(dat_combat)  # norm combat adjusted data
  }else{
    print("Datasets are NOT normalized.")
    dat_whole_norm <- dat 
    dat_batch_norm <- dat
    dat_combat_norm <- dat_combat
  }
  train_lst <- lapply(batch_names, function(k){dat_batch_norm[, batch==k]})
  y_sgbatch_train <- lapply(batch_names, function(k){group[batch==k]})
  
  
  ####  Training 
  #l_type="lasso"
  unadj_mod_lst <- combat_mod_lst <- sgbatch_mod_lst <- list()
  cs_zmat_lst <- cs_weights_seq <- reg_ssl_res <- reg_a_beta <- reg_s_beta <- list()
  for(l_type in learner_types){
    learner_fit <- getPredFunctions(l_type)
    print(paste("Model:", l_type))
    
    ##  Training on original train set
    pred_unadj_res <- trainPipe(train_set=dat_whole_norm, train_label=group, test_set=NULL, 
                                lfit=learner_fit, use_ref_combat=use_ref_combat)
    unadj_mod_lst[[l_type]] <- pred_unadj_res$mod
    
    ##  Training on train set after batch adjustment (Merged)
    pred_combat_res <- trainPipe(train_set=dat_combat_norm, train_label=group, test_set=NULL, 
                                 lfit=learner_fit, use_ref_combat=use_ref_combat)
    combat_mod_lst[[l_type]] <- pred_combat_res$mod
      
    ##  Training within each batch from train set
    pred_sgbatch_res <- lapply(batch_names, function(k){
      trainPipe(train_set=dat_batch_norm[, batch==k], train_label=group[batch==k], test_set=NULL, 
                lfit=learner_fit, use_ref_combat=use_ref_combat)
    })
    names(pred_sgbatch_res) <- paste0("Batch", batch_names)
    sgbatch_mod_lst[[l_type]] <- lapply(pred_sgbatch_res, function(res){res$mod})

    
    ## Ensemble weights
    # cs
    cs_zmat_lst[[l_type]] <- CS_zmatrix(study_lst=train_lst, label_lst=y_sgbatch_train, 
                                        lfit=learner_fit, perf_name="mxe", use_ref_combat=use_ref_combat)
    cs_weights_seq[[l_type]] <- CS_weight(cs_zmat_lst[[l_type]])
    # reg-a
    reg_ssl_res[[l_type]] <- Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train, 
                                          lfit=learner_fit, use_ref_combat=use_ref_combat)
    reg_a_beta[[l_type]] <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res[[l_type]]$coef), 
                                         n_seq=table(batch)[batch_names])
    # reg-s
    stacked_pred <- do.call(rbind, reg_ssl_res[[l_type]]$pred)
    stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
    reg_s_beta[[l_type]] <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
    reg_s_beta[[l_type]] <- reg_s_beta[[l_type]] / sum(reg_s_beta[[l_type]])
  }
  rm(pred_unadj_res, pred_combat_res, pred_sgbatch_res)
  
  
  ####  Prediction & Ensemble
  b=1
  perf_df_lst <- tst_scores_modlst <- list()
  # save original test set
  dat_testOri <- dat_test; group_testOri <- group_test; covar_testOri <- covar_test  
  rm(dat_test, group_test, covar_test)  # save original data for bootstrap
  
  while(b<=B){
    print(sprintf("Test: %s; Bootstrap: %s", test_name, b))
    
    ## draw bootstrap sample
    boot_ind <- sample(1:ncol(dat_testOri), ncol(dat_testOri), replace=TRUE)
    dat_test <- dat_testOri[, boot_ind]
    group_test <- group_testOri[boot_ind]
    covar_test <- covar_testOri[boot_ind]
    # print(table(group_test))
    
    ## normalize
    if(norm_data){dat_test_norm <- normalizeData(dat_test)  # norm test set
    }else{dat_test_norm <- dat_test}
    if(any(is.na(dat_test_norm))){
      b <- b - 1
    }else{
      #l_type="lasso"
      for(l_type in learner_types){
        ##  Obtain prediction on each test set
        unadj_tst_prob <- predWrapper(unadj_mod_lst[[l_type]], dat_test_norm, l_type)
        combat_tst_prob <-  predWrapper(combat_mod_lst[[l_type]], dat_test_norm, l_type)
        pred_test_lst <- lapply(sgbatch_mod_lst[[l_type]], function(curr_mod){
          predWrapper(curr_mod, dat_test_norm, l_type)
        })
        
        ##  Aggregate with different weights
        pred_mat <- do.call(cbind, pred_test_lst)
        # Avg: simple average
        pred_avg <- rowMeans(pred_mat)
        # n-Avg: sample-size-weighted average
        pred_N_avg <- pred_mat %*% (as.matrix(table(batch)[batch_names]) / sum(table(batch)))
        # CS-Avg: replicability weights
        pred_cs_avg <- pred_mat %*% cs_weights_seq[[l_type]]
        # Reg-a: use each function to predict on one study, bind predictions and do regression
        pred_reg_a <- pred_mat %*% reg_a_beta[[l_type]]
        # Reg-s:
        pred_reg_s <- pred_mat %*% reg_s_beta[[l_type]]
        
        ##  Evaluate performance
        tst_scores <- c(list(Batch=unadj_tst_prob), pred_test_lst,
                        list(ComBat=combat_tst_prob, Avg=pred_avg, n_Avg=pred_N_avg, 
                             CS_Avg=pred_cs_avg, Reg_a=pred_reg_a, Reg_s=pred_reg_s))
        perf_df <- lapply(perf_measures, function(perf_name){
          as.data.frame(t(sapply(tst_scores, function(preds){
            if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  # avoid Inf in computing cross-entropy loss
            rocr_pred <- prediction(preds, group_test)
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
        
        ##  Cache results
        perf_df_lst[[l_type]] <- perf_df
        tst_scores_modlst[[l_type]] <- tst_scores
      }
      
      
      ####  Ensemble across models
      print(paste("Ensemble across models."))
      
      preds_crossmod <- lapply(tst_scores_modlst, function(x){
        do.call(cbind, x[paste0("Batch", batch_names)])
      })
      pred_mat_crossmod <- do.call(cbind, preds_crossmod)
      pred_mat_crossmod_reg <- lapply(batch_names, function(i){
        do.call(cbind, lapply(preds_crossmod, function(pred)pred[,paste0("Batch", i)]))
      })
      pred_mat_crossmod_reg <- do.call(cbind, pred_mat_crossmod_reg)
      
      # Avg
      cm_avg <- rowMeans(pred_mat_crossmod)
      # n-Avg
      navg_weights <- rep((as.matrix(table(batch)[batch_names]) / length(group)), length(learner_types))
      cm_N_avg <- pred_mat_crossmod %*% (navg_weights / sum(navg_weights))
      # CS-Avg
      cm_cs_weights_seq <- CS_weight_crossmod(cs_zmat_lst)
      cm_cs_avg <- pred_mat_crossmod %*% cm_cs_weights_seq
      # Reg-a
      cm_reg_ssl_res <- crossmod_Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train,
                                              learner_lst=learner_types, use_ref_combat=use_ref_combat)
      cm_reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, cm_reg_ssl_res$coef), 
                                    n_seq=table(batch)[batch_names])
      cm_reg_a <- pred_mat_crossmod_reg %*% cm_reg_a_beta
      # Reg-s
      cm_stacked_pred <- do.call(rbind, cm_reg_ssl_res$pred)
      cm_stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
      cm_reg_s_beta <- nnls(A=cm_stacked_pred, b=as.numeric(cm_stacked_label))$x
      cm_reg_s <- pred_mat_crossmod_reg %*% (cm_reg_s_beta / sum(cm_reg_s_beta))
      
      ## calculate performance
      tst_cm_scores <- list(Avg=cm_avg, n_Avg=cm_N_avg, CS_Avg=cm_cs_avg, Reg_a=cm_reg_a, Reg_s=cm_reg_s)
      perf_crossmod_df <- lapply(perf_measures, function(perf_name){
        as.data.frame(t(sapply(tst_cm_scores, function(preds){
          if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  # avoid Inf in computing cross-entropy loss
          rocr_pred <- prediction(preds, as.numeric(as.character(group_test)))
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
      
      perf_df_lst[["crossmod"]] <- cbind(perf_df_lst$rf[,1:4], perf_crossmod_df)
      tst_scores_modlst[["crossmod"]] <- tst_cm_scores
      
      for(i in 1:length(perf_measures)){
        summary_df <- melt(lapply(perf_df_lst, function(perf_res){perf_res[perf_measures[i], -c(2:3)]}))
        summary_df$iteration <- b
        
        ## write out performances
        first_file <- !file.exists(sprintf('test%s_%s.csv', test_name, perf_measures[i]))
        write.table(summary_df, sprintf('test%s_%s.csv', test_name, perf_measures[i]),
                    append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
      }
      
      rm(boot_ind, dat_test, group_test, covar_test)
    }
    b <- b + 1
  }
}
end_time <- Sys.time()
print(end_time - start_time)
