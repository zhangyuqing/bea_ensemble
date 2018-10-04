rm(list=ls())
setwd("~/yuqingz/multistudy_batcheffect/diff_expr/")
#setwd("C:/Users/zhang/Dropbox/Work/MultiStudy_BatchEffect/differential_expression")
sapply(c("sva", "MCMCpack", "BatchQC", "ROCR", "ggplot2", "limma", "nnls"), require, character.only=TRUE)
source("helper_DE.R")
load("sim_data.RData")  #rm(test_expr, y_test)


####  Parameters 
#command_args <- commandArgs(trailingOnly=TRUE)
# batch
N_batch <- 5 #as.numeric(command_args[1])
hyper_pars <- list(hyper_mu=c(-0.2, -0.1, 0, 0.1, 0.2),  #rep(0,5),
                   hyper_sd=sqrt(rep(0.001, 5)),
                   hyper_alpha=c(4.5, 12, 24.5, 42, 64.5),  #c(2.2, 2.8, 7, 22, 82),  
                   hyper_beta=c(0.175, 1.1, 3.525, 8.2, 15.875))  #c(0.12, 0.36, 3, 21, 162))  
if(any(sapply(hyper_pars,length)<N_batch)){stop("Not enough hyper parameters for batch effect!")}

# pipeline
iterations <- 100  #5
exp_name <- "avg"  #command_args[1]
set.seed(1)



####  Run pipeline
ID <- 0
while(ID < iterations){
  ID <- ID + 1
  print(paste("Simulation:", ID))
  
  
  ####  Split training set (AEGIS-2) in batches
  batches_ind <- splitBatch(condition=y_train, N_batch=N_batch)
  #tmp=do.call(rbind,lapply(1:N_batch,function(i){table(y_train[batches_ind[[i]]])}));print(tmp);print(apply(tmp,2,sum))
  y_sgbatch_train <- lapply(1:N_batch, function(k){y_train[batches_ind[[k]]]})
  if(any(sapply(y_sgbatch_train,table)<5)){ID <- ID - 1; next}
  batch <- rep(0, ncol(train_expr))
  for(i in 1:N_batch){batch[batches_ind[[i]]] <- i}
  
  
  ####  Simulate batch effect 
  sim_batch_res <- simBatch(dat=train_expr, condition=y_train, batches_ind=batches_ind, 
                            batch=batch, hyper_pars=hyper_pars)
  if(ID==1){print(round(do.call(cbind, sim_batch_res$batch_par), 3))}
  train_expr_batch <- sim_batch_res$new_dat
  
  
  ####  ComBat + DE on whole data
  mod_whole <- model.matrix(~y_train)
  train_expr_combat <- try(ComBat(train_expr_batch, batch=batch, mod=mod_whole))
  if(class(train_expr_combat)=="try-error"){ID <- ID - 1; next}
  
  res_combat <- t(apply(train_expr_combat, 1, function(x, group){
    x <- scale(x, center=TRUE, scale=TRUE)
    lm_f <- lm(x ~ group)
    return(summary(lm_f)$coefficients[2, ])
  }, group=y_train))
  
  res_combat_sorted <- res_combat[order(res_combat[, 4]), ]
  combat_DE_lst <- rownames(res_combat_sorted)[1:N_DE]
  combat_DE_stats <- calculate_stats(called_lst=combat_DE_lst, truth=ground_truth, N_genes, N_DE)
  
    
  ####  Ensemble-based DE
  ## Obtain batch-wise coefficients and standard error
  beta_lst <- std_lst <- list()
  for(b in 1:N_batch){
    batch_dat <- train_expr_batch[, batch==b]
    batch_y <- y_sgbatch_train[[b]]
    beta_batch <- std_batch <- c()
    for(ii in 1:nrow(batch_dat)){
      x <- scale(batch_dat[ii, ], center=TRUE, scale=TRUE)
      lm_batch <- lm(x ~ batch_y)
      beta_batch <- c(beta_batch, summary(lm_batch)$coefficients[2, "Estimate"])
      std_batch <- c(std_batch, summary(lm_batch)$coefficients[2, "Std. Error"])
    }
    beta_lst[[b]] <- beta_batch; std_lst[[b]] <- std_batch
  }
  beta_mat <- do.call(cbind, beta_lst); std_mat <- do.call(cbind, std_lst)

  
  ## Simple average
  weights_avg <- rep(1/N_batch, N_batch)
  avg_pvalues <- ensmbl_DE_avg(beta_mat, std_mat, weights_seq=weights_avg)
  avg_DE_lst <- paste0("gene", order(avg_pvalues)[1:N_DE])
  ensmbl_avg_stats <- calculate_stats(called_lst=avg_DE_lst, truth=ground_truth, N_genes, N_DE)
  
  
  ## Batch size weighted average
  # weights_n_avg <- table(batch) / sum(table(batch))
  # n_avg_pvalues <- ensmbl_DE_avg(beta_mat, std_mat, weights_seq=weights_n_avg)
  # n_avg_DE_lst <- paste0("gene", order(n_avg_pvalues)[1:N_DE])
  # ensmbl_n_avg_stats <- calculate_stats(called_lst=n_avg_DE_lst, truth=ground_truth, N_genes, N_DE)
  
  
  ####  Output tp and fp statistics
  res_out <- c(combat_DE_stats, ensmbl_avg_stats)  #, ensmbl_n_avg_stats)
  names(res_out) <- c("combat_tpr", "combat_fpr", "avg_tpr", "avg_fpr")  #, "N-avg_tpr", "N-avg_fpr")
    
  first.file <- !file.exists(sprintf("simDE_%s.csv", exp_name))
  write.table(as.data.frame(t(res_out)), sprintf("simDE_%s.csv", exp_name), 
              append=!first.file, col.names=first.file, row.names=FALSE, sep=",")
}
