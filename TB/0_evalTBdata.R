rm(list=ls())
setwd("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB/")
sapply(c("SummarizedExperiment", "glmnet", "plyr", "MLmetrics"), require, character.only=TRUE)
data_name <- "logfpkm" #"logtpm"
source("~/Dropbox/Work/MultiStudy_BatchEffect/TB/helper_new.R")

rds_obj <- readRDS("combined.rds")
G_seq <- c(seq(50, 450, 50), seq(500, 5000, 500)) #c(300, seq(500, 5000, 500))  

## organize train and test set
tb_africa <- rds_obj[, rds_obj$SequencingBatch=="Africa"]
tb_africa_sub <- tb_africa[, tb_africa$Label %in% c("Non-progressor", "Progressor")]
tb_g6 <- rds_obj[, rds_obj$SequencingBatch=="G6"]

trn_set <- assays(tb_africa_sub)[[data_name]]
y_trn <- revalue(as.factor(as.character(tb_africa_sub$Label)), c("Non-progressor"="0", "Progressor"="1"))
tst_set <- assays(tb_g6)[[data_name]]
y_tst <- revalue(as.factor(as.character(tb_g6$Label)), c("Non-progressor"="0", "Progressor"="1"))

## remove genes with constant values across samples (no variation)
# keep_trn <- rowVars(trn_set)!=0
# keep_tst <- rowVars(tst_set)!=0
# trn_set <- trn_set[keep_trn & keep_tst, ]
# tst_set <- tst_set[keep_trn & keep_tst, ]

perf_stats <- list()
ii <- 1
#G <- 500
for(G in G_seq){
  print(paste("G =", G))
  
  ## select top G highly variable genes
  var_trn <- rowVars(trn_set); #var_tst <- rowVars(tst_set)
  genes_sel <- rownames(trn_set)[order(var_trn, decreasing=TRUE)[1:G]]
  curr_trn_set <- trn_set[genes_sel, ]; curr_tst_set <- tst_set[genes_sel, ]
  
  ## normalize data
  # curr_trn_set <- t(apply(curr_trn_set, 1, scale, center=TRUE, scale=TRUE))
  # curr_tst_set <- t(apply(curr_tst_set, 1, scale, center=TRUE, scale=TRUE))
  
  ## train & predict
  pred_res <- predLasso_pp(curr_trn_set, curr_tst_set, y_trn)
    
  ## calculate performances
  perf_df <- matrix(NA, nrow=2, ncol=4, dimnames=list(c("trn", "tst"), c("LogLoss", "AUC", "F1 score", "Accuracy")))
  perf_df["trn", "LogLoss"] <- MLmetrics::LogLoss(pred_res$pred_trn_prob, as.numeric(as.character(y_trn)))
  perf_df["tst", "LogLoss"] <- MLmetrics::LogLoss(pred_res$pred_tst_prob, as.numeric(as.character(y_tst)))
  perf_df["trn", "AUC"] <- MLmetrics::AUC(pred_res$pred_trn_prob, as.numeric(as.character(y_trn)))
  perf_df["tst", "AUC"] <- MLmetrics::AUC(pred_res$pred_tst_prob, as.numeric(as.character(y_tst)))
  perf_df["trn", "F1 score"] <- MLmetrics::F1_Score(as.numeric(as.character(y_trn)), pred_res$pred_trn_class)
  perf_df["tst", "F1 score"] <- MLmetrics::F1_Score(as.numeric(as.character(y_tst)), pred_res$pred_tst_class)
  perf_df["trn", "Accuracy"] <- MLmetrics::Accuracy(as.numeric(as.character(y_trn)), pred_res$pred_trn_class)
  perf_df["tst", "Accuracy"] <- MLmetrics::Accuracy(as.numeric(as.character(y_tst)), pred_res$pred_tst_class)
  
  perf_stats[[ii]] <- perf_df
  ii <- ii + 1
}
names(perf_stats) <- G_seq

perf_stats_trn <- lapply(perf_stats, function(ps){return(ps["trn", ])})
perf_stats_trn <- data.frame(G=G_seq, do.call(rbind, perf_stats_trn))
perf_stats_tst <- lapply(perf_stats, function(ps){return(ps["tst", ])})
perf_stats_tst <- data.frame(G=G_seq, do.call(rbind, perf_stats_tst))
print(perf_stats_trn)
print(perf_stats_tst)
G_seq[which.max(perf_stats_tst$F1.score)]

## 
