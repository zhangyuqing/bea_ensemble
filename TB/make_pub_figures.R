########  Figure 1: PCA of training data  ########
rm(list=ls())
setwd("~/Documents/MSBE/pub_figures/")
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "plyr", "ggpubr", "SummarizedExperiment", "MCMCpack"), 
       require, character.only=TRUE)
rds_obj <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB/combined.rds")
source("~/Dropbox/Work/MultiStudy_BatchEffect/TB/helper_crossmod.R")
command_args <- c("3", "4")
data_name <- "logfpkm"
set.seed(123)

tb_africa <- rds_obj[, rds_obj$SequencingBatch=="Africa"]
tb_africa_sub <- tb_africa[, tb_africa$Label %in% c("Non-progressor", "Progressor")]
tb_g6 <- rds_obj[, rds_obj$SequencingBatch=="G6"]

train_expr <- assays(tb_africa_sub)[[data_name]]
y_train <- revalue(as.factor(as.character(tb_africa_sub$Label)), 
                   c("Non-progressor"="0", "Progressor"="1"))
test_expr <- assays(tb_g6)[[data_name]]
y_test <- revalue(as.factor(as.character(tb_g6$Label)), c("Non-progressor"="0", "Progressor"="1"))
identical(rownames(train_expr), rownames(test_expr))
rm(tb_africa, tb_africa_sub, tb_g6)

tmp <- train_expr[rowVars(train_expr)!=0, ]

pca_ori_res <- prcomp(t(tmp), center=TRUE, scale.=TRUE)
pca_ori_plt_obj <- data.frame(PC1=pca_ori_res$x[, 1], PC2=pca_ori_res$x[, 2], 
                              Condition=revalue(y_train, c("0"="Non-progressor", "1"="Progressor")))
plt_ori <- ggplot(pca_ori_plt_obj, aes(x=PC1, y=PC2, shape=Condition)) +
  geom_point(size=2) +
  scale_shape_manual(values=c(16, 17)) +
  #theme_bw() +
  theme(legend.position = "none")+
  labs(x=sprintf("PC1: %s Variance", scales::percent(pca_ori_res$sdev[1] / sum(pca_ori_res$sdev))),
       y=sprintf("PC2: %s Variance", scales::percent(pca_ori_res$sdev[2] / sum(pca_ori_res$sdev))),
       title="Original Africa study")

var_trn <- rowVars(train_expr)
genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:1000]]
train_expr <- train_expr[genes_sel, ]
test_expr <- test_expr[genes_sel, ]

## Degree of batch effect (strength of signal)
N_batch <- 3
N_sample_size <- 20   # number of samples per batch (case + control)
max_batch_mean <- as.numeric(command_args[1]) 
max_batch_var <- as.numeric(command_args[2]) 
hyper_pars <- list(hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  
                   hyper_sd=sqrt(rep(0.01, N_batch)),
                   hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                     v=rep(0.01, N_batch))$alpha,  
                   hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                    v=rep(0.01, N_batch))$beta)  

## Subset training set in batches
batches_ind <- subsetBatch(condition=y_train, N_sample_size=N_sample_size, N_batch=N_batch)
y_sgbatch_train <- lapply(1:N_batch, function(k){y_train[batches_ind[[k]]]})
batch <- rep(0, ncol(train_expr)); for(i in 1:N_batch){batch[batches_ind[[i]]] <- i}
curr_train_expr <- train_expr[, do.call(c, batches_ind)]
curr_y_train <- y_train[do.call(c, batches_ind)]
batch <- batch[do.call(c, batches_ind)]  
batches_ind <- lapply(1:N_batch, function(i){which(batch==i)})

## Remove genes with only 0 values in any batch in current training set
g_keep <- lapply(1:N_batch, function(j){
  which(apply(curr_train_expr[, batch==j], 1, function(x){!all(x==0)}))
})
g_keep <- Reduce(intersect, g_keep)  
curr_train_expr <- curr_train_expr[g_keep, ]
curr_test_expr <- test_expr[g_keep, ]

## Simulate batch effect 
sim_batch_res <- simBatch(dat=curr_train_expr, condition=curr_y_train, 
                          batches_ind=batches_ind, batch=batch, hyper_pars=hyper_pars)
train_expr_batch <- sim_batch_res$new_dat

pca_res <- prcomp(t(train_expr_batch), center=TRUE, scale.=TRUE)
pca_plt_obj <- data.frame(PC1=pca_res$x[, 1], PC2=pca_res$x[, 2],
                          Condition=revalue(curr_y_train, c("0"="Non-progressor", "1"="Progressor")),
                          Batch=as.factor(batch))
plt_sim <- ggplot(pca_plt_obj, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=Batch, shape=Condition), size=2) +
  scale_shape_manual(values=c(16, 17)) +
  #theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical")+
  #legend.direction = "vertical") +
  #guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  labs(x=sprintf("PC1: %s Variance", scales::percent(pca_res$sdev[1] / sum(pca_res$sdev))),
       y=sprintf("PC2: %s Variance", scales::percent(pca_res$sdev[2] / sum(pca_res$sdev))),
       title="Example training set with 3 simulated batches")

png("Fig1.png", width=3000, height=1200, units="px", res=300)
ggarrange(plt_ori, plt_sim, ncol=2, widths=c(1,1.4))
dev.off()




########  Figure 2: ensemble single learner - RF  ########
rm(list=ls())
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "plyr", "ggpubr"), 
       require, character.only=TRUE)
results_dir <- "~/Documents/MSBE/TB_crossmod_reproduce/"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "nnet", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", 
                      "Neural Networks", "Support Vector Machines"); 
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 5)  
batch_var_vec <- c(1, 3, 5)
perf_measures <- c("mxe", "auc")  #, "acc", "f")
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")  #, "Accuracy", "F1 score")
names(perf_measures_plt) <- perf_measures
subset_colnames_single <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s") 
batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                           1, paste, collapse="_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                                   1, paste, collapse="_"))
curr_mod <- "rf"
plt_lst <- list()
for(curr_perf in perf_measures){
  curr_file_lst <-  sort(apply(expand.grid(method_names, paste0(curr_perf, "_batchN", N_sample_size, "_", 
                                                                curr_files_mv_suffix, ".csv")), 
                               1, paste, collapse="_"))
  print(curr_file_lst)
  sgmod_res_lst <- base_lst <- list()
  #bl = batch_levels[1]
  for(bl in batch_levels){
    curr_file_lst_bl <- grep(bl, curr_file_lst, value=TRUE)
    curr_res <- read.csv(paste0(results_dir, grep(curr_mod, curr_file_lst_bl, value=TRUE)), header=TRUE)
    tmp <- curr_res[, subset_colnames_single]
    tmp <- data.frame(melt(tmp), 
                      Type=rep(c("With batch effect,\nno adjustment", "Merging", rep("Ensemble", 3)), 
                               each=nrow(tmp)))
    sgmod_res_lst[[bl]] <- tmp
    base_lst[[bl]] <- curr_res[, "NoBatch"]
  }
  plt_df <- melt(sgmod_res_lst)
  
  plt_df$variable <- factor(plt_df$variable, levels=subset_colnames_single)
  plt_df$variable <- revalue(plt_df$variable, c("Batch" = "No adjustment", 
                                                "ComBat" = "Merge + ComBat",
                                                "n_Avg" = "Batch size weights",  
                                                "CS_Avg" = "Replicability weights", 
                                                "Reg_s" = "Stacking regression weights"))
  plt_df$L1 <- revalue(plt_df$L1, c("m0_v1"="Mean difference 0\nVariance fold change 1", 
                                    "m0_v3"="Mean difference 0\nVariance fold change 3",
                                    "m0_v5"="Mean difference 0\nVariance fold change 5",
                                    "m5_v1"="Mean difference 5\nVariance fold change 1", 
                                    "m5_v3"="Mean difference 5\nVariance fold change 3",
                                    "m5_v5"="Mean difference 5\nVariance fold change 5"))
  plt_df$Type <- factor(plt_df$Type, levels=c("With batch effect,\nno adjustment", "Merging", "Ensemble"))
  
  plt_lst[[curr_perf]] <- ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=median(do.call(c, base_lst)), 
                   linetype="Average AUC on data\nwithout simulated batch\neffect"), color="red") +
    scale_linetype_manual(name="limit", values=2,
                          guide=guide_legend(override.aes=list(color=c("red")))) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    facet_wrap(~L1, ncol=length(batch_var_vec)) +
    #ylim(0.55, max(plt_df$value))+
    theme(axis.text.x=element_text(angle=30, hjust=1),
          axis.title.x=element_blank(),
          legend.title=element_blank(),
          plot.margin = margin(0.1, 0.1, 0.1, 0.4, "cm")) +
    labs(y=perf_measures_plt[curr_perf]) 
}

png("Fig2.png", width=7, height=7, units="in", res=300)
#ggarrange(plt_lst[["auc"]], plt_lst[["mxe"]], ncol=2, common.legend=TRUE, legend="right")
plt_lst[["auc"]]
dev.off()  


## include cross-mod results (Fig2_v2)
rm(list=ls())
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "plyr", "ggpubr"), 
       require, character.only=TRUE)
results_dir <- "~/Documents/MSBE/TB_crossmod_reproduce/"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "nnet", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Neural Networks", "Support Vector Machines"); 
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 5)  
batch_var_vec <- c(1,3,5)
perf_measures <- c("mxe", "auc", "acc", "f")
perf_measures_plt <- c("Mean cross-entropy loss", "AUC", "Accuracy", "F1 score")
names(perf_measures_plt) <- perf_measures
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s") 

curr_perf <- "auc"
i <- 4
curr_mod <- "rf"
batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                           1, paste, collapse="_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <-  sort(apply(expand.grid(c("crossmod", method_names), paste0(curr_perf, "_batchN", N_sample_size, "_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))
print(curr_file_lst)
crossmod_res_lst <- base_lst <- list()
bl = batch_levels[1]
for(bl in batch_levels){
  curr_file_lst_bl <- curr_file_lst[grep(bl, curr_file_lst)]
  crossmod_res <- read.csv(paste0(results_dir, curr_file_lst_bl[1]), header=TRUE)[, subset_colnames_crossmod[-c(1:3)]]
  crossmod_res <- data.frame(melt(crossmod_res), Type="Ensemble (across learners)")
  
  tmp <- read.csv(paste0(results_dir, curr_file_lst_bl[i]), header=TRUE)[, c("Batch", "ComBat", "n_Avg","CS_Avg","Reg_s")]
  colnames(tmp) <- paste0(curr_mod, c(".Batch", ".ComBat", ".n_Avg",".CS_Avg",".Reg_s"))
  tmp <- data.frame(melt(tmp), Type=rep(c("With batch effect, no adjustment", "Merge + ComBat", rep("Ensemble (single learner)",3)), each=nrow(tmp)))
  crossmod_res_lst[[bl]] <- rbind(crossmod_res, tmp)
  
  base_lst[[bl]] <- read.csv(paste0(results_dir, curr_file_lst_bl[i]), header=TRUE)[, "NoBatch"]
}
plt_df <- melt(crossmod_res_lst)
plt_df$variable <- factor(plt_df$variable, levels=c(paste0(strsplit(curr_file_lst_bl[i],"_")[[1]][1], c(".Batch", ".ComBat", ".n_Avg",".CS_Avg",".Reg_s")), subset_colnames_crossmod[-c(1:3)]))
plt_df$Type <- factor(plt_df$Type, levels=c("With batch effect, no adjustment", "Merge + ComBat", 
                                            "Ensemble (single learner)", "Ensemble (across learners)"))
plt_df$variable <- revalue(plt_df$variable, c("rf.Batch" = "No adjustment (RF only)", 
                                              "rf.ComBat" = "Merge + ComBat (RF only)",
                                              "rf.n_Avg" = "Weighted average (RF only)",  
                                              "rf.CS_Avg" = "Rewards replicability (RF only)", 
                                              "rf.Reg_s" = "Regression stacking (RF only)",
                                              "n_Avg" = "Weighted average (across learners)",
                                              "CS_Avg" = "Rewards replicability (across learners)",
                                              "Reg_s" = "Regression stacking (across learners)"))
plt_df$L1 <- revalue(plt_df$L1, c("m0_v1"="Mean difference 0\nVariance fold change 1", 
                                  "m0_v3"="Mean difference 0\nVariance fold change 3",
                                  "m0_v5"="Mean difference 0\nVariance fold change 5",
                                  "m5_v1"="Mean difference 5\nVariance fold change 1", 
                                  "m5_v3"="Mean difference 5\nVariance fold change 3",
                                  "m5_v5"="Mean difference 5\nVariance fold change 5"))

png("Fig2_v2.png", width=8, height=7, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=median(do.call(c, base_lst)), 
                       linetype="Average AUC on data\nwithout simulated batch\neffect (RF only)"), color="red") +
        scale_linetype_manual(name="limit", values=2,
                              guide=guide_legend(override.aes=list(color=c("red")))) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "dark green")) +
        facet_wrap(~L1, ncol=length(batch_var_vec)) +
        ylim(0.5,max(plt_df$value))+
        #theme_bw() +
        theme(axis.text.x=element_text(angle=30, hjust=1, size=6),
              axis.title.x=element_blank(),
              legend.title=element_blank(),
              #legend.box = 'vertical',
              #legend.direction = 'vertical', 
              legend.position = "right",
              plot.margin = margin(0.1, 0.1, 0.1, 0.6, "cm")) +
        #guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
        labs(y=perf_measures_plt[curr_perf])) #, title=sprintf("Batch size: %s", N_sample_size))) 
dev.off()  




########  Figure 4: real data application  ########
rm(list=ls())
setwd("~/Documents/MSBE/TB_realdata_Testswap/")
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "dplyr", "plyr","ggpubr"), 
       require, character.only=TRUE)

study_names <- c("Brazil_1", "India", "Africa")
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- TRUE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")

#i=2; j=1
plt_crossmod <- res_lst <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_crossmod[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "ComBat"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "ComBat", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Weighted Average", "CS_Avg"="Rewards replicability",
                                                        "Reg_s"="Regression stacking"))
    
    sumstats_crossmod <- dplyr::filter(sumstats, Model=="crossmod")
    plt_crossmod[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      #facet_wrap(~Model, nrow=2, ncol=2) +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="auc", paste("Test set:", gsub("_1", "", curr_testset)), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}

png("Fig5.png", width=7, height=8, units="in", res=300)
ggarrange(plt_crossmod[["auc"]][["Africa"]], plt_crossmod[["mxe"]][["Africa"]], 
          plt_crossmod[["auc"]][["Brazil_1"]], plt_crossmod[["mxe"]][["Brazil_1"]], 
          plt_crossmod[["auc"]][["India"]], plt_crossmod[["mxe"]][["India"]], 
          ncol=2, nrow=3, common.legend=TRUE, legend="right")
dev.off()


## average across test studies
plt_pooled <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_pooled <- do.call(rbind, res_lst[[curr_perf]])
  sumstats_pooled <- res_pooled %>% 
    dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
    dplyr::group_by(Method, Model) %>% 
    dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
  sumstats_pooled$Type <- rep("Original Data", nrow(sumstats_pooled))
  sumstats_pooled$Type[sumstats_pooled$Method=="ComBat"] <- "ComBat"
  sumstats_pooled$Type[sumstats_pooled$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
  sumstats_pooled$Type <- factor(sumstats_pooled$Type, levels=c("Original Data", "ComBat", "Ensemble"))
  sumstats_pooled$Method <- factor(sumstats_pooled$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
  sumstats_pooled$Method <- plyr::revalue(sumstats_pooled$Method, 
                                          c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                            "n_Avg"="Weighted Average", "CS_Avg"="Rewards replicability",
                                            "Reg_s"="Regression stacking"))
  sumstats_pooled_crossmod <- dplyr::filter(sumstats_pooled, Model=="crossmod")
  plt_pooled[[curr_perf]] <- ggplot(sumstats_pooled_crossmod, aes(x=Method, y=Avg, color=Type)) + 
    geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
    geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
    geom_point() +
    labs(y=perf_measures_names[perf_measures[i]], 
         title=ifelse(curr_perf=="auc", "Averaging over test studies", "")) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1),
          legend.title=element_blank())
}

png("Fig5_b.png", width=7.5, height=3.5, units="in", res=300)
ggarrange(plt_pooled[["auc"]], plt_pooled[["mxe"]], 
          ncol=2, nrow=1, common.legend=TRUE, legend="right")
dev.off()




########  Table 3: calculate percentage in all bootstrap when ensemble out-performs combat  ########
rm(list=ls())
setwd("~/Documents/MSBE/TB_realdata_Testswap/")
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "dplyr", "plyr","ggpubr"), 
       require, character.only=TRUE)

study_names <- c("Africa", "Brazil_1", "India")
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- TRUE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")

#j=1; b=1
curr_perf <- "auc"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- filter(res, Iteration==b, Model=="crossmod")
    curr_count <- as.numeric(c(
      Avg=curr_iter[curr_iter$Method=="Avg", "value"] > curr_iter[curr_iter$Method=="ComBat", "value"],
      n_Avg=curr_iter[curr_iter$Method=="n_Avg", "value"] > curr_iter[curr_iter$Method=="ComBat", "value"],
      CS_Avg=curr_iter[curr_iter$Method=="CS_Avg", "value"] > curr_iter[curr_iter$Method=="ComBat", "value"],
      Reg_a=curr_iter[curr_iter$Method=="Reg_a", "value"] > curr_iter[curr_iter$Method=="ComBat", "value"],
      Reg_s=curr_iter[curr_iter$Method=="Reg_s", "value"] > curr_iter[curr_iter$Method=="ComBat", "value"])
    )
    total_count <- total_count + curr_count
  }
  names(total_count) <- c("Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}

do.call(rbind, freq_lst)[, c("n_Avg", "CS_Avg", "Reg_s")]




########  Table 2: collect simulated batch moments  ########
rm(list=ls())
setwd("~/Dropbox/Work/MultiStudy_BatchEffect/TB/")
rds_obj <- readRDS("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB/combined.rds")
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "limma", 
             "nnls",  "glmnet", "rpart", "genefilter", "nnet", "e1071", "RcppArmadillo", 
             "foreach", "parallel", "doParallel", "MLmetrics", "ranger", "scales"), 
           require, character.only=TRUE))
source("helper_crossmod.R")
command_args <- c("5", "8")
data_name <- "logfpkm"   #"logtpm"
tb_africa <- rds_obj[, rds_obj$SequencingBatch=="Africa"]
tb_africa_sub <- tb_africa[, tb_africa$Label %in% c("Non-progressor", "Progressor")]
tb_g6 <- rds_obj[, rds_obj$SequencingBatch=="G6"]

train_expr <- assays(tb_africa_sub)[[data_name]]
y_train <- revalue(as.factor(as.character(tb_africa_sub$Label)), 
                   c("Non-progressor"="0", "Progressor"="1"))
test_expr <- assays(tb_g6)[[data_name]]
y_test <- revalue(as.factor(as.character(tb_g6$Label)), c("Non-progressor"="0", "Progressor"="1"))
identical(rownames(train_expr), rownames(test_expr))
rm(tb_africa, tb_africa_sub, tb_g6)

# select 300 genes with largest variance in training set (Africa)
var_trn <- rowVars(train_expr)
genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:300]]
train_expr <- train_expr[genes_sel, ]
test_expr <- test_expr[genes_sel, ]

N_batch <- 3
N_sample_size <- 20   # number of samples per batch (case + control)
max_batch_mean <- as.numeric(command_args[1]) 
max_batch_var <- as.numeric(command_args[2]) 
hyper_pars <- list(hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  
                   hyper_sd=sqrt(rep(0.01, N_batch)),
                   hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                     v=rep(0.01, N_batch))$alpha,  
                   hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                                    v=rep(0.01, N_batch))$beta)  

## Subset training set in batches
batches_ind <- subsetBatch(condition=y_train, N_sample_size=N_sample_size, N_batch=N_batch)
y_sgbatch_train <- lapply(1:N_batch, function(k){y_train[batches_ind[[k]]]})
batch <- rep(0, ncol(train_expr)); for(i in 1:N_batch){batch[batches_ind[[i]]] <- i}
curr_train_expr <- train_expr[, do.call(c, batches_ind)]
curr_y_train <- y_train[do.call(c, batches_ind)]
batch <- batch[do.call(c, batches_ind)]  
batches_ind <- lapply(1:N_batch, function(i){which(batch==i)})

## Remove genes with only 0 values in any batch in current training set
g_keep <- lapply(1:N_batch, function(j){
  which(apply(curr_train_expr[, batch==j], 1, function(x){!all(x==0)}))
})
g_keep <- Reduce(intersect, g_keep)  
curr_train_expr <- curr_train_expr[g_keep, ]
curr_test_expr <- test_expr[g_keep, ]

## Simulate batch effect 
sim_batch_res <- simBatch(dat=curr_train_expr, condition=curr_y_train, 
                          batches_ind=batches_ind, batch=batch, hyper_pars=hyper_pars)
train_expr_batch <- sim_batch_res$new_dat

res <- cbind(
  round(sapply(1:N_batch, function(i){mean(apply(train_expr_batch[, batches_ind[[i]]], 1, mean))}), 2),
  round(sapply(1:N_batch, function(i){mean(apply(train_expr_batch[, batches_ind[[i]]], 1, var))}), 2)
)
res




########  Figure 3: sample size comparison   ########
rm(list=ls())
results_dir <- "~/Documents/MSBE/TB_crossmod_reproduce"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines"); 
names(method_names_plt) <- method_names
perf_measures <- c("mxe", "auc")
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")
names(perf_measures_plt) <- perf_measures
batch_mean_vec <- 3
batch_var_vec <- 1:5
batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_perf <- "auc"
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s") 

curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <- sort(c(sort(apply(expand.grid(c("crossmod", method_names), paste0(curr_perf, "_batchN20_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_")),
                        sort(apply(expand.grid(c("crossmod", method_names), paste0(curr_perf, "_batchN40_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))))
print(curr_file_lst)

#bl=batch_levels[1]
crossmod_res_lst <- list() 
for(bl in batch_levels){
  curr_file_lst_bl <- grep(bl, curr_file_lst, value=TRUE)
  #curr_fname = curr_file_lst_bl[1]
  crossmod_res <- lapply(curr_file_lst_bl, function(curr_fname){
    res <- read.csv(file.path(results_dir, curr_fname), header=TRUE)
    curr_f_id <- strsplit(curr_fname, "_")[[1]]
    if(curr_f_id[1]=="crossmod"){
      res <- res[, subset_colnames_crossmod[-c(1:3)]]
      res <- data.frame(melt(res), Type="Ensemble", N_sample=curr_f_id[3])
    }else if(curr_f_id[1] %in% c("lasso", "rf")){
      res <- res[, c("NoBatch", "Batch", "ComBat")]
      colnames(res) <- paste0(curr_f_id[1], c(".NoBatch", ".Batch", ".ComBat"))
      res <- data.frame(melt(res), Type=rep(c("Baseline", "Batch", "ComBat"), each=nrow(res)), N_sample=curr_f_id[3])
    }else{
      res <- NULL
    }
    return(res)
  })
  crossmod_res <- crossmod_res[sapply(crossmod_res, function(res){!is.null(res)})]
  crossmod_res_lst[[bl]] <- do.call(rbind, crossmod_res)
}
names(crossmod_res_lst) <- batch_levels
plt_df <- melt(crossmod_res_lst)
plt_df$variable <- factor(plt_df$variable, levels=c("lasso.NoBatch", "rf.NoBatch", "lasso.Batch", "rf.Batch", 
                                                    "lasso.ComBat", "rf.ComBat", subset_colnames_crossmod[-c(1:3)]))
plt_df$Type <- factor(plt_df$Type, levels=c("Baseline", "Batch", "ComBat", "Ensemble"))

plt_df$variable <- revalue(plt_df$variable, c("lasso.NoBatch"="No batch effect, LASSO",
                                              "rf.NoBatch"="No batch effect, RF",
                                              "lasso.Batch"="With batch, no adjustment, LASSO",
                                              "rf.Batch"="With batch, no adjustment, RF",
                                              "lasso.ComBat"="Merge + ComBat, LASSO",
                                              "rf.ComBat"="Merge + ComBat, RF",
                                              "n_Avg"="Weighted average (across learners)",
                                              "CS_Avg"="Replicability weights (across learners)",
                                              "Reg_s"="Regression stacking (across learners)"))
plt_df$Type <- revalue(plt_df$Type, c("Baseline"="No batch effect",
                                      "Batch"="With batch, no adjustment",
                                      "ComBat"="Merge + ComBat",
                                      "Ensemble" = "Ensemble (across learners)"))
plt_df$N_sample <- revalue(plt_df$N_sample, c("batchN20"="Batch size 20", "batchN40"="Batch size 40"))
plt_df$L1 <- revalue(plt_df$L1, c("m3_v1"="Mean difference 3\nVariance fold change 1", 
                                  "m3_v2"="Mean difference 3\nVariance fold change 2",
                                  "m3_v3"="Mean difference 3\nVariance fold change 3",
                                  "m3_v4"="Mean difference 3\nVariance fold change 4",
                                  "m3_v5"="Mean difference 3\nVariance fold change 5"))
                            
png("Fig3.png", width=12, height=7, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        scale_fill_manual(values=c("white","#999999", "#E69F00", "#56B4E9")) +
        facet_grid(N_sample~L1) +
        theme(axis.text.x=element_text(angle=40, hjust=1),
              axis.title.x=element_blank(), 
              legend.title=element_blank(),
              legend.position="bottom",
              plot.margin = margin(0.1, 0.1, 0.1, 1, "cm")) +
        labs(y=perf_measures_plt[curr_perf])) 
dev.off()  




########  Supplementary figure 2: full simulation results   ########
rm(list=ls())
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "plyr", "ggpubr"), 
       require, character.only=TRUE)
results_dir <- "~/Documents/MSBE/TB_crossmod_reproduce/"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines"); 
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 3, 5)  
batch_var_vec <- 1:5 #c(1,3,5)
perf_measures <- c("mxe", "auc")   #, "acc", "f")
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")   #, "Accuracy", "F1 score")
names(perf_measures_plt) <- perf_measures
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s") 
curr_perf <- "auc"

#i = 2; curr_mod = 'lasso'; ylim_l = 0.4
#i = 3; curr_mod = 'rf'; ylim_l = 0.5
i = 4; curr_mod = 'svm'; ylim_l = 0.3
plt_label_vec <- c("svm.Batch" = "No adjustment",  # (LASSO only)", 
                   "svm.ComBat" = "Merge + ComBat",  # (LASSO only)",
                   "svm.Avg" = "Simple average",  # (LASSO only)",
                   "svm.n_Avg" = "Weighted average",  # (LASSO only)",  
                   "svm.CS_Avg" = "Rewards replicability",  # (LASSO only)", 
                   "svm.Reg_a" = "Regression aggregate",  # (LASSO only)",
                   "svm.Reg_s" = "Regression stacking",  # (LASSO only)",
                   "Avg" = "Simple average (across learners)",
                   "n_Avg" = "Weighted average (across learners)",  
                   "CS_Avg" = "Rewards replicability (across learners)", 
                   "Reg_a" = "Regression aggregate (across learners)",
                   "Reg_s" = "Regression stacking (across learners)")

batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                           1, paste, collapse="_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <-  sort(apply(expand.grid(c("crossmod", method_names), paste0(curr_perf, "_batchN", N_sample_size, "_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))
print(curr_file_lst)
crossmod_res_lst <- base_lst <- list()
bl = batch_levels[1]
for(bl in batch_levels){
  curr_file_lst_bl <- curr_file_lst[grep(bl, curr_file_lst)]
  crossmod_res <- read.csv(paste0(results_dir, curr_file_lst_bl[1]), header=TRUE)[, subset_colnames_crossmod[-c(1:3)]]
  crossmod_res <- data.frame(melt(crossmod_res), Type="Ensemble (across learners)")
  
  tmp <- read.csv(paste0(results_dir, curr_file_lst_bl[i]), header=TRUE)[, subset_colnames_crossmod[-1]]
  colnames(tmp) <- paste0(curr_mod, '.', subset_colnames_crossmod[-1])
  tmp <- data.frame(melt(tmp), Type=rep(c("With batch effect, no adjustment", "Merge + ComBat", rep("Ensemble (single learner)",5)), each=nrow(tmp)))
  crossmod_res_lst[[bl]] <- rbind(crossmod_res, tmp)
  
  base_lst[[bl]] <- read.csv(paste0(results_dir, curr_file_lst_bl[i]), header=TRUE)[, "NoBatch"]
}
plt_df <- melt(crossmod_res_lst)
plt_df$variable <- factor(plt_df$variable, levels=c(paste0(strsplit(curr_file_lst_bl[i],"_")[[1]][1], '.', subset_colnames_crossmod[-1]), subset_colnames_crossmod[-c(1:3)]))
plt_df$Type <- factor(plt_df$Type, levels=c("With batch effect, no adjustment", "Merge + ComBat", 
                                            "Ensemble (single learner)", "Ensemble (across learners)"))
plt_df$variable <- revalue(plt_df$variable, plt_label_vec)
plt_df$L1 <- revalue(plt_df$L1, c("m0_v1"="Mean difference 0\nVariance fold change 1", 
                                  "m0_v2"="Mean difference 0\nVariance fold change 2",
                                  "m0_v3"="Mean difference 0\nVariance fold change 3",
                                  "m0_v4"="Mean difference 0\nVariance fold change 4",
                                  "m0_v5"="Mean difference 0\nVariance fold change 5",
                                  #"m0_v6"="Mean difference 0\nVariance fold change 6",
                                  #"m0_v7"="Mean difference 0\nVariance fold change 7",
                                  #"m0_v8"="Mean difference 0\nVariance fold change 8",
                                  "m3_v1"="Mean difference 3\nVariance fold change 1", 
                                  "m3_v2"="Mean difference 3\nVariance fold change 2",
                                  "m3_v3"="Mean difference 3\nVariance fold change 3",
                                  "m3_v4"="Mean difference 3\nVariance fold change 4",
                                  "m3_v5"="Mean difference 3\nVariance fold change 5",
                                  #"m3_v6"="Mean difference 3\nVariance fold change 6",
                                  #"m3_v7"="Mean difference 3\nVariance fold change 7",
                                  #"m3_v8"="Mean difference 3\nVariance fold change 8",
                                  "m5_v1"="Mean difference 5\nVariance fold change 1", 
                                  "m5_v2"="Mean difference 5\nVariance fold change 2",
                                  "m5_v3"="Mean difference 5\nVariance fold change 3",
                                  "m5_v4"="Mean difference 5\nVariance fold change 4",
                                  "m5_v5"="Mean difference 5\nVariance fold change 5"))#,
                                  #"m5_v6"="Mean difference 5\nVariance fold change 6",
                                  #"m5_v7"="Mean difference 5\nVariance fold change 7",
                                  #"m5_v8"="Mean difference 5\nVariance fold change 8"))

png(sprintf("FigS2_v15_%s.png", curr_mod), width=9, height=9, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=median(do.call(c, base_lst)), 
                       linetype="Average AUC on data\nwithout simulated batch\neffect (RF only)"), color="red") +
        scale_linetype_manual(name="limit", values=2,
                              guide=guide_legend(override.aes=list(color=c("red")))) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "dark green")) +
        facet_wrap(~L1, ncol=length(batch_var_vec)) +
        ylim(ylim_l, max(plt_df$value))+
        theme(axis.text.x=element_text(angle=30, hjust=1, size=6),
              axis.title.x=element_blank(),
              legend.title=element_blank(),
              legend.position = "bottom",
              plot.margin = margin(0.1, 0.1, 0.1, 0.3, "cm")) +
        labs(y=perf_measures_plt[curr_perf]))
dev.off()  
