rm(list=ls())
setwd("./")
if(!dir.exists("./figures")){dir.create("./figures")}
all(sapply(c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", "reshape2", 
             "ggpubr", "gridExtra","limma", "nnls",  "glmnet", "rpart", "genefilter", "nnet", 
             "e1071", "RcppArmadillo", "foreach", "parallel", "doParallel", "MLmetrics", 
             "ranger", "scales"), 
           require, character.only=TRUE))


########  Supplementary figure S2: PCA of training data  ########
source("./code/helper.R")
load("./data/combined_sub.RData")
command_args <- c("3", "4")
set.seed(123)

tmp <- train_expr[rowVars(train_expr)!=0, ]
pca_ori_res <- prcomp(t(tmp), center=TRUE, scale.=TRUE)
pca_ori_plt_obj <- data.frame(PC1=pca_ori_res$x[, 1], PC2=pca_ori_res$x[, 2], 
                              Condition=revalue(y_train, c("0"="Non-progressor", "1"="Progressor")))
plt_ori <- ggplot(pca_ori_plt_obj, aes(x=PC1, y=PC2, shape=Condition)) +
  geom_point(size=2) +
  scale_shape_manual(values=c(16, 17)) +
  theme(legend.position = "none")+
  labs(x=sprintf("PC1: %s Variance", scales::percent(pca_ori_res$sdev[1] / sum(pca_ori_res$sdev))),
       y=sprintf("PC2: %s Variance", scales::percent(pca_ori_res$sdev[2] / sum(pca_ori_res$sdev))),
       title="Original study A")

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
  theme(legend.position = "right",
        legend.box = "vertical")+
  labs(x=sprintf("PC1: %s Variance", scales::percent(pca_res$sdev[1] / sum(pca_res$sdev))),
       y=sprintf("PC2: %s Variance", scales::percent(pca_res$sdev[2] / sum(pca_res$sdev))),
       title="Example training set with 3 simulated batches")

png("./figures/FigS2.png", width=3000, height=1200, units="px", res=300)
ggarrange(plt_ori, plt_sim, ncol=2, widths=c(1,1.4))
dev.off()




########  Supplementary table S2: collect simulated batch moments  ########
rm(list=ls())
load("./data/combined_sub.RData")
source("./code/helper.R")
command_args <- c("3", "5") # mean and variance fold change

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
rownames(res) <- paste("batch", 1:3)
sprintf("Curr mean FC: %s, var FC: %s", command_args[1], command_args[2])
t(data.frame(mean=res[,1], var=res[,2]))




########  Supplementary figure S3: LASSO full simulation results   ########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_sim/"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines"); 
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 3, 5)  
batch_var_vec <- 1:5 
perf_measures <- c("mxe", "auc")   
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")   
names(perf_measures_plt) <- perf_measures
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s") 
curr_perf <- "auc"

i = 2; curr_mod = 'lasso'
plt_label_vec <- c("lasso.Batch" = "No adjustment",   
                   "lasso.ComBat" = "Merge + ComBat",  
                   "lasso.Avg" = "Simple average",  
                   "lasso.n_Avg" = "Batch size",   
                   "lasso.CS_Avg" = "Cross-study",  
                   "lasso.Reg_a" = "Regression aggregate",  
                   "lasso.Reg_s" = "Regression stacking", 
                   "Avg" = "Simple average (across learners)",
                   "n_Avg" = "Batch size weights (across learners)",  
                   "CS_Avg" = "Cross-study weights (across learners)", 
                   "Reg_a" = "Regression aggregate (across learners)",
                   "Reg_s" = "Regression stacking (across learners)")

batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                           1, paste, collapse="_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <-  sort(apply(expand.grid(c("crossmod", method_names), 
                                         paste0(curr_perf, "_batchN", N_sample_size, "_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))
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
                                  "m3_v1"="Mean difference 3\nVariance fold change 1", 
                                  "m3_v2"="Mean difference 3\nVariance fold change 2",
                                  "m3_v3"="Mean difference 3\nVariance fold change 3",
                                  "m3_v4"="Mean difference 3\nVariance fold change 4",
                                  "m3_v5"="Mean difference 3\nVariance fold change 5",
                                  "m5_v1"="Mean difference 5\nVariance fold change 1", 
                                  "m5_v2"="Mean difference 5\nVariance fold change 2",
                                  "m5_v3"="Mean difference 5\nVariance fold change 3",
                                  "m5_v4"="Mean difference 5\nVariance fold change 4",
                                  "m5_v5"="Mean difference 5\nVariance fold change 5"))#,

png("./figures/FigS3.png", width=9, height=9, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=median(do.call(c, base_lst)), 
                       linetype="Average AUC on data\nwithout simulated batch\neffect (RF only)"), color="red") +
        scale_linetype_manual(name="limit", values=2,
                              guide=guide_legend(override.aes=list(color=c("red")))) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "dark green")) +
        facet_wrap(~L1, ncol=length(batch_var_vec)) +
        ylim(0.45, 0.75)+
        theme(axis.text.x=element_text(angle=30, hjust=1, size=6),
              axis.title.x=element_blank(),
              legend.title=element_blank(),
              legend.position = "bottom",
              plot.margin = margin(0.1, 0.1, 0.1, 0.3, "cm")) +
        labs(y=perf_measures_plt[curr_perf]))
dev.off()  




########  Supplementary figure S4: RF full simulation results   ########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_sim/"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines"); 
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 3, 5)  
batch_var_vec <- 1:5 
perf_measures <- c("mxe", "auc")   
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")   
names(perf_measures_plt) <- perf_measures
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s") 
curr_perf <- "auc"

i = 3; curr_mod = 'rf'
plt_label_vec <- c("rf.Batch" = "No adjustment",  
                   "rf.ComBat" = "Merge + ComBat",  
                   "rf.Avg" = "Simple average",  
                   "rf.n_Avg" = "Batch size weights",    
                   "rf.CS_Avg" = "Cross-study weights", 
                   "rf.Reg_a" = "Regression aggregate",  
                   "rf.Reg_s" = "Regression stacking",  
                   "Avg" = "Simple average (across learners)",
                   "n_Avg" = "Batch size weights (across learners)",  
                   "CS_Avg" = "Cross-study weights (across learners)", 
                   "Reg_a" = "Regression aggregate (across learners)",
                   "Reg_s" = "Regression stacking (across learners)")

batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                           1, paste, collapse="_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <-  sort(apply(expand.grid(c("crossmod", method_names), 
                                         paste0(curr_perf, "_batchN", N_sample_size, "_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))
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
                                  "m3_v1"="Mean difference 3\nVariance fold change 1", 
                                  "m3_v2"="Mean difference 3\nVariance fold change 2",
                                  "m3_v3"="Mean difference 3\nVariance fold change 3",
                                  "m3_v4"="Mean difference 3\nVariance fold change 4",
                                  "m3_v5"="Mean difference 3\nVariance fold change 5",
                                  "m5_v1"="Mean difference 5\nVariance fold change 1", 
                                  "m5_v2"="Mean difference 5\nVariance fold change 2",
                                  "m5_v3"="Mean difference 5\nVariance fold change 3",
                                  "m5_v4"="Mean difference 5\nVariance fold change 4",
                                  "m5_v5"="Mean difference 5\nVariance fold change 5"))#,

png("./figures/FigS4.png", width=9, height=9, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=median(do.call(c, base_lst)), 
                       linetype="Average AUC on data\nwithout simulated batch\neffect (RF only)"), color="red") +
        scale_linetype_manual(name="limit", values=2,
                              guide=guide_legend(override.aes=list(color=c("red")))) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "dark green")) +
        facet_wrap(~L1, ncol=length(batch_var_vec)) +
        ylim(0.45, 0.75)+
        theme(axis.text.x=element_text(angle=30, hjust=1, size=6),
              axis.title.x=element_blank(),
              legend.title=element_blank(),
              legend.position = "bottom",
              plot.margin = margin(0.1, 0.1, 0.1, 0.3, "cm")) +
        labs(y=perf_measures_plt[curr_perf]))
dev.off()  




########  Supplementary figure S5: SVM full simulation results   ########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_sim/"
file_lst <- grep(".csv", dir(results_dir), fixed=TRUE, value=TRUE)  # all results files
method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines"); 
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 3, 5)  
batch_var_vec <- 1:5 
perf_measures <- c("mxe", "auc")   
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")   
names(perf_measures_plt) <- perf_measures
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s") 
curr_perf <- "auc"

i = 4; curr_mod = 'svm'
plt_label_vec <- c("svm.Batch" = "No adjustment",   
                   "svm.ComBat" = "Merge + ComBat",  
                   "svm.Avg" = "Simple average",  
                   "svm.n_Avg" = "Batch size weights",   
                   "svm.CS_Avg" = "Cross-study weights", 
                   "svm.Reg_a" = "Regression aggregate", 
                   "svm.Reg_s" = "Regression stacking",  
                   "Avg" = "Simple average (across learners)",
                   "n_Avg" = "Batch size weights (across learners)",  
                   "CS_Avg" = "Cross-study weights (across learners)", 
                   "Reg_a" = "Regression aggregate (across learners)",
                   "Reg_s" = "Regression stacking (across learners)")

batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 
                           1, paste, collapse="_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <-  sort(apply(expand.grid(c("crossmod", method_names), 
                                         paste0(curr_perf, "_batchN", N_sample_size, "_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))
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
                                  "m3_v1"="Mean difference 3\nVariance fold change 1", 
                                  "m3_v2"="Mean difference 3\nVariance fold change 2",
                                  "m3_v3"="Mean difference 3\nVariance fold change 3",
                                  "m3_v4"="Mean difference 3\nVariance fold change 4",
                                  "m3_v5"="Mean difference 3\nVariance fold change 5",
                                  "m5_v1"="Mean difference 5\nVariance fold change 1", 
                                  "m5_v2"="Mean difference 5\nVariance fold change 2",
                                  "m5_v3"="Mean difference 5\nVariance fold change 3",
                                  "m5_v4"="Mean difference 5\nVariance fold change 4",
                                  "m5_v5"="Mean difference 5\nVariance fold change 5"))#,

png("./figures/FigS5.png", width=9, height=9, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        geom_hline(aes(yintercept=median(do.call(c, base_lst)), 
                       linetype="Average AUC on data\nwithout simulated batch\neffect (RF only)"), color="red") +
        scale_linetype_manual(name="limit", values=2,
                              guide=guide_legend(override.aes=list(color=c("red")))) +
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "dark green")) +
        facet_wrap(~L1, ncol=length(batch_var_vec)) +
        ylim(0.45, 0.75)+
        theme(axis.text.x=element_text(angle=30, hjust=1, size=6),
              axis.title.x=element_blank(),
              legend.title=element_blank(),
              legend.position = "bottom",
              plot.margin = margin(0.1, 0.1, 0.1, 0.3, "cm")) +
        labs(y=perf_measures_plt[curr_perf]))
dev.off()  




########  Supplementary figure S6: number of batches comparison   ########
rm(list=ls())
source("./code/helper.R")
results_dir_3b <- "./results_sim/"
results_dir_5b <- "~/Documents/MSBE/TB_crossmod_reproduce_5batches/"
file_lst <- grep(".csv", dir(results_dir_3b), fixed=TRUE, value=TRUE)  # all results files

method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines"); 
names(method_names_plt) <- method_names
perf_measures <- c("mxe", "auc")
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")
names(perf_measures_plt) <- perf_measures
N_sample_size <- 20
batch_mean_vec <- 3
batch_var_vec <- c(2,4,6,8)
batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_perf <- "auc"
subset_colnames_crossmod <- c("NoBatch", "Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s") 

curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
curr_file_lst <- sort(sort(apply(expand.grid(c("crossmod", method_names), paste0(curr_perf, "_batchN20_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_")))
print(curr_file_lst)

#bl=batch_levels[1]
crossmod_res_lst <- list() 
for(bl in batch_levels){
  curr_file_lst_bl <- grep(bl, curr_file_lst, value=TRUE)
  #curr_fname = curr_file_lst_bl[1]
  crossmod_res <- lapply(curr_file_lst_bl, function(curr_fname){
    res_3b <- read.csv(file.path(results_dir_3b, curr_fname), header=TRUE)
    res_5b <- read.csv(file.path(results_dir_5b, curr_fname), header=TRUE)
    
    curr_f_id <- strsplit(curr_fname, "_")[[1]]
    if(curr_f_id[1]=="crossmod"){
      res <- rbind(data.frame(res_3b[, subset_colnames_crossmod[-c(1:3)]], Nbatch="3 simulated batches"),
                   data.frame(res_5b[, subset_colnames_crossmod[-c(1:3)]], Nbatch="5 simulated batches"))
      res <- data.frame(melt(res), Type="Ensemble")
    }else if(curr_f_id[1] %in% c("lasso", "rf")){
      res <- rbind(data.frame(res_3b[, c("NoBatch", "Batch", "ComBat")], Nbatch="3 simulated batches"),
                   data.frame(res_5b[, c("NoBatch", "Batch", "ComBat")], Nbatch="5 simulated batches"))
      colnames(res)[1:3] <- paste0(curr_f_id[1], c(".NoBatch", ".Batch", ".ComBat"))
      res <- data.frame(melt(res), Type=rep(c("Baseline", "Batch", "ComBat"), each=nrow(res)))
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
                                              "n_Avg"="Batch size weights (across learners)",
                                              "CS_Avg"="Cross-study weights (across learners)",
                                              "Reg_s"="Regression stacking (across learners)"))
plt_df$Type <- revalue(plt_df$Type, c("Baseline"="No batch effect",
                                      "Batch"="With batch, no adjustment",
                                      "ComBat"="Merge + ComBat",
                                      "Ensemble" = "Ensemble (across learners)"))
plt_df$L1 <- revalue(plt_df$L1, c("m3_v2"="Mean difference 3\nVariance fold change 2",
  "m3_v4"="Mean difference 3\nVariance fold change 4",
  "m3_v6"="Mean difference 3\nVariance fold change 6",
  "m3_v8"="Mean difference 3\nVariance fold change 8"))

png("./figures/FigS6.png", width=10, height=7, units="in", res=300)
print(ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
        geom_boxplot() +
        scale_fill_manual(values=c("white","#999999", "#E69F00", "#56B4E9")) +
        facet_grid(Nbatch~L1) +
        theme(axis.text.x=element_text(angle=40, hjust=1),
              axis.title.x=element_blank(), 
              legend.title=element_blank(),
              legend.position="bottom",
              plot.margin = margin(0.1, 0.1, 0.1, 1, "cm")) +
        labs(y=perf_measures_plt[curr_perf])) 
dev.off()  




########  Supplementary figure S7: batch size comparison   ########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_sim/"
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
curr_file_lst <- sort(c(sort(apply(expand.grid(c("crossmod", method_names), 
                                               paste0(curr_perf, "_batchN20_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_")),
                        sort(apply(expand.grid(c("crossmod", method_names), 
                                               paste0(curr_perf, "_batchN40_", curr_files_mv_suffix, ".csv")), 1, paste, collapse="_"))))
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
                                              "n_Avg"="Batch size weights (across learners)",
                                              "CS_Avg"="Cross-study weights (across learners)",
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

png("./figures/FigS7.png", width=11, height=7, units="in", res=300)
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




########  Supplementary figure S8: 6-study mxe and weights (rf only)   ########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_real_6studies/"
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
learner_types <- c("lasso", "rf", "svm", "crossmod")

plt <- function(weights_df){
  mlt_df <- melt(weights_df)
  colnames(mlt_df)[1:2] <- c("Type", "Study")
  mlt_df$Study <- factor(mlt_df$Study, levels=c("A", "C", "D", "E", "F", "G"))
  colors <- c("A"="#999999", "C"="#E69F00", "D"="#56B4E9", 
              "E"="#009E73", "F"="#F0E442", "G"="#0072B2")
  p <- ggplot(mlt_df, aes(x=Type, y=value, fill=Study)) +
    geom_bar(stat="identity") + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1)) +
    scale_fill_manual(values = colors) +
    labs(y="Weights")
  return(p)
}

#STUDY <- "GSE37250_SA"
fig_lst <- list()
for(STUDY in study_names){
  other_study_names <- setdiff(study_names, STUDY)
  
  #i=2
  plt_lst <- res_lst <- list()
  for(i in 1:length(perf_measures)){
    curr_perf <- perf_measures[i]
    
    file_prefix <- sprintf("test%s", STUDY)
    res_lst[[curr_perf]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "ComBat"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "ComBat", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    #STUDY_NAME <- ifelse(STUDY=="GSE39941_M", "GSE39941", STUDY)
    for(MODEL in learner_types){
      sumstats_curr <- dplyr::filter(sumstats, Model==MODEL)
      plt_lst[[curr_perf]][[MODEL]] <- ggplot(sumstats_curr, aes(x=Method, y=Avg, color=Type)) +
        geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
        geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
        geom_point() +
        labs(y=perf_measures_names[perf_measures[i]],
             title=paste("Test set:", study_label[STUDY])) +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(angle=30, hjust=1),
              legend.title=element_blank())
    }
  }
  
  load(sprintf("./results_real_6studies/test%s_weights.RData", STUDY))
  
  lasso_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                                `Cross-study weights`=cs_weights_seq$lasso,  
                                `Stacked regression weights`=reg_s_beta$lasso))
  rf_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                             `Cross-study weights`=cs_weights_seq$rf, 
                             `Stacked regression weights`=reg_s_beta$rf))
  svm_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                              `Cross-study weights`=cs_weights_seq$svm,  
                              `Stacked regression weights`=reg_s_beta$svm))
  colnames(lasso_weights) <- colnames(rf_weights) <- colnames(svm_weights) <- study_label[other_study_names]
  
  lasso_w_plt <- plt(lasso_weights)
  rf_w_plt <- plt(rf_weights)
  svm_w_plt <- plt(svm_weights)
  
  curr_fig <- ggarrange(plt_lst[["mxe"]][["rf"]], rf_w_plt, ncol=2, nrow=1)
  annotate_figure(curr_fig, top = text_grob(paste("Test set:", gsub("_1", "", STUDY))))
  fig_lst[[STUDY]] <- curr_fig
}

png("./figures/FigS8.png", width=7, height=15, units="in", res=300)
ggarrange(fig_lst[["Africa"]], fig_lst[["GSE39941_M"]],
          fig_lst[["India"]], fig_lst[["US"]],
          fig_lst[["GSE37250_SA"]], fig_lst[["GSE37250_M"]],
          nrow=6, ncol=1)
dev.off()








########  Supplementary figure S9: AUC ##########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_real_4studies/"
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"

## calculate frequency table
#j=1; b=1
curr_perf <- "auc"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.max(curr_iter$value)] <- total_count[which.max(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label


plt_4s <- res_lst <- list()
i=1;j=1
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_4s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size weights"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study weights"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Regression stacking"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    
    plt_4s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=0.6), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf==curr_perf, paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}

#### 6 studies
results_dir <- "./results_real_6studies/"
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"

curr_perf <- "auc"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.max(curr_iter$value)] <- total_count[which.max(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label

plt_6s <- res_lst <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_6s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size weights"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study weights"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Regression stacking"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    
    plt_6s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=1), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf==curr_perf, paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}

lgd <- get_legend(plt_4s[["auc"]][["GSE37250_SA"]])
plt_4studies <- ggarrange(plt_4s[["auc"]][["India"]], plt_4s[["auc"]][["US"]], 
                          plt_4s[["auc"]][["GSE37250_SA"]], plt_4s[["auc"]][["GSE37250_M"]], 
                          ncol=4, nrow=1, common.legend=TRUE, legend="none")
plt_4studies <- annotate_figure(plt_4studies, 
                                top=text_grob("4 studies", face="bold", size=18, hjust=1, x=0.15))
plt_6studies <- ggarrange(plt_6s[["auc"]][["India"]], plt_6s[["auc"]][["US"]], 
                          plt_6s[["auc"]][["GSE37250_SA"]], plt_6s[["auc"]][["GSE37250_M"]], 
                          plt_6s[["auc"]][["Africa"]], plt_6s[["auc"]][["GSE39941_M"]],
                          ncol=6, nrow=1, common.legend=TRUE, legend="none")
plt_6studies <- annotate_figure(plt_6studies, 
                                top=text_grob("6 studies", face="bold", size=18, hjust=1, x=0.1))
p <- arrangeGrob(plt_4studies, plt_6studies, lgd,                              
                 ncol=3, nrow=3, layout_matrix=rbind(c(2,2,2), c(NA,NA,NA), c(1,1,3)),
                 heights=c(1,0.1,1))

png("./figures/FigS9.png", width=13, height=6, units="in", res=300)
print(as_ggplot(p))
dev.off()




########  Supplementary figure S10: real data choice of genes ##########
rm(list=ls())
source("./code/helper.R")
#### 4 studies
results_dir <- "~/Documents/MSBE/TB_realdata_v4_10gene/"
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.min(curr_iter$value)] <- total_count[which.min(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label
plt_4s <- res_lst <- list()
i=1;j=1
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_4s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    plt_4s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=1), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}
#### 6 studies
results_dir <- "~/Documents/MSBE/TB_realdata_v3_10gene/"
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.min(curr_iter$value)] <- total_count[which.min(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label
plt_6s <- res_lst <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_6s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    plt_6s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=1), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}
lgd <- get_legend(plt_4s[["mxe"]][["GSE37250_SA"]])
plt_4studies <- ggarrange(plt_4s[["mxe"]][["India"]], plt_4s[["mxe"]][["US"]], 
                          plt_4s[["mxe"]][["GSE37250_SA"]], plt_4s[["mxe"]][["GSE37250_M"]], 
                          ncol=4, nrow=1, common.legend=TRUE, legend="none")
plt_4studies <- annotate_figure(plt_4studies, 
                                top=text_grob("4 studies", face="bold", size=18, hjust=1, x=0.15))
plt_6studies <- ggarrange(plt_6s[["mxe"]][["India"]], plt_6s[["mxe"]][["US"]], 
                          plt_6s[["mxe"]][["GSE37250_SA"]], plt_6s[["mxe"]][["GSE37250_M"]], 
                          plt_6s[["mxe"]][["Africa"]], plt_6s[["mxe"]][["GSE39941_M"]],
                          ncol=6, nrow=1, common.legend=TRUE, legend="none")
plt_6studies <- annotate_figure(plt_6studies, 
                                top=text_grob("6 studies", face="bold", size=18, hjust=1, x=0.1))
p <- arrangeGrob(plt_4studies, plt_6studies, lgd,                              
                 ncol=3, nrow=3, layout_matrix=rbind(c(2,2,2), c(NA,NA,NA), c(1,1,3)),
                 heights=c(1,0.1,1))
plt_10gene <- as_ggplot(p)
plt_10gene <- annotate_figure(plt_10gene, top=text_grob("Top 10 genes", face="bold", size=18, hjust=1, x=0.5))


#### 4 studies
results_dir <- "~/Documents/MSBE/TB_realdata_v4_100gene/"
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.min(curr_iter$value)] <- total_count[which.min(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label
plt_4s <- res_lst <- list()
i=1;j=1
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_4s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    plt_4s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=0.6), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}
#### 6 studies
results_dir <- "~/Documents/MSBE/TB_realdata_v3_100gene/"
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.min(curr_iter$value)] <- total_count[which.min(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label
plt_6s <- res_lst <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_6s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    plt_6s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=1), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}
lgd <- get_legend(plt_4s[["mxe"]][["GSE37250_SA"]])
plt_4studies <- ggarrange(plt_4s[["mxe"]][["India"]], plt_4s[["mxe"]][["US"]], 
                          plt_4s[["mxe"]][["GSE37250_SA"]], plt_4s[["mxe"]][["GSE37250_M"]], 
                          ncol=4, nrow=1, common.legend=TRUE, legend="none")
plt_4studies <- annotate_figure(plt_4studies, 
                                top=text_grob("4 studies", face="bold", size=18, hjust=1, x=0.15))
plt_6studies <- ggarrange(plt_6s[["mxe"]][["India"]], plt_6s[["mxe"]][["US"]], 
                          plt_6s[["mxe"]][["GSE37250_SA"]], plt_6s[["mxe"]][["GSE37250_M"]], 
                          plt_6s[["mxe"]][["Africa"]], plt_6s[["mxe"]][["GSE39941_M"]],
                          ncol=6, nrow=1, common.legend=TRUE, legend="none")
plt_6studies <- annotate_figure(plt_6studies, 
                                top=text_grob("6 studies", face="bold", size=18, hjust=1, x=0.1))
p <- arrangeGrob(plt_4studies, plt_6studies, lgd,                              
                 ncol=3, nrow=3, layout_matrix=rbind(c(2,2,2), c(NA,NA,NA), c(1,1,3)),
                 heights=c(1,0.1,1))
plt_100gene <- as_ggplot(p)
plt_100gene <- annotate_figure(plt_100gene, top=text_grob("Top 100 genes", face="bold", size=18, hjust=1, x=0.5))



#### 4 studies
results_dir <- "~/Documents/MSBE/TB_realdata_v4_10000gene/"
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.min(curr_iter$value)] <- total_count[which.min(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label
plt_4s <- res_lst <- list()
i=1;j=1
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_4s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    plt_4s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=0.6), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}
#### 6 studies
results_dir <- "~/Documents/MSBE/TB_realdata_v3_10000gene/"
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  total_count <- rep(0, 5)
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    total_count[which.min(curr_iter$value)] <- total_count[which.min(curr_iter$value)] + 1 
  }
  names(total_count) <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label
plt_6s <- res_lst <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_6s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    sumstats <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    plt_6s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=1), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}
lgd <- get_legend(plt_4s[["mxe"]][["GSE37250_SA"]])
plt_4studies <- ggarrange(plt_4s[["mxe"]][["India"]], plt_4s[["mxe"]][["US"]], 
                          plt_4s[["mxe"]][["GSE37250_SA"]], plt_4s[["mxe"]][["GSE37250_M"]], 
                          ncol=4, nrow=1, common.legend=TRUE, legend="none")
plt_4studies <- annotate_figure(plt_4studies, 
                                top=text_grob("4 studies", face="bold", size=18, hjust=1, x=0.15))
plt_6studies <- ggarrange(plt_6s[["mxe"]][["India"]], plt_6s[["mxe"]][["US"]], 
                          plt_6s[["mxe"]][["GSE37250_SA"]], plt_6s[["mxe"]][["GSE37250_M"]], 
                          plt_6s[["mxe"]][["Africa"]], plt_6s[["mxe"]][["GSE39941_M"]],
                          ncol=6, nrow=1, common.legend=TRUE, legend="none")
plt_6studies <- annotate_figure(plt_6studies, 
                                top=text_grob("6 studies", face="bold", size=18, hjust=1, x=0.1))
p <- arrangeGrob(plt_4studies, plt_6studies, lgd,                              
                 ncol=3, nrow=3, layout_matrix=rbind(c(2,2,2), c(NA,NA,NA), c(1,1,3)),
                 heights=c(1,0.1,1))
plt_10000gene <- as_ggplot(p)
plt_10000gene <- annotate_figure(plt_10000gene, top=text_grob("Top 10000 genes", face="bold", size=18, hjust=1, x=0.5))

png("./figures/FigS10.png", width=13, height=22, units="in", res=300)
ggarrange(plt_10gene, NA, plt_100gene, NA, plt_10000gene, nrow=5, heights=c(1,0.1,1,0.1,1))
dev.off()




########  Supplementary Table S3: biomarker difference   ########
rm(list=ls())
load("./data/TB_real_data.RData")
source("./code/helper.R")
n_highvar_genes <- 1000
study_names <- names(dat_lst)
study_names
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
s = "Africa"
sel_genes_lst <- up_genes_lst <- down_genes_lst <- list()
for(s in study_names){
  ## Get training & test set
  test_name <- s
  train_name <- setdiff(study_names, test_name)
  
  # training set
  dat <- do.call(cbind, dat_lst[train_name])
  batch <- rep(train_name, times=sapply(dat_lst[train_name], ncol))
  batches_ind <- lapply(1:length(train_name), function(i){which(batch==i)})
  batch_names <- levels(factor(batch))
  group <- do.call(c, label_lst[train_name])
  y_sgbatch_train <- lapply(batch_names, function(k){group[batch==k]})
  
  # test
  dat_test <- dat_lst[[test_name]]
  group_test <- label_lst[[test_name]]
  
  ## feature reduction - select highly variable genes in training data
  genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_highvar_genes]
  dat <- dat[genes_sel_names, ]
  dat_test <- dat_test[genes_sel_names, ]
  
  dat_batch_whole_norm <- normalizeData(dat)  # norm training set as a whole
  # norm training set within each batch
  dat_batch_norm <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))  
  for(k in batch_names){dat_batch_norm[, batch==k] <- normalizeData(dat[, batch==k])}
  
  ## select top DE genes in training & plot their expression
  group <- factor(group)
  dsg <- model.matrix(~0+group)
  ff <- lmFit(dat_batch_norm, dsg)
  contrast_matrix <- makeContrasts(group1-group0, levels=dsg)
  ff <- contrasts.fit(ff, contrast_matrix)
  ff <- eBayes(ff)
  res <- topTable(ff, coef=1, adjust="BH", number=50)
  sel_genes_lst[[s]] <- rownames(res)
  up_genes_lst[[s]] <- rownames(res)[res$logFC>0]
  down_genes_lst[[s]] <- rownames(res)[res$logFC<0]
}

sel_up_genes <- Reduce(intersect, up_genes_lst)
sel_down_genes <- Reduce(intersect, down_genes_lst)

dat_pooled <- do.call(cbind, dat_lst)
batch_pooled <- rep(names(dat_lst), times=sapply(dat_lst, ncol))
dat_norm <- matrix(NA, nrow=nrow(dat_pooled), ncol=ncol(dat_pooled), dimnames=dimnames(dat_pooled))  
for(k in study_names){dat_norm[, batch_pooled==k] <- normalizeData(dat_pooled[, batch_pooled==k])}
label_pooled <- do.call(c, label_lst)

mean_res <- lapply(study_names, function(s){
  c(up_LTBI=mean(dat_norm[sel_up_genes, batch_pooled==s & label_pooled==0]),
    up_TB=mean(dat_norm[sel_up_genes, batch_pooled==s & label_pooled==1]),
    down_LTBI=mean(dat_norm[sel_down_genes, batch_pooled==s & label_pooled==0]),
    down_TB=mean(dat_norm[sel_down_genes, batch_pooled==s & label_pooled==1]))
})
names(mean_res) <- study_names
mean_res <- t(round(do.call(rbind, mean_res)[c("India","US","GSE37250_SA","GSE37250_M","Africa","GSE39941_M"),], 3))
colnames(mean_res) <- c("D", "E", "F", "G", "A", "C")
mean_res




########  Supplementary figure S11: biomarker difference   ########
plt_hist <- function(curr_dat){
  ggplot(curr_dat, aes(x=Expression, y=..density.., fill=Condition)) +
    geom_histogram(bins=30, alpha=0.3) +
    geom_vline(data=ddply(curr_dat, .(Condition, Type), summarize, M=mean(Expression)),
               aes(xintercept=M, color=Condition), linetype="dashed") +
    facet_wrap(~Type, ncol=2) +
    scale_fill_manual(name="Condition", values=c("blue", "orange")) +
    scale_color_manual(values=c("blue", "orange")) +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          legend.position="bottom") +
    xlim(-3, 3)
}
#s="Africa"
histogram_list <- list()
for(s in study_names){
  df_up <- melt(list(LTBI=c(dat_norm[sel_up_genes, batch_pooled==s & label_pooled==0]),
                     TB=c(dat_norm[sel_up_genes, batch_pooled==s & label_pooled==1])))
  df_down <- melt(list(LTBI=c(dat_norm[sel_down_genes, batch_pooled==s & label_pooled==0]),
                       TB=c(dat_norm[sel_down_genes, batch_pooled==s & label_pooled==1])))
  colnames(df_up) <- colnames(df_down) <- c("Expression", "Condition")
  df_up$Type <- "Up regulated genes"; df_down$Type <- "Down regulated genes"
  df <- rbind(df_up, df_down)
  
  fig <- plt_hist(df)
  histogram_list[[s]] <- fig
}

lgd <- get_legend(fig)
histogram_list <- lapply(study_names, function(s){
  ht <- histogram_list[[s]] + theme(legend.position="none")
  ht <- annotate_figure(ht, top=text_grob(study_label[s], size=15, face="bold"))
  return(ht)
})
names(histogram_list) <- study_names

png("./figures/FigS11.png", width=10, height=12, units="in", res=300)
grid.arrange(histogram_list[["India"]], histogram_list[["US"]],
             histogram_list[["GSE37250_SA"]], histogram_list[["GSE37250_M"]],
             histogram_list[["Africa"]], histogram_list[["GSE39941_M"]],
             lgd, 
             nrow=4, ncol=2, widths=c(2,2), heights=c(2,2,2,0.5),  
             layout_matrix=rbind(c(1, 2), c(3, 4), c(5, 6), c(NA, 7)))
dev.off()




########  Supplementary figure S12: 4-study mxe and weights (rf only)   ########
rm(list=ls())
source("./code/helper.R")
results_dir <- "./results_real_4studies/"
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
learner_types <- c("lasso", "rf", "svm", "crossmod")

plt <- function(weights_df){
  mlt_df <- melt(weights_df)
  colnames(mlt_df)[1:2] <- c("Type", "Study")
  mlt_df$Study <- factor(mlt_df$Study, levels=c("D", "E", "F", "G"))
  colors <- c("D"="#56B4E9", "E"="#009E73", "F"="#F0E442", "G"="#0072B2")
  p <- ggplot(mlt_df, aes(x=Type, y=value, fill=Study)) +
    geom_bar(stat="identity") + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1)) +
    scale_fill_manual(values = colors) +
    labs(y="Weights")
  return(p)
}

STUDY <- "GSE37250_SA"
fig_lst <- list()
for(STUDY in study_names){
  other_study_names <- setdiff(study_names, STUDY)
  
  #i=2
  plt_lst <- res_lst <- list()
  for(i in 1:length(perf_measures)){
    curr_perf <- perf_measures[i]
    
    file_prefix <- sprintf("test%s", STUDY)
    res_lst[[curr_perf]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "ComBat"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "ComBat", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    for(MODEL in learner_types){
      sumstats_curr <- dplyr::filter(sumstats, Model==MODEL)
      plt_lst[[curr_perf]][[MODEL]] <- ggplot(sumstats_curr, aes(x=Method, y=Avg, color=Type)) +
        geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
        geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
        geom_point() +
        labs(y=perf_measures_names[perf_measures[i]],
             title=paste("Test set:", study_label[STUDY])) +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(angle=30, hjust=1),
              legend.title=element_blank())
    }
  }
  
  load(sprintf("./results_real_4studies/test%s_weights.RData", STUDY))
  
  lasso_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                                `Cross-study weights`=cs_weights_seq$lasso, 
                                `Stacked regression weights`=reg_s_beta$lasso))
  rf_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                             `Cross-study weights`=cs_weights_seq$rf, 
                             `Stacked regression weights`=reg_s_beta$rf))
  svm_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                              `Cross-study weights`=cs_weights_seq$svm, 
                              `Stacked regression weights`=reg_s_beta$svm))
  colnames(lasso_weights) <- colnames(rf_weights) <- colnames(svm_weights) <- study_label[other_study_names]
  
  lasso_w_plt <- plt(lasso_weights)
  rf_w_plt <- plt(rf_weights)
  svm_w_plt <- plt(svm_weights)
  
  curr_fig <- ggarrange(plt_lst[["mxe"]][["rf"]], rf_w_plt, ncol=2, nrow=1)
  annotate_figure(curr_fig, top = text_grob(paste("Test set:", gsub("_1", "", STUDY))))
  fig_lst[[STUDY]] <- curr_fig
}

png("./figures/FigS12.png", width=8, height=13, units="in", res=300)
ggarrange(fig_lst[["India"]], fig_lst[["US"]], fig_lst[["GSE37250_SA"]], fig_lst[["GSE37250_M"]], nrow=4, ncol=1)
dev.off()




########  Supplementary figure S13: 4-study mxe and weights (rf only) - upsampling India ##########
rm(list=ls())
source("./code/helper.R")
results_dir <- "~/Documents/MSBE/TB_realdata_v4_upsample"
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
learner_types <- c("lasso", "rf", "svm", "crossmod")

plt <- function(weights_df){
  mlt_df <- melt(weights_df)
  colnames(mlt_df)[1:2] <- c("Type", "Study")
  mlt_df$Study <- factor(mlt_df$Study, levels=c("D", "E", "F", "G"))
  colors <- c("D"="#56B4E9", "E"="#009E73", "F"="#F0E442", "G"="#0072B2")
  p <- ggplot(mlt_df, aes(x=Type, y=value, fill=Study)) +
    geom_bar(stat="identity") + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1)) +
    scale_fill_manual(values = colors) +
    labs(y="Weights")
  return(p)
}

STUDY <- "GSE37250_SA"
fig_lst <- list()
for(STUDY in study_names){
  other_study_names <- setdiff(study_names, STUDY)
  
  #i=2
  plt_lst <- res_lst <- list()
  for(i in 1:length(perf_measures)){
    curr_perf <- perf_measures[i]
    
    file_prefix <- sprintf("test%s", STUDY)
    res_lst[[curr_perf]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "ComBat"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "ComBat", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    for(MODEL in learner_types){
      sumstats_curr <- dplyr::filter(sumstats, Model==MODEL)
      plt_lst[[curr_perf]][[MODEL]] <- ggplot(sumstats_curr, aes(x=Method, y=Avg, color=Type)) +
        geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
        geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
        geom_point() +
        labs(y=perf_measures_names[perf_measures[i]],
             title=paste("Test set:", study_label[STUDY])) +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(angle=30, hjust=1),
              legend.title=element_blank())
    }
  }
  
  load(sprintf("~/Documents/MSBE/TB_realdata_v4_upsample/test%s_weights.RData", STUDY))
  
  lasso_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                                `Cross-study weights`=cs_weights_seq$lasso, 
                                `Stacked regression weights`=reg_s_beta$lasso))
  rf_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                             `Cross-study weights`=cs_weights_seq$rf, 
                             `Stacked regression weights`=reg_s_beta$rf))
  svm_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                              `Cross-study weights`=cs_weights_seq$svm, 
                              `Stacked regression weights`=reg_s_beta$svm))
  colnames(lasso_weights) <- colnames(rf_weights) <- colnames(svm_weights) <- study_label[other_study_names]
  
  lasso_w_plt <- plt(lasso_weights)
  rf_w_plt <- plt(rf_weights)
  svm_w_plt <- plt(svm_weights)
  
  curr_fig <- ggarrange(plt_lst[["mxe"]][["rf"]], rf_w_plt, ncol=2, nrow=1)
  annotate_figure(curr_fig, top = text_grob(paste("Test set:", gsub("_1", "", STUDY))))
  fig_lst[[STUDY]] <- curr_fig
}

png("./figures/FigS13.png", width=8, height=13, units="in", res=300)
ggarrange(fig_lst[["India"]], fig_lst[["US"]], 
          fig_lst[["GSE37250_SA"]], fig_lst[["GSE37250_M"]],
          nrow=4, ncol=1)
dev.off()




########  Supplementary figure S14: 4-study mxe and weights (rf only) - downsampling Africa ##########
rm(list=ls())
source("./code/helper.R")
results_dir <- "~/Documents/MSBE/TB_realdata_v4_downsample/"
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
learner_types <- c("lasso", "rf", "svm", "crossmod")

plt <- function(weights_df){
  mlt_df <- melt(weights_df)
  colnames(mlt_df)[1:2] <- c("Type", "Study")
  mlt_df$Study <- factor(mlt_df$Study, levels=c("A", "C", "D", "E", "F", "G"))
  colors <- c("A"="#999999", "C"="#E69F00", "D"="#56B4E9", 
              "E"="#009E73", "F"="#F0E442", "G"="#0072B2")
  p <- ggplot(mlt_df, aes(x=Type, y=value, fill=Study)) +
    geom_bar(stat="identity") + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=30, hjust=1)) +
    scale_fill_manual(values = colors) +
    labs(y="Weights")
  return(p)
}

STUDY <- "GSE37250_SA"
fig_lst <- list()
for(STUDY in study_names){
  other_study_names <- setdiff(study_names, STUDY)
  
  #i=2
  plt_lst <- res_lst <- list()
  for(i in 1:length(perf_measures)){
    curr_perf <- perf_measures[i]
    
    file_prefix <- sprintf("test%s", STUDY)
    res_lst[[curr_perf]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]]) <- c("Method", "value", "Model", "Iteration")
    
    sumstats <- res_lst[[curr_perf]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")) %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025))
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method=="ComBat"] <- "ComBat"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "ComBat", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    #STUDY_NAME <- ifelse(STUDY=="GSE39941_M", "GSE39941", STUDY)
    for(MODEL in learner_types){
      sumstats_curr <- dplyr::filter(sumstats, Model==MODEL)
      plt_lst[[curr_perf]][[MODEL]] <- ggplot(sumstats_curr, aes(x=Method, y=Avg, color=Type)) +
        geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
        geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
        geom_point() +
        labs(y=perf_measures_names[perf_measures[i]],
             title=paste("Test set:", study_label[STUDY])) +
        #title=ifelse(curr_perf=="auc", paste("Learner:", MODEL), "")) +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_text(angle=30, hjust=1),
              legend.title=element_blank())
    }
  }
  
  load(sprintf("~/Documents/MSBE/TB_realdata_v4_downsample/test%s_weights.RData", STUDY))
  
  lasso_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                                `Cross-study weights`=cs_weights_seq$lasso, 
                                `Stacked regression weights`=reg_s_beta$lasso))
  rf_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                             `Cross-study weights`=cs_weights_seq$rf, 
                             `Stacked regression weights`=reg_s_beta$rf))
  svm_weights <- t(data.frame(`Batch-size weights`=navg_weights,
                              `Cross-study weights`=cs_weights_seq$svm, 
                              `Stacked regression weights`=reg_s_beta$svm))
  colnames(lasso_weights) <- colnames(rf_weights) <- colnames(svm_weights) <- study_label[other_study_names]
  
  lasso_w_plt <- plt(lasso_weights)
  rf_w_plt <- plt(rf_weights)
  svm_w_plt <- plt(svm_weights)
  
  curr_fig <- ggarrange(plt_lst[["mxe"]][["rf"]], rf_w_plt, ncol=2, nrow=1)
  annotate_figure(curr_fig, top = text_grob(paste("Test set:", gsub("_1", "", STUDY))))
  fig_lst[[STUDY]] <- curr_fig
}

png("./figures/FigS14.png", width=7, height=15, units="in", res=300)
ggarrange(fig_lst[["Africa"]], fig_lst[["GSE39941_M"]], 
          fig_lst[["India"]], fig_lst[["US"]], 
          fig_lst[["GSE37250_SA"]], fig_lst[["GSE37250_M"]], 
          nrow=6, ncol=1)
dev.off()

