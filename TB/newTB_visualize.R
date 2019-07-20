rm(list=ls())
setwd("~/Documents/MSBE/newTB/")
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", 
         "dplyr", "ggpubr"), require, character.only=TRUE)

study_names <- c("Africa", "India", "US")
norm_data <- TRUE 
use_ref_combat <- FALSE 
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")

#i=1; j=1
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
                                                        "n_Avg"="Batch size weights", "CS_Avg"="Cross-study weights",
                                                        "Reg_s"="Regression stacking"))
    
    sumstats_crossmod <- dplyr::filter(sumstats, Model=="crossmod")
    plt_crossmod[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="auc", paste("Test set:", gsub("_1", "", curr_testset)), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}

png("newTB.png", width=9, height=3.5, units="in", res=300)
ggarrange(plt_crossmod[["auc"]][["Africa"]], #plt_crossmod[["mxe"]][["Africa"]], 
          plt_crossmod[["auc"]][["India"]], #plt_crossmod[["mxe"]][["India"]], 
          plt_crossmod[["auc"]][["US"]], #plt_crossmod[["mxe"]][["US"]], 
          ncol=3, nrow=1, common.legend=TRUE, legend="bottom")
dev.off()



## calculate percentage in all bootstrap when ensemble out-performs combat
rm(list=ls())
setwd("~/Documents/MSBE/newTB/")
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "dplyr", "plyr","ggpubr"), 
       require, character.only=TRUE)

study_names <- c("Africa", "India", "US")
norm_data <- TRUE 
use_ref_combat <- FALSE 
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "n_Avg", "CS_Avg", "Reg_s")
B <- 100

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
  for(b in 1:B){
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
  freq_lst[[study_names[j]]] <- total_count / B
}

do.call(rbind, freq_lst)[, c("n_Avg", "CS_Avg", "Reg_s")]


## PCA
rm(list=ls())
setwd("~/Documents/MSBE/newTB/")
load("new_combined_unmatched.RData")
dat_combined <- t(do.call(cbind, dat_lst))
dat_combined <- dat_combined[, colVars(dat_combined)!=0]
batch <- rep(names(dat_lst), sapply(dat_lst, ncol))
  
pca_obj <- prcomp(dat_combined, center=TRUE, scale.=TRUE)
prop_var <- summary(pca_obj)$importance["Proportion of Variance", 1:2]
png("newTB_PCA.png", width=6, height=5, units="in", res=300)
ggplot(data.frame(pca_obj$x[,1:2], Batch=factor(batch)), aes(x=PC1, y=PC2, color=Batch)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", scales::percent(prop_var[1])),
       y=sprintf("PC2: %s Variance", scales::percent(prop_var[2])))
dev.off()
