rm(list=ls())
setwd("~/Documents/MSBE/TB/")
sapply(c("ggplot2", "gridExtra", "reshape2"), require, character.only=TRUE)

if(!dir.exists("figures")){dir.create("figures")}
file_lst <- dir()[grep(".csv", dir(), fixed=TRUE)]  # all results files

method_names <- c("lasso", "rf")
Nbatch <- 5
batch_mean_vec <- c(2, 4)  #c(0, 2, 4)
batch_var_vec <- c(1, 2, 4)
subset_colnames <- c("NoBatch", "Batch", "ComBat", "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s")
perf_measures <- c("mxe", "auc", "acc", "f")
perf_measures_plt <- c("Mean cross-entropy loss", "AUC", "Accuracy", "F1 score")
names(perf_measures_plt) <- perf_measures

#curr_mod <- "rf"; curr_perf <- "mxe"
for(curr_mod in method_names){
  for(curr_perf in perf_measures){
    curr_file_lst <- file_lst[grepl(paste0('_',curr_mod,'_'), file_lst) & grepl(paste0("^",curr_perf), file_lst)]
    curr_file_lst <- curr_file_lst[4:9]
    print(curr_file_lst)
    res_lst <- lapply(curr_file_lst, function(curr_file){read.csv(curr_file, header=TRUE)[, subset_colnames]})
    names(res_lst) <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)), 1, paste, collapse="_"))
    
    png(sprintf("./figures/%s_%s.png", curr_mod, curr_perf), width=8, height=6, units="in", res=300)
    print(ggplot(melt(res_lst), aes(x=variable, y=value)) + #, color=variable)) +
            geom_boxplot() +
            facet_wrap(~L1, ncol=3) +
            theme(axis.text.x=element_text(angle=45, hjust=1),
                  axis.title.x=element_blank()) +
            labs(y=perf_measures_plt[curr_perf])) 
    dev.off()  
  }
}



