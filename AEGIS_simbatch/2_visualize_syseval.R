rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/AEGIS_simbatch/SysEval")
sapply(c("ggplot2", "gridExtra", "reshape2"), require, character.only=TRUE)

if(!dir.exists("figures")){dir.create("figures")}
file_lst <- dir()[grep(".csv", dir(), fixed=TRUE)]  # all results files
method_vec <- sort(c("lasso", "rf", "nnet", "svm"))
Nbatch_vec <- c(2, 3, 4, 5)
batch_mean_vec <- c(0, 1, 2, 3)
batch_var_vec <- c(1, 2, 4, 8)
Nsize_vec <- c(30, 40, 50, 60, 70)
subset_colnames <- c("NoBatch", "Batch", "ComBat", 
                     "Avg", "n_Avg", "CS_Avg", "Reg_a", "Reg_s", "Method")


####  Number of batches  ####
nbatch_file_lst <- file_lst[grep("nbatch", file_lst)]
print(length(nbatch_file_lst))

nbatch_res_lst <- nbatch_subset_res_lst <- list()
for(i in seq_along(Nbatch_vec)){
  # read in result files for this N_batch level
  curr_nbatch_files <- nbatch_file_lst[grep(paste0("batchN",Nbatch_vec[i]), nbatch_file_lst)]
  res_lst <- lapply(curr_nbatch_files, function(curr_file){
    read.csv(curr_file, header=TRUE)
  })
  
  res_allmethods <- data.frame(do.call(rbind, res_lst), Method=rep(method_vec, sapply(res_lst, nrow)))
  res_subset <- res_allmethods[, subset_colnames]
  
  nbatch_res_lst[[i]] <- melt(res_allmethods)
  nbatch_subset_res_lst[[i]] <- melt(res_subset)
}

res_all <- data.frame(do.call(rbind, nbatch_res_lst), 
                      Nbatch=rep(paste("N.batch", Nbatch_vec), sapply(nbatch_res_lst, nrow)))
re_level <- levels(res_all$variable)[c(1:2,5:10,3:4,11:13)]
res_all$variable <- factor(res_all$variable, levels=re_level)
png("./figures/nbatch_whole.png", width=12, height=9, units="in", res=300)
print(ggplot(res_all, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ Nbatch) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x=element_blank()) +
        scale_y_continuous(name="Mean cross-entropy loss", 
                           limits=c(min(res_allmethods[,-ncol(res_allmethods)])-0.02,5)) +
        labs(title="Number of batches")) 
dev.off()  

res_all_subset <- data.frame(do.call(rbind, nbatch_subset_res_lst), 
                             Nbatch=rep(paste("N.batch", Nbatch_vec), sapply(nbatch_subset_res_lst, nrow)))
png("./figures/nbatch_subset.png", width=10, height=9, units="in", res=300)
print(ggplot(res_all_subset, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ Nbatch) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x=element_blank()) +
        scale_y_continuous(name="Mean cross-entropy loss", 
                           limits=c(min(res_subset[,-ncol(res_subset)])-0.02,5))+
        labs(title="Number of batches")) 
dev.off()  
rm(curr_nbatch_files, nbatch_file_lst, i, nbatch_res_lst, nbatch_subset_res_lst, 
   re_level, res_all, res_all_subset, res_allmethods, res_subset, res_lst)


####  Sample size  ####
size_file_lst <- file_lst[grep("^size", file_lst)]
print(length(size_file_lst))

size_res_lst <- size_subset_res_lst <- list()
for(j in seq_along(Nsize_vec)){
  # read in result files for this N_size level
  curr_size_files <- size_file_lst[grep(paste0("size", Nsize_vec[j]), size_file_lst)]
  res_lst <- lapply(curr_size_files, function(curr_file){
    read.csv(curr_file, header=TRUE)
  })
  
  res_allmethods <- data.frame(do.call(rbind, res_lst), Method=rep(method_vec, sapply(res_lst, nrow)))
  res_subset <- res_allmethods[, subset_colnames]
  
  size_res_lst[[j]] <- melt(res_allmethods)
  size_subset_res_lst[[j]] <- melt(res_subset)
}

res_all <- data.frame(do.call(rbind, size_res_lst), 
                      Nsize=rep(paste("N.sample", Nsize_vec*2), sapply(size_res_lst, nrow)))
re_level <- levels(res_all$variable)[c(1:2,6:11,3:5)]
res_all$variable <- factor(res_all$variable, levels=re_level)
res_all$Nsize <- factor(res_all$Nsize, levels=rep(paste("N.sample", Nsize_vec*2)))
                        
png("./figures/size_whole.png", width=14, height=9, units="in", res=300)
print(ggplot(res_all, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ Nsize) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x=element_blank()) +
        scale_y_continuous(name="Mean cross-entropy loss", 
                           limits=c(min(res_allmethods[,-ncol(res_allmethods)])-0.02,5)) +
        labs(title="Sample size")) 
dev.off()  

res_all_subset <- data.frame(do.call(rbind, size_subset_res_lst), 
                             Nsize=rep(paste("N.sample", Nsize_vec*2), sapply(size_subset_res_lst, nrow)))
res_all_subset$Nsize <- factor(res_all_subset$Nsize, levels=rep(paste("N.sample", Nsize_vec*2)))

png("./figures/size_subset.png", width=12, height=9, units="in", res=300)
print(ggplot(res_all_subset, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ Nsize) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x=element_blank()) +
        scale_y_continuous(name="Mean cross-entropy loss", 
                           limits=c(min(res_subset[,-ncol(res_subset)])-0.02,5))+
        labs(title="Sample size")) 
dev.off()  

rm(curr_size_files, size_file_lst, j, size_res_lst, size_subset_res_lst, 
   re_level, res_all, res_all_subset, res_allmethods, res_subset, res_lst)


####  Strength of signal: mean batch effect  ####




