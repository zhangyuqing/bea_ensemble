rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/AEGIS_simbatch/")
sapply(c("ggplot2", "gridExtra", "reshape"), require, character.only = TRUE)
perf_measure_name <- "mxe"
use_ref_combat <- "F"

####  Read in results
res_lasso <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_reduceF_useref%s.csv', perf_measure_name, 'lasso', use_ref_combat), header=TRUE)
#res_lasso <- res_lasso[1:500, ]
res_rf <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_reduceF_useref%s.csv', perf_measure_name, 'rf', use_ref_combat), header=TRUE)
res_nnet <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_reduceF_useref%s.csv', perf_measure_name, 'nnet', use_ref_combat), header=TRUE)
res_svm <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_reduceF_useref%s.csv', perf_measure_name, 'svm', use_ref_combat), header=TRUE)


####  Multiple methods
res_allmethods <- data.frame(rbind(res_lasso, res_rf, res_nnet, res_svm),
                             Method=c(rep("Lasso",nrow(res_lasso)), rep("Random Forest",nrow(res_rf)),
                                      rep("Neural Nets",nrow(res_nnet)), rep("SVM",nrow(res_svm))))
res_mlt <- melt(res_allmethods)
png(sprintf('./figures/cmpMod_boxplot_%s_useref%s_whole.png', perf_measure_name, use_ref_combat), width=6, height=8, units="in", res=300)
print(ggplot(res_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss")) 
dev.off()  

res_subset <- res_allmethods[, c(1,8:14)]
res_subset_mlt <- melt(res_subset)
png(sprintf('./figures/cmpMod_boxplot_%s_useref%s_subset.png', perf_measure_name, use_ref_combat), width=4, height=8, units="in", res=300)
print(ggplot(res_subset_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss", limits=c(min(res_subset[,-8])-0.02,1.5))) 
dev.off()  


