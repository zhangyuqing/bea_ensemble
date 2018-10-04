rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/AEGIS_simbatch/")
sapply(c("ggplot2", "gridExtra", "reshape"), require, character.only = TRUE)
perf_measure_name <- "mxe"

####  Read in results
res_lasso <- read.csv(sprintf('batchCSL_AEGIS_%s_%s.csv', perf_measure_name, 'lasso'), header=TRUE)
res_lasso <- res_lasso[1:500, ]
res_rf <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_redFALSE.csv', perf_measure_name, 'rf'), header=TRUE)
res_nnet <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_redFALSE.csv', perf_measure_name, 'nnet'), header=TRUE)
res_svm <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_redFALSE.csv', perf_measure_name, 'svm'), header=TRUE)

res_Red_lasso <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_redTRUE.csv', perf_measure_name, 'lasso'), header=TRUE)
res_Red_rf <- read.csv(sprintf('batchCSL_AEGIS_%s_%s_redTRUE.csv', perf_measure_name, 'rf'), header=TRUE)


####  Multiple methods
res_allmethods <- data.frame(rbind(res_lasso, res_rf, res_nnet, res_svm),
                             Method=c(rep("Lasso",nrow(res_lasso)), rep("Random Forest",nrow(res_rf)),
                                      rep("Neural Nets",nrow(res_nnet)), rep("SVM",nrow(res_svm))))
res_mlt <- melt(res_allmethods)
png(sprintf('./figures/cmpMod_boxplot_%s_whole.png', perf_measure_name), width=6, height=8, units="in", res=300)
print(ggplot(res_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss")) 
dev.off()  

res_subset <- res_allmethods[, c(1,8:14)]
res_subset_mlt <- melt(res_subset)
png(sprintf('./figures/cmpMod_boxplot_%s_subset.png', perf_measure_name), width=4, height=8, units="in", res=300)
print(ggplot(res_subset_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Method ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss", limits=c(min(res_subset[,-8])-0.02,1.5))) 
dev.off()  


####  Reduce batch size
# lasso
res_size_lasso <- data.frame(rbind(res_lasso, res_Red_lasso),
                             Data=c(rep("Original", nrow(res_lasso)), rep("Reduced size", nrow(res_Red_lasso))))
res_size_lasso_mlt <- melt(res_size_lasso)
png(sprintf('./figures/cmpSize_boxplot_lasso_%s_whole.png', perf_measure_name), width=6, height=8, units="in", res=300)
print(ggplot(res_size_lasso_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Data ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss") +
        ggtitle("Lasso")) 
dev.off() 

res_size_lasso_subset <- res_size_lasso[, c(1,8:14)]
res_size_lasso_subset_mlt <- melt(res_size_lasso_subset)
png(sprintf('./figures/cmpSize_boxplot_lasso_%s_subset.png', perf_measure_name), width=4, height=8, units="in", res=300)
print(ggplot(res_size_lasso_subset_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Data ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss", limits=c(min(res_size_lasso_subset[,-8])-0.02,1.5)) +
        ggtitle("Lasso")) 
dev.off()  


# random forest
res_size_rf <- data.frame(rbind(res_rf, res_Red_rf),
                          Data=c(rep("Original", nrow(res_rf)), rep("Reduced size", nrow(res_Red_rf))))
res_size_rf_mlt <- melt(res_size_rf)
png(sprintf('./figures/cmpSize_boxplot_rf_%s_whole.png', perf_measure_name), width=6, height=8, units="in", res=300)
print(ggplot(res_size_rf_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Data ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss") +
        ggtitle("Random Forest"))
dev.off() 

res_size_rf_subset <- res_size_rf[, c(1,8:14)]
res_size_rf_subset_mlt <- melt(res_size_rf_subset)
png(sprintf('./figures/cmpSize_boxplot_rf_%s_subset.png', perf_measure_name), width=4, height=8, units="in", res=300)
print(ggplot(res_size_rf_subset_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        facet_grid(Data ~ .) +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss", limits=c(min(res_size_rf_subset[,-8])-0.02,1.5)) +
        ggtitle("Random Forest")) 
dev.off()  
