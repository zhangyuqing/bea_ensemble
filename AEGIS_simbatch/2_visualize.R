rm(list=ls())
setwd("/restricted/projectnb/combat/work/yuqingz/multistudy_batcheffect/AEGIS_simbatch/")
sapply(c("ggplot2", "gridExtra", "reshape"), require, character.only = TRUE)
learner_type <- "Lasso"
perf_measure_name <- "mxe"
res <- read.csv(sprintf('batchCSL_AEGIS_%s_%s.csv', perf_measure_name, tolower(learner_type)), header=TRUE)


####  Boxplot of prediction accuracy
## All panels
#boxplot
res_mlt <- melt(res)
png(sprintf('%s_%s_boxplot_whole.png', perf_measure_name, learner_type), width=8, height=5, units="in", res=300)
print(ggplot(res_mlt, aes(x=variable, y=value)) +
  geom_boxplot() +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Mean cross-entropy loss") +
  ggtitle(learner_type))
dev.off()  

# line of median accuracy
res_avg <- apply(res, 2, mean)
res_avg_df <- data.frame(method=factor(names(res_avg), levels=names(res_avg)), 
                         value=as.numeric(res_avg))
png(sprintf('%s_%s_line_whole.png', perf_measure_name, learner_type), width=8, height=5, units="in", res=300)
print(ggplot(data=res_avg_df, aes(x=method, y=value, group=1)) +
  geom_line()+
  geom_point() +
  labs(title=learner_type, x="", y="Mean cross-entropy loss")) 
dev.off()


## Subset of panels
res_subset <- res[, c(1,8:13)]
#boxplot
res_subset_mlt <- melt(res_subset)
png(sprintf('%s_%s_boxplot_subset.png', perf_measure_name, learner_type), width=5, height=5, units="in", res=300)
print(ggplot(res_subset_mlt, aes(x=variable, y=value)) +
        geom_boxplot() +
        scale_x_discrete(name="") +
        scale_y_continuous(name="Mean cross-entropy loss", limits=c(min(res_subset)-0.02,1.5)) +
        ggtitle(learner_type))
dev.off()  

# line of median accuracy
res_subset_avg <- apply(res_subset, 2, mean)
res_subset_avg_df <- data.frame(method=factor(names(res_subset_avg), levels=names(res_subset_avg)), 
                                value=as.numeric(res_subset_avg))
png(sprintf('%s_%s_line_subset.png', perf_measure_name, learner_type), width=5, height=5, units="in", res=300)
print(ggplot(data=res_subset_avg_df, aes(x=method, y=value, group=1)) +
        geom_line()+
        geom_point() +
        labs(title=learner_type, x="", y="Mean cross-entropy loss")) 
dev.off()



####  Correlation of prediction scores
load("test_pred_scores.RData")
corr_mat_lst <- lapply(pred_mat_lst, function(pred_mat){cor(pred_mat, method="spearman")})
corr_mat_include <- sapply(corr_mat_lst, function(cor_mat){!any(is.na(cor_mat))})
corr_mat_lst <- corr_mat_lst[corr_mat_include]
corr_mat_avg <- Reduce("+", corr_mat_lst) / length(corr_mat_lst)
corr_mat_avg <- round(corr_mat_avg, 2)

corr_mat_avg_mlt <- melt(corr_mat_avg)
corr_mat_avg_mlt$X1 <- factor(corr_mat_avg_mlt$X1, levels=paste0("Batch",5:1))
png(sprintf('%s_%s_corr_heatmap.png', perf_measure_name, learner_type), width=5, height=4.5, units="in", res=300)
ggplot(corr_mat_avg_mlt, aes(X2, X1)) + 
  geom_tile(aes(fill=value), colour="white") + 
  geom_text(aes(X2, X1, label=value), color="black", size=4) +
  scale_fill_gradient(low="white", high="red", limit=c(0,1)) +
  labs(title="Average Spearman correlation of predictions", x="", y="") +
  theme_minimal() +
  coord_fixed()
dev.off()
