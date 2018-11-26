rm(list=ls())
setwd("~/Dropbox/Work/MultiStudy_BatchEffect/AEGIS_simbatch/")
sapply(c("GEOquery", "ggplot2", "limma"), require, character.only=TRUE)
source("helper.R")

####  Load data from GEO
aegis_obj <- getGEO(GEO="GSE66499")[[1]]
table(substr(aegis_obj$title,1,2))

####  AEGIS-1
aegis1_data <- aegis_obj[, substr(aegis_obj$title,1,2)=="A1"]
# keep samples with unique identifies
unique_title <- names(table(substr(aegis1_data$title,1,7)))[table(substr(aegis1_data$title,1,7))==1]
aegis1_data <- aegis1_data[, substr(aegis1_data$title,1,7) %in% unique_title]
#which(table(substr(aegis1_data$title,1,7))!=1)
condition_1 <- pData(aegis1_data)[, "cancer_yes_no:ch1"]
condition_1 <- ifelse(condition_1=="Y", "Cancer", "Normal")
table(condition_1)
# PCA
pca_res_1 <- prcomp(t(exprs(aegis1_data)), center=TRUE, scale.=TRUE)
plt_obj_1 <- data.frame(PC1=pca_res_1$x[, 1], PC2=pca_res_1$x[, 2], Condition=condition_1)
png("PCAplot_AEGIS1.png", width=6, height=5, units="in", res=300)
ggplot(plt_obj_1, aes(x=PC1, y=PC2, color=Condition)) + 
  geom_point() +
  scale_color_manual(values=c('#E69F00', '#56B4E9')) +
  ggtitle("AEGIS-1 samples")
dev.off()


####  AEGIS-2
aegis2_data <- aegis_obj[, substr(aegis_obj$title,1,2)=="A2"]
condition_2 <- pData(aegis2_data)[, "cancer_yes_no:ch1"]
condition_2 <- ifelse(condition_2=="Y", "Cancer", "Normal")
table(condition_2)
### 267 cancer, 74 normal, lung cancer prevalence 78%
### Data processing: RMA normalized, all run in a single microarray batch; identify outliers
# PCA
pca_res_2 <- prcomp(t(exprs(aegis2_data)), center=TRUE, scale.=TRUE)
plt_obj_2 <- data.frame(PC1=pca_res_2$x[, 1], PC2=pca_res_2$x[, 2], Condition=condition_2)
png("PCAplot_AEGIS2.png", width=6, height=5, units="in", res=300)
ggplot(plt_obj_2, aes(x=PC1, y=PC2, color=Condition)) + 
  geom_point() +
  scale_color_manual(values=c('#E69F00', '#56B4E9')) +
  ggtitle("AEGIS-2 samples")
dev.off()


####  Pooled PCA
pca_res_cmb <- prcomp(rbind(t(exprs(aegis1_data)), t(exprs(aegis2_data))), center=TRUE, scale.=TRUE)
condition_cmb <- c(condition_1, condition_2)
batch <- c(rep("AEGIS-1",ncol(aegis1_data)), rep("AEGIS-2",ncol(aegis2_data)))
plt_obj_cmb <- data.frame(PC1=pca_res_cmb$x[, 1], PC2=pca_res_cmb$x[, 2], 
                          Condition=condition_cmb, Batch=batch)
# color by condition, marker shape by batch
png("PCAplot_combined.png", width=6, height=5, units="in", res=300)
ggplot(plt_obj_cmb, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=Condition, shape=Batch)) +
  scale_shape_manual(values=c(16, 17))+
  scale_color_manual(values=c('#E69F00', '#56B4E9')) +
  ggtitle("Combined AEGIS 1 and 2")
dev.off()


####  Feature reduction: DE with limma
de_fit <- lmFit(exprs(aegis2_data), design=model.matrix(~factor(condition_2, levels=c("Normal", "Cancer"))))
de_fit <- eBayes(de_fit)
selected_genes <- rownames(topTable(de_fit, number=100))


####  Construct hyper pars for batch effect
## For mean batch effect:
# what is the range of biological signal (fold change of genes)?
de_res <- topTable(de_fit, number=100)
print(round(range(exp(de_res$logFC)), 3))
# what is the range of gene-wise mean?
sel_dat <- exprs(aegis2_data)[selected_genes, ]
mean_seq <- apply(sel_dat, 1, mean)
print(median(mean_seq)); print(range(mean_seq))

# for variance batch effect - what are variances of genes in one condition?


####  Save data
test_expr_sel <- exprs(aegis1_data)[selected_genes, ]
train_expr_sel <- exprs(aegis2_data)[selected_genes, ]

test_expr_whole <- exprs(aegis1_data)
train_expr_whole <- exprs(aegis2_data)

y_test <- rep(0, length(condition_1)); y_test[condition_1=="Cancer"] <- 1; y_test=as.factor(y_test)
y_train <- rep(0, length(condition_2)); y_train[condition_2=="Cancer"] <- 1; y_train=as.factor(y_train)

save(train_expr_sel, test_expr_sel, train_expr_whole, test_expr_whole, y_train, y_test, 
     file="AEGIS_data_new.RData")
