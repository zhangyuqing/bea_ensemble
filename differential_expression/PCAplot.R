pca_res <- prcomp(t(train_expr_batch), center=TRUE, scale.=TRUE)
pca_plt_obj <- data.frame(PC1=pca_res$x[, 1], PC2=pca_res$x[, 2],
                          Condition=y_train, Batch=as.factor(batch))
png(sprintf("PCAplot_DE_ID%s.png", ID), width=6, height=5, units="in", res=300)
ggplot(pca_plt_obj, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=Batch, shape=Condition)) +
  scale_shape_manual(values=c(16, 17))+
  ggtitle("dataset with simulated batch effect")
dev.off()
rm(pca_res, pca_plt_obj)
