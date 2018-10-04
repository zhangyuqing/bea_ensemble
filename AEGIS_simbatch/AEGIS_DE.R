rm(list=ls())
setwd("C:/Users/zhang/Dropbox/Work/MultiStudy_BatchEffect/AEGIS_simbatch/")
lapply(c("GEOquery", "ggplot2", "limma", "BatchQC"), require, character.only=TRUE)


####  Load data from GEO
aegis_obj <- getGEO(GEO="GSE66499")[[1]]
#table(substr(aegis_obj$title,1,2))


####  AEGIS-1
aegis1_data <- aegis_obj[, substr(aegis_obj$title,1,2)=="A1"]
# keep samples with unique identifies
unique_title <- names(table(substr(aegis1_data$title,1,7)))[table(substr(aegis1_data$title,1,7))==1]
aegis1_data <- aegis1_data[, substr(aegis1_data$title,1,7) %in% unique_title]
#which(table(substr(aegis1_data$title,1,7))!=1)
condition_1 <- pData(aegis1_data)[, "cancer_yes_no:ch1"]
condition_1 <- ifelse(condition_1=="Y", "Cancer", "Normal")
table(condition_1)


####  AEGIS-2
aegis2_data <- aegis_obj[, substr(aegis_obj$title,1,2)=="A2"]
condition_2 <- pData(aegis2_data)[, "cancer_yes_no:ch1"]
condition_2 <- ifelse(condition_2=="Y", "Cancer", "Normal")
table(condition_2)


####  DE in AEGIS-2
de_fit <- lmFit(exprs(aegis2_data), design=model.matrix(~factor(condition_2, levels=c("Normal", "Cancer"))))
de_fit <- eBayes(de_fit)
de_res_aegis2 <- topTable(de_fit, number=nrow(aegis2_data))
length(which(de_res_aegis2$adj.P.Val<0.05))

de_res_aegis2_sel <- topTable(de_fit, number=100)
selected_genes <- rownames(de_res_aegis2_sel)
range(de_res_aegis2_sel$adj.P.Val)
hist(de_res_aegis2_sel$adj.P.Val, xlab="FDR adjusted P values",
     main="FDR adjusted P values in top 100 genes, AEGIS-2 training")
rm(de_fit)


####  DE in AEGIS-1
de_fit <- lmFit(exprs(aegis1_data), design=model.matrix(~factor(condition_1, levels=c("Normal", "Cancer"))))
de_fit <- eBayes(de_fit)
de_res_aegis1 <- topTable(de_fit, number=nrow(aegis1_data))
length(which(de_res_aegis1$adj.P.Val<0.05))
range(de_res_aegis1[selected_genes, "adj.P.Val"])
hist(de_res_aegis1[selected_genes, "adj.P.Val"], xlab="FDR adjusted P values",
     main="FDR adjusted P values in selected genes, AEGIS-1 test")
