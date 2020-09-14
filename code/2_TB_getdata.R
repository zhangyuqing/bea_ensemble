rm(list=ls())
sapply(c("GEOquery", "annotate", "hugene11sttranscriptcluster.db",
         "SummarizedExperiment", "limma", "BatchQC", "ggplot2"), require, character.only=TRUE)
set.seed(123)

##  Download data from GEO
gse <- getGEO("GSE73408", GSEMatrix=TRUE)[[1]]

##  Annotation
# annotate genes
x <- hugene11sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
gene_symbols <- sapply(featureNames(gse), function(s) 
  ifelse(s %in% names(xx), xx[[s]], NA)
)
#length(gene_symbols); sum(!is.na(gene_symbols))

gse <- gse[!is.na(gene_symbols), ]
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
dat <- exprs(gse)
rownames(dat) <- gene_symbols

# annotate samples and create group
sample_info <- read.csv("./data/new_data_info.csv", as.is=TRUE)
identical(sample_info$ID, colnames(dat))
colnames(dat) <- sample_info$Label

dat <- dat[, -grep("_PNA_", colnames(dat))]
group <- rep(0, ncol(dat))
group[grep("_TB_", colnames(dat))] <- 1


##  Clean together with the other studies
rds_obj <- readRDS("./data/combined.rds")
overlapping_genes <- intersect(rownames(rds_obj), rownames(dat))
length(overlapping_genes)

dat <- dat[overlapping_genes, ]

# take india, africa studies
batch_filter <- rds_obj$SequencingBatch %in% c("Africa", "India")
group_filter <- rds_obj$Label %in% c("Non-progressor", "Active")
rds_obj <- rds_obj[overlapping_genes, batch_filter & group_filter]

# take a subset of Africa to make the data balanced
# africa_id <- which(rds_obj$SequencingBatch=="Africa")
# old_africa_group <- rds_obj$Label[africa_id]
# keep_id <- sample(which(old_africa_group=="Non-progressor"),
#                   sum(old_africa_group=="Active"), replace=FALSE)
# rm_id <- setdiff(which(old_africa_group=="Non-progressor"), keep_id)
# rds_obj <- rds_obj[, -africa_id[rm_id]]
# new_africa_group <- rds_obj$Label[rds_obj$SequencingBatch=="Africa"]
# print(sprintf("Prevalence in Africa: %s",
#               round(table(new_africa_group)['Active']/sum(table(new_africa_group)),2)))

# take a subset of India to make the data balanced
# india_id <- which(rds_obj$SequencingBatch=="India")
# old_india_group <- rds_obj$Label[india_id]
# keep_id <- sample(which(old_india_group=="Active"),
#                   sum(old_india_group=="Non-progressor"), replace=FALSE)
# rm_id <- setdiff(which(old_india_group=="Active"), keep_id)
# rds_obj <- rds_obj[, -india_id[rm_id]]
# new_india_group <- rds_obj$Label[rds_obj$SequencingBatch=="India"]
# print(sprintf("Prevalence in India: %s",
#               round(table(new_india_group)['Active']/sum(table(new_india_group)),2)))

# get africa study
rds_a <- rds_obj[, rds_obj$SequencingBatch=="Africa"]
dat_a <- assay(rds_a, "logfpkm")
group_a_org <- as.character(colData(rds_a)$Label)
group_a <- rep(0, length(group_a_org)); group_a[group_a_org=="Active"] <- 1

# get india study
rds_i <- rds_obj[, rds_obj$SequencingBatch=="India"]
dat_i <- assay(rds_i, "logfpkm")
group_i_org <- as.character(colData(rds_i)$Label)
group_i <- rep(0, length(group_i_org)); group_i[group_i_org=="Active"] <- 1


keepA <- rowVars(dat_a) > 0 & rowSums(dat_a!=0) > 2 
keepI <-  rowVars(dat_i) > 0 & rowSums(dat_i!=0) > 2 
keepUS <- rowVars(dat) > 0 & rowSums(dat!=0) > 2
dat_lst <- list(US=dat, Africa=dat_a, India=dat_i)
dat_lst <- lapply(dat_lst, function(m) return(m[keepA & keepI & keepUS, ]))
label_lst <- list(US=group, Africa=group_a, India=group_i)




##  Download new data from GEO
library("illuminaHumanv4.db")
# GSE37250
gse <- getGEO("GSE37250", GSEMatrix=TRUE)[[1]]
meta_info <- pData(gse)[, c(8, 10, 11, 12, 35, 36, 37)]

# GSE37250: South Africa
data1_ind <- meta_info$`hiv status:ch1`=="HIV negative" &
  meta_info$`geographical region:ch1`=="South Africa" & 
  meta_info$`disease state:ch1` %in% c("active tuberculosis", "latent TB infection") 
data1 <- exprs(gse[, data1_ind])
meta1 <- meta_info[data1_ind, ]
group1 <- factor(meta1$`disease state:ch1`, 
                 levels=c("latent TB infection", "active tuberculosis"))
group1 <- plyr::revalue(group1, c("latent TB infection"=0, "active tuberculosis"=1))

# GSE37250: Malawi
data2_ind <- meta_info$`hiv status:ch1`=="HIV negative" &
  meta_info$`geographical region:ch1`=="Malawi" & 
  meta_info$`disease state:ch1` %in% c("active tuberculosis", "latent TB infection") 
data2 <- exprs(gse[, data2_ind])
meta2 <- meta_info[data2_ind, ]
group2 <- factor(meta2$`disease state:ch1`, 
                 levels=c("latent TB infection", "active tuberculosis"))
group2 <- plyr::revalue(group2, c("latent TB infection"=0, "active tuberculosis"=1))

rm(gse, meta_info)

# GSE39941
gse <- getGEO("GSE39941", GSEMatrix=TRUE)[[1]]
meta_info <- pData(gse)[, c(8, 10, 11, 12, 35, 36, 37)]

# GSE39941: Malawi
data3_ind <- meta_info$`hiv status:ch1`=="HIV negative" &
  meta_info$`geographical region:ch1`=="Malawi" & 
  meta_info$`disease status:ch1` %in% c("active tuberculosis", "latent TB infection") 
data3 <- exprs(gse[, data3_ind])
meta3 <- meta_info[data3_ind, ]
group3 <- factor(meta3$`disease status:ch1`, 
                 levels=c("latent TB infection", "active tuberculosis"))
group3 <- plyr::revalue(group3, c("latent TB infection"=0, "active tuberculosis"=1))

rm(gse, meta_info)
rm(meta1, meta2, meta3)


##  Annotate new data
annotateNewData <- function(dat){
  probeIDs <- rownames(dat)
  x <- illuminaHumanv4SYMBOL
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  mapped_symbols <- sapply(probeIDs, function(id) xx[[id]])
  keep <- !is.na(mapped_symbols)
  mapped_symbols <- mapped_symbols[keep]
  dat <- dat[keep, ]
  rownames(dat) <- mapped_symbols
  return(dat)
}

data1 <- annotateNewData(data1)
data2 <- annotateNewData(data2)
data3 <- annotateNewData(data3)
# identical(rownames(data1), rownames(data2))


##  Previous data
#load("~/Documents/MSBE/newTB/new_combined_unmatched.RData")
data4 <- dat_lst[["US"]]; group4 <- label_lst[["US"]]
data5 <- dat_lst[["Africa"]]; group5 <- label_lst[["Africa"]]
data6 <- dat_lst[["India"]]; group6 <- label_lst[["India"]]

overlapped_genes <- intersect(rownames(data1), rownames(data4))
data1 <- data1[overlapped_genes, ]; data2 <- data2[overlapped_genes, ]; data3 <- data3[overlapped_genes, ]
data4 <- data4[overlapped_genes, ]; data5 <- data5[overlapped_genes, ]; data6 <- data6[overlapped_genes, ]

old_names <- names(dat_lst)
rm(dat_lst, label_lst)


##  Merge
dat_lst <- list(data1, data2, data3, data4, data5, data6)
label_lst <- list(group1, group2, group3, group4, group5, group6)
sapply(dat_lst, dim)
sapply(label_lst, table)
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", old_names)
names(dat_lst) <- study_names
names(label_lst) <- study_names


# ##  PCA
# dat_merged <- t(do.call(cbind, dat_lst))
# condition <- factor(do.call(c, lapply(label_lst, function(label) as.character(label))))
# batch <- factor(rep(1:length(dat_lst), sapply(dat_lst, ncol)))
# pca_res <- prcomp(dat_merged, center=TRUE, scale.=TRUE)
# pca_df <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], 
#                      Condition=condition, Batch=batch)
# ggplot(pca_df, aes(x=PC1, y=PC2, color=Batch, shape=Condition)) +
#   geom_point() +
#   labs(x=sprintf("PC1: %s percent variance", 
#                  scales::percent(pca_res$sdev[1]^2/sum(pca_res$sdev^2))),
#        y=sprintf("PC2: %s percent variance", 
#                  scales::percent(pca_res$sdev[2]^2/sum(pca_res$sdev^2))))
# 
# 
# 
# ##  Signal-to-noise analysis
# # biological signal
# s="GSE37250_SA"
# res_bio <- list()
# for(s in study_names){
#   ## Get training & test set
#   test_name <- s
#   train_name <- setdiff(study_names, test_name)
#   
#   # training set
#   dat <- do.call(cbind, dat_lst[train_name])
#   batch <- rep(1:(length(study_names)-1), times=sapply(dat_lst[train_name], ncol))
#   batch_names <- levels(factor(batch))
#   group <- do.call(c, label_lst[train_name])
#   
#   # test
#   dat_test <- dat_lst[[test_name]]
#   group_test <- label_lst[[test_name]]
#   
#   # biological signal in test set
#   design_bio <- model.matrix(~group_test)
#   fit_bio <- lmFit(dat_test, design_bio)
#   fit_bio <- eBayes(fit_bio)
#   res_bio[[test_name]] <- topTable(fit_bio, coef=2, adjust="BH", number=nrow(dat_test))
# }
# 
# lapply(res_bio, function(res){summary(exp(res$logFC))})
# 
# logFC_res <- lapply(res_bio, function(res){exp(res$logFC)})
# boxplot(logFC_res)
# 
# 
# ## batch level across negative samples in training set
# mean_batch <- var_batch <- matrix(NA, nrow=length(dat_lst), ncol=length(dat_lst),
#                                   dimnames=list(names(dat_lst), names(dat_lst)))
# for(i in seq_along(dat_lst)){
#   for(j in seq_along(dat_lst)){
#     if(i!=j){
#       dat1 <- dat_lst[[i]][, label_lst[[i]]==0]
#       dat2 <- dat_lst[[j]][, label_lst[[j]]==0]
#       
#       mean_batch[names(dat_lst)[i], names(dat_lst)[j]] <- abs(mean(dat1) - mean(dat2))
#       
#       rowvars1 <- mean(matrixStats::rowVars(dat1))
#       rowvars2 <- mean(matrixStats::rowVars(dat2))
#       var_batch[names(dat_lst)[i], names(dat_lst)[j]] <- 
#         ifelse(rowvars1 > rowvars2, rowvars1/rowvars2, rowvars2/rowvars1)
#     }
#   }
# }
# 
# round(mean_batch, 2)
# round(var_batch, 2)


##  Save data
label_lst <- lapply(label_lst, function(lb){as.numeric(as.character(lb))})
save(dat_lst, label_lst, file="./data/TB_real_data.RData")
