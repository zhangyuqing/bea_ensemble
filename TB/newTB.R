rm(list=ls())
sapply(c("GEOquery", "annotate", "hugene11sttranscriptcluster.db", 
         "SummarizedExperiment"), require, character.only=TRUE)
set.seed(123)

##  Download data from GEO
gse <- getGEO("GSE73408", GSEMatrix=TRUE)[[1]]
#exprs(gse)[1:5,1:5]


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
sample_info <- read.csv("~/Dropbox/Work/MultiStudy_BatchEffect/TB/new_data_info.csv", 
                        as.is=TRUE)
identical(sample_info$ID, colnames(dat))
colnames(dat) <- sample_info$Label

dat <- dat[, -grep("_PNA_", colnames(dat))]
group <- rep(0, ncol(dat))
group[grep("_TB_", colnames(dat))] <- 1


##  Clean together with the other studies
rds_obj <- readRDS(file.path("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB", 
                             "combined.rds"))
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
# #rm_id <- africa_id[1:round(length(africa_id)-table(rds_obj$Label[africa_id])['Active']/0.5)]  #0.588)]
# rds_obj <- rds_obj[, -rm_id]
# new_africa_group <- rds_obj$Label[rds_obj$SequencingBatch=="Africa"]
# print(sprintf("Prevalence in Africa: %s", 
#               round(table(new_africa_group)['Active']/sum(table(new_africa_group)),2)))

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

save(dat_lst, label_lst, file="~/Documents/MSBE/newTB/new_combined_unmatched.RData")

