########  Observe real data
rm(list=ls())
setwd("C:/Users/zhang/Dropbox/Work/MultiStudy_BatchEffect/differential_expression")
library(limma)
load("AEGIS_data.RData")

ff <- lmFit(train_expr, design=model.matrix(~y_train))
ff <- eBayes(ff)
tb <- topTable(ff, number=nrow(train_expr))
up_ind <- rownames(tb)[tb$logFC>0] 

round(mean(train_expr[up_ind, y_train==0]), 2)  # 8.32
round(mean(train_expr[up_ind, y_train==1]), 2)  # 8.55
round(mean(apply(train_expr[up_ind, y_train==0], 1, sd)), 2)  # 0.35
round(mean(apply(train_expr[up_ind, y_train==1], 1, sd)), 2)  # 0.37


########  Simulate data
rm(list=ls())
setwd("C:/Users/zhang/Dropbox/Work/MultiStudy_BatchEffect/differential_expression")
set.seed(1)

N_genes <- 2000
N_DE <- 100
ground_truth <- paste0("gene", 1:N_DE)
N_samples <- c(50, 50)
y_train <- c(rep(0, N_samples[1]), rep(1, N_samples[2]))
y_train <- factor(y_train)
  
train_expr <- matrix(rep(NA, N_genes*sum(N_samples)), nrow=N_genes, ncol=sum(N_samples),
                     dimnames=list(paste0("gene",1:N_genes), paste0("sample",1:sum(N_samples))))
for(i in 1:nrow(train_expr)){
  if(i<=N_DE){
    # DE genes
    train_expr[i, y_train==0] <- rnorm(N_samples[1], mean=8.32, sd=0.36) # control
    train_expr[i, y_train==1] <- rnorm(N_samples[2], mean=8.55, sd=0.36) # case
  }else{
    # NULL genes
    train_expr[i, ] <- rnorm(sum(N_samples), mean=8.32, sd=0.36)
  }
}

save(train_expr, y_train, N_genes, N_DE, N_samples, ground_truth, file="sim_data.RData")
