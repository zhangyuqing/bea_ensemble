rm(list=ls())
setwd("~/yuqingz/multistudy_batcheffect/diff_expr/")
sapply(c("ggplot2", "gridExtra", "reshape2"), require, character.only=TRUE)

res <- read.csv("simDE_avg.csv", header=TRUE)
res_mlt <- melt(res)
res_mlt <- data.frame(do.call(rbind, strsplit(as.character(res_mlt$variable),"_")), res_mlt)
colnames(res_mlt) <- c("Method", "Stat", "M_S", "Value")
  
png("boxplots_avg.png", width=7, height=4, units="in", res=300)
plt1 <- ggplot(res_mlt[res_mlt$Stat=="tpr", ], aes(x=Method, y=Value)) +
  geom_boxplot() +
  facet_grid(~Stat) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="")
plt2 <- ggplot(res_mlt[res_mlt$Stat=="fpr", ], aes(x=Method, y=Value)) +
  geom_boxplot() +
  facet_grid(~Stat) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="")
print(grid.arrange(plt1, plt2, ncol=2))
dev.off()  
