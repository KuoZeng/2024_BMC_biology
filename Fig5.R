library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(tibble)
library(stringr)
library(cols4all)
library(ggrepel)
library(tidyverse)
library(randomForest)
library(rfPermute)
library(doParallel)

# ---- fig5a ----
load("02-result/step00_promoter-new.RData")
count <- peak_promoter

count <- count %>% column_to_rownames("peak")

count <- count[,!str_detect(colnames(count),"C")]
group <- data.frame("Sample"=colnames(count),
                    "Group"=c(rep(c("A","B"),each=8)))

train <- log2(edgeR::cpm(count)+1)
train <- as.data.frame(t(train))
group$Group <- as.factor(group$Group)

set.seed(20230815)
rf = randomForest(train, group$Group, importance=TRUE, proximity=TRUE, ntree = 500)
print(rf)

imp= as.data.frame(rf$importance)
imp = imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]
imp <- imp %>% rownames_to_column("peak")

name <- rio::import("02-result/04-DEG/04-DEG-AvsB-deseq2.txt")
colnames(name)[1] <- "peak"
imp1 <- left_join(imp,name[,c(1,35)],by="peak")
imp1 <- imp1[!duplicated(imp1[,1]),]

#write.table(imp1,"02-result/032-importance-AvsB.txt",sep = "\t",row.names = F,quote = F)


distance.matrix <- dist(1-rf$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

mds.values <- mds.stuff$points
mds.data <- data.frame(S=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Group= group$Group,
                       Sample=paste0("S",c(rep(c(1:8),2))))

mds_all <- ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(color=Group,shape=Sample),size=1)+
  theme_classic() +
  theme(axis.line = element_line(colour="black"),
        axis.text = element_text(hjust=0.5, vjust=0.5,size = 6,colour="black"),
        axis.title= element_text(size=8),
        legend.text= element_text(size=6,color = "black"),
        legend.title =  element_text(size=8,color = "black"),
        legend.key.size = unit(10, "pt"),
        plot.title = element_text(color="BLACK", size=10,hjust = 0.5))+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  stat_ellipse(aes(color = Group), geom = 'polygon', level = 0.68, alpha = 0, show.legend = FALSE)+
  scale_color_manual(values = c4a("dark2", 3)[1:2])+scale_shape_manual(values = c(3,4,7,8,15:18));mds_all

ggsave(mds_all,filename = "02-result/fig5a.pdf",width = 89,height = 59,units = "mm")


# ---- fig5b ----
rm(list= ls())
load("02-result/step00_promoter-new.RData")
count <- peak_promoter

count <- count %>% column_to_rownames("peak")

count <- count[,!str_detect(colnames(count),"B")]
group <- data.frame("Sample"=colnames(count),
                    "Group"=c(rep(c("A","C"),each=8)))

train <- log2(edgeR::cpm(count)+1)
train <- as.data.frame(t(train))
group$Group <- as.factor(group$Group)

set.seed(20230815)
rf = randomForest(train, group$Group, importance=TRUE, proximity=TRUE, ntree = 500)
print(rf)

imp= as.data.frame(rf$importance)
imp = imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]
imp <- imp %>% rownames_to_column("peak")

name <- rio::import("02-result/04-DEG/04-DEG-AvsC-deseq2.txt")
colnames(name)[1] <- "peak"
imp1 <- left_join(imp,name[,c(1,35)],by="peak")
imp1 <- imp1[!duplicated(imp1[,1]),]

#write.table(imp1,"02-result/032-importance-AvsC.txt",sep = "\t",row.names = F,quote = F)

distance.matrix <- dist(1-rf$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

mds.values <- mds.stuff$points
mds.data <- data.frame(S=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Group= group$Group,
                       Sample=paste0("S",c(rep(c(1:8),2))))

mds_all <- ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(color=Group,shape=Sample),size=1)+
  theme_classic() +
  theme(axis.line = element_line(colour="black"),
        axis.text = element_text(hjust=0.5, vjust=0.5,size = 6,colour="black"),
        axis.title= element_text(size=8),
        legend.text= element_text(size=6,color = "black"),
        legend.title =  element_text(size=8,color = "black"),
        legend.key.size = unit(10, "pt"),
        plot.title = element_text(color="BLACK", size=10,hjust = 0.5))+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  stat_ellipse(aes(color = Group), geom = 'polygon', level = 0.68, alpha = 0, show.legend = FALSE)+
  scale_color_manual(values = c4a("dark2", 3)[c(1,3)])+scale_shape_manual(values = c(3,4,7,8,15:18));mds_all


ggsave(mds_all,filename = "02-result/fig5b.pdf",width = 89,height = 59,units = "mm")

