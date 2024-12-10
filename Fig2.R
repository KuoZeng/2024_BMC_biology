
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


# ---- fig2a ----

dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
meta <- dat[,1:2]

dat <- dat %>% column_to_rownames("Sample")
meta <- meta %>% column_to_rownames("Sample")

train <- dat[,-1]
meta$Group <- as.factor(meta$Group)

set.seed(20230720)
rf = randomForest(train, meta$Group, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
imp <- imp %>% rownames_to_column("Gene_name")
imp <- arrange(imp, desc(MeanDecreaseAccuracy))

#write.table(imp,"02-result/031-importance.txt",sep = "\t",row.names = F,quote = F)


distance.matrix <- dist(1-rf$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.values <- mds.stuff$points

mds.data <- data.frame(S=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Group= meta$Group,
                       Sample=paste0("S",c(rep(c(1:8),3))))

mds.data$Group <- as.factor(mds.data$Group)

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
  scale_color_manual(values = c4a("dark2", 3))+scale_shape_manual(values = c(3,4,7,8,15:18));mds_all

ggsave(mds_all,filename = "02-result/Fig2a.pdf",width = 89,height = 59,units = "mm")



# ---- fig2b ----
rm(list=ls())
dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
dat <- dat[ !dat$Group %in% "C",]  
meta <- dat[,1:2]

dat <- dat %>% column_to_rownames("Sample")
meta <- meta %>% column_to_rownames("Sample")

train <- dat[,-1]
meta$Group <- as.factor(meta$Group)

set.seed(20230720)
rf = randomForest(train, meta$Group, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

imp= as.data.frame(rf$importance)
imp = imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]
imp <- imp %>% rownames_to_column("ID")

#write.table(imp,"02-result/031-importance-AvsB.txt",sep = "\t",row.names = F,quote = F)

distance.matrix <- dist(1-rf$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)

mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

mds.values <- mds.stuff$points

mds.data <- data.frame(S=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Group= meta$Group,
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

ggsave(mds_all,filename = "02-result/Fig2b.pdf",width = 89,height = 59,units = "mm")


# ---- fig2c ----
rm(list=ls())

dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
dat <- dat[ !dat$Group %in% "B",]  
rownames(dat) <- dat$Sample
dat <- dat[,-1]

meta <- data.frame("Sample"=rownames(dat),"Group"=c(rep(c("A","C"),each=8)))

train <- dat[,-1]
meta$Group <- as.factor(meta$Group)

set.seed(20230720)
rf = randomForest(train, meta$Group, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)

imp= as.data.frame(rf$importance)
imp = imp[order(imp$MeanDecreaseAccuracy,decreasing = T),]
imp <- imp %>% rownames_to_column("ID")

#write.table(imp,"02-result/031-importance-AvsC.txt",sep = "\t",row.names = F,quote = F)


distance.matrix <- dist(1-rf$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

mds.values <- mds.stuff$points
mds.data <- data.frame(S=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2],
                       Group= meta$Group,
                       Sample=paste0("S",c(rep(c(1:8),2))))

mds_all <- ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_point(aes(color=Group,shape=Sample),size=1)+
  theme_classic() +
  theme(axis.line = element_line(colour="black"),
        axis.text = element_text(hjust=0.5, vjust=0.5,size = 6,colour="black"),
        axis.title= element_text(size=8),
        legend.text= element_text(size=6,color = "black"),
        # legend.key = element_rect(fill = 'transparent'), 
        legend.title =  element_text(size=8,color = "black"),
        legend.key.size = unit(10, "pt"),
        plot.title = element_text(color="BLACK", size=10,hjust = 0.5))+
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  stat_ellipse(aes(color = Group), geom = 'polygon', level = 0.68, alpha = 0, show.legend = FALSE)+
  scale_color_manual(values = c4a("dark2", 3)[c(1,3)])+scale_shape_manual(values = c(3,4,7,8,15:18));mds_all

ggsave(mds_all,filename = "02-result/Fig2c.pdf",width = 89,height = 59,units = "mm")
