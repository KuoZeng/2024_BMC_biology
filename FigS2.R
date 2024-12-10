library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(tibble)
library(stringr)
library(cols4all)
library(ggrepel)
library(tidyverse)
library(pheatmap)

# ---- figS2a ----

datAB <- rio::import("02-result/04-DEG/04-DEG-AvsB-deseq2.txt")
dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
dat <- dat %>% column_to_rownames("Sample")

sub_abu = dat[,datAB$Row.names]
sub_abu <- as.data.frame(t(sub_abu))

sub_abu <- sub_abu %>% rownames_to_column("ID")
colnames(datAB)[1] <- "ID"
sub_abu <- merge(sub_abu,datAB,by="ID")
sub_abu <- sub_abu %>% column_to_rownames("gene_name")
sub_abu <- sub_abu[,-1]

bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

sub_abu2 <- sub_abu[,c(1:16)]
col_anno2 <- data.frame(Group=c(rep(c("A","B"),each=8)),
                       row.names=colnames(sub_abu2))
col_anno2$Group <- as.factor(col_anno2$Group)
ann_colors2 <-  list(Group = c("A"=c4a("dark2", 3)[1],"B"=c4a("dark2", 3)[2])) 


pheatmap(sub_abu2,scale='row' ,cluster_cols=F, fontsize=3,show_rownames=T,
         show_colnames=F,
         filename="02-result/04-DEG/figS2a.pdf",
         cellheight = 3.5, cellwidth = 8, annotation_col=col_anno2, border_color=NA, 
         annotation_colors =ann_colors2,
         color = c(colorRampPalette(colors = c("#2fa1dd","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f87669"))(length(bk)/2)),
         legend_breaks=seq(-4,4,2),breaks=bk,treeheight_row=15,
         gaps_col = c( 8))


# ---- figS2b ----
rm(list= ls())
datAC <- rio::import("02-result/04-DEG/04-DEG-AvsC-deseq2.txt")
dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
dat <- dat %>% column_to_rownames("Sample")

sub_abu = dat[,datAC$Row.names]
sub_abu <- as.data.frame(t(sub_abu))

sub_abu <- sub_abu %>% rownames_to_column("ID")
colnames(datAC)[1] <- "ID"
sub_abu <- merge(sub_abu,datAC,by="ID")
sub_abu <- sub_abu %>% column_to_rownames("gene_name")
sub_abu <- sub_abu[,-1]

bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

sub_abu2 <- sub_abu[,c(1:16)]
col_anno2 <- data.frame(Group=c(rep(c("A","C"),each=8)),
                       row.names=colnames(sub_abu2))
col_anno2$Group <- as.factor(col_anno2$Group)
ann_colors2 <-  list(Group = c("A"=c4a("dark2", 3)[1],"C"=c4a("dark2", 3)[3])) 


pheatmap(sub_abu2,scale='row' ,cluster_cols=F, fontsize=3,show_rownames=T,
         show_colnames=F,
         filename="02-result/04-DEG/figS2b.pdf",
         cellheight = 5, cellwidth = 8, annotation_col=col_anno2, border_color=NA, 
         annotation_colors =ann_colors2,
         color = c(colorRampPalette(colors = c("#2fa1dd","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f87669"))(length(bk)/2)),
         legend_breaks=seq(-4,4,2),breaks=bk,treeheight_row=15,
         gaps_col = c( 8))


# ---- figS2c ----
rm(list= ls())
dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
dat <- dat[ !dat$Group %in% "C",] 
meta <- dat[,1:2]

dat <- dat %>% column_to_rownames("Sample")
meta <- meta %>% column_to_rownames("Sample")

train <- dat[,-1]
meta$Group <- as.factor(meta$Group)

set.seed(123)
boruta <- Boruta(x=train, y=meta$Group, pValue=0.01, mcAdj=T, 
                 maxRuns=100)
boruta
boruta$timeTaken
Boruta::plotImpHistory(boruta)

boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}
boruta.variable.imp <- boruta.imp(boruta)
head(boruta.variable.imp)

boruta.finalVars <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta")


library(pheatmap)
sub_abu = dat[,boruta.finalVars$Item]
sub_abu <- as.data.frame(t(sub_abu))

infor <- rio::import("02-result/01-information.txt")
sub_abu <- sub_abu %>% rownames_to_column("ID")

sub_abu <- merge(sub_abu,infor[,1:2],by="ID")
biomarker <- sub_abu[,c("ID","gene_name")]

sub_abu <- sub_abu %>% column_to_rownames("gene_name")
sub_abu <- sub_abu[,-1]

bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

col_anno <- data.frame(Group=c(rep(c("A","B"),each=8)),
                       row.names=colnames(sub_abu))
col_anno$Group <- as.factor(col_anno$Group)
ann_colors <-  list(Group = c("A"=c4a("dark2", 3)[1],"B"= c4a("dark2", 3)[2])) 

pheatmap(sub_abu,scale='row' ,cluster_cols=F, fontsize=8,show_rownames=T,
         show_colnames=T,
         filename="02-result/03-boruta/figS2c.pdf",
         cellheight = 15, cellwidth = 12, annotation_col=col_anno, border_color=NA,  annotation_colors =ann_colors,
         color = c(colorRampPalette(colors = c("#2fa1dd","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f87669"))(length(bk)/2)),
         legend_breaks=seq(-4,4,2),breaks=bk,treeheight_row=15,
         gaps_col = c( 8))


# ---- figS2d ----
rm(list= ls())
dat <- rio::import("02-result/021-PLS-DA-fpkm.txt")
dat <- dat[ !dat$Group %in% "B",]  
rownames(dat) <- dat$Sample
dat <- dat[,-1]

meta <- data.frame("Sample"=rownames(dat),"Group"=c(rep(c("A","C"),each=8)))

train <- dat[,-1]
meta$Group <- as.factor(meta$Group)

set.seed(123)
boruta <- Boruta(x=train, y=meta$Group, pValue=0.01, mcAdj=T, 
                 maxRuns=100)
boruta
boruta$timeTaken

table(boruta$finalDecision)

boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}
boruta.variable.imp <- boruta.imp(boruta)
head(boruta.variable.imp)

boruta.finalVars <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta")

library(pheatmap)
sub_abu = dat[,boruta.finalVars$Item]
sub_abu <- as.data.frame(t(sub_abu))

infor <- rio::import("02-result/01-information.txt")
sub_abu <- sub_abu %>% rownames_to_column("ID")

sub_abu <- merge(sub_abu,infor[,1:2],by="ID")

biomarker <- sub_abu[,c("ID","gene_name")]
write.table(biomarker,"02-result/03-boruta/importance-AvsC.txt",sep = "\t",row.names = F,quote = F)

sub_abu <- sub_abu %>% column_to_rownames("gene_name")
sub_abu <- sub_abu[,-1]

bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

col_anno <- data.frame(Group=c(rep(c("A","C"),each=8)),
                       row.names=colnames(sub_abu))
col_anno$Group <- as.factor(col_anno$Group)
ann_colors <-  list(Group = c("A"=c4a("dark2", 3)[1],"C"= c4a("dark2", 3)[3])) 

pheatmap(sub_abu,scale='row' ,cluster_cols=F, fontsize=8,show_rownames=T,
         show_colnames=T,
         filename="02-result/03-boruta/figS2d.pdf",
         cellheight = 15, cellwidth = 12, annotation_col=col_anno, border_color=NA,  annotation_colors =ann_colors,
         color = c(colorRampPalette(colors = c("#2fa1dd","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f87669"))(length(bk)/2)),
         legend_breaks=seq(-4,4,2),breaks=bk,treeheight_row=15,
         gaps_col = c( 8))

