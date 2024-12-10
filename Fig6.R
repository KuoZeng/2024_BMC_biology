library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(tibble)
library(stringr)
library(cols4all)
library(ggrepel)
library(tidyverse)

# ---- fig6a ----

load("02-result/step00_promoter-new.RData")
count <- peak_promoter
count <- count %>% column_to_rownames("peak")

count <- count[,!str_detect(colnames(count),"C")]
group <- data.frame("Sample"=colnames(count),
                    "Group"=c(rep(c("A","B"),each=8)))
exprSet=count 
group_list <- group$Group

source('03-code/run_DEG_RNA-seq.R')

table(group_list)  
exprSet[1:4,1:4] 
dim(exprSet)

table(group_list)
run_DEG_RNAseq(exprSet,group_list,
               g1="A",g2="B",
               pro='02-result/04-DEG/04-DEG-AvsB')

# ---- fig6b ----
rm(list = ls())
load("02-result/step00_promoter-new.RData")
count <- peak_promoter
count <- count %>% column_to_rownames("peak")

count <- count[,!str_detect(colnames(count),"B")]
group <- data.frame("Sample"=colnames(count),
                    "Group"=c(rep(c("A","C"),each=8)))

exprSet=count 
group_list <- group$Group

source('03-code/run_DEG_RNA-seq.R')

table(group_list)  
exprSet[1:4,1:4] 
dim(exprSet)

table(group_list)
run_DEG_RNAseq(exprSet,group_list,
               g1="A",g2="C",
               pro='02-result/04-DEG/04-DEG-AvsC')

# ---- fig6c ----
rm(list = ls())
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
library(Boruta)

load("02-result/step00_promoter-new.RData")
count <- peak_promoter

count <- count %>% column_to_rownames("peak")

count <- count[,!str_detect(colnames(count),"C")]
group <- data.frame("Sample"=colnames(count),
                    "Group"=c(rep(c("A","B"),each=8)))

train <- log2(edgeR::cpm(count)+1)
train <- as.data.frame(t(train))
group$Group <- as.factor(group$Group)

set.seed(123)
boruta <- Boruta(x=train, y=group$Group, pValue=0.01, mcAdj=T, 
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


library(ImageGP)

(plot_ab <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90,manual_color_vector = rev(c4a("accent", 4)) ))

ggsave(filename = "02-result/03-RF-boruta/fig6c.pdf",width = 89,height = 59,
       units = "mm")

# ---- fig6d ----
rm(list = ls())
load("02-result/step00_promoter-new.RData")
count <- peak_promoter

count <- count %>% column_to_rownames("peak")

count <- count[,!str_detect(colnames(count),"B")]
group <- data.frame("Sample"=colnames(count),
                    "Group"=c(rep(c("A","C"),each=8)))

train <- log2(edgeR::cpm(count)+1)
train <- as.data.frame(t(train))
group$Group <- as.factor(group$Group)

set.seed(123)
boruta <- Boruta(x=train, y=group$Group, pValue=0.01, mcAdj=T, 
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

library(ImageGP)

(plot_ac <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90,manual_color_vector = rev(c4a("accent", 4)) ))

ggsave(filename = "02-result/03-RF-boruta/fig6d.pdf",width = 89,height = 59,
       units = "mm")