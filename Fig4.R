library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(tibble)
library(stringr)
library(cols4all)
library(ggrepel)
library(tidyverse)

## ---- fig4a ----
load("02-result/step00_promoter-new.RData")

count <- peak_promoter

count <- count %>% column_to_rownames("peak")

exprSet = count

dat = log2(edgeR::cpm(exprSet)+1)
dat[1:4,1:4] 
colnames(exprSet)
dat=t(dat)

data_pca <- prcomp(dat, scale. = TRUE)
pca_sum=summary(data_pca)

pca_sample <- as.data.frame(predict(data_pca))

pca_sample$Group <- c(rep(c("A","B","C"),each=8))

pca_sample$Sample <- paste0("S",c(rep(c(1:8),3)))
pca_sum   ##

(all <- ggplot(data = pca_sample, aes(x = PC1, y = PC2)) +
    theme_classic()+
    theme(axis.line = element_line(colour="black"),
          axis.text = element_text(hjust=0.5, vjust=0.5,size = 6,colour="black"),
          axis.title= element_text(size=8),
          legend.text= element_text(size=6,color = "black"),
          legend.title =  element_text(size=8,color = "black"),
          legend.key.size = unit(10, "pt"),
          plot.title = element_text(color="BLACK", size=10,hjust = 0.5)) +
    labs(x = 'PC1(13.53%)', y = 'PC2(6.56%)')+
    geom_point(aes(color =  Group,shape=Sample),size=1) + 
    scale_shape_manual(values = c(3,4,7,8,15:18))+
    scale_color_manual(values = c4a("dark2", 3)) 
)
ggsave("02-result/fig4a.pdf",width = 89,height = 59,units = "mm")

