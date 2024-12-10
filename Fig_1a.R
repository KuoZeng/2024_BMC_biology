library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)
library(tibble)
library(stringr)
library(cols4all)
library(ggrepel)

## PCA
count <- rio::import("02-result/01-count.txt")
count <- count %>% column_to_rownames("ID")

exprSet = count

dat = log2(edgeR::cpm(exprSet)+1)
dat[1:4,1:4] 
colnames(exprSet)
dat=t(dat)

data_pca <- prcomp(dat, scale. = TRUE)
pca_sum=summary(data_pca)

pca_sample <- as.data.frame(predict(data_pca))

pca_sample$Group <- c(rep(c("A","B","C"),each=8))
rownames(pca_sample) <- str_split(rownames(pca_sample),"_",simplify = T)[,1]
pca_sample$Sample <- paste0("S",c(rep(c(1:8),3)))

pca_sum  
# c4a_gui()
(all <- ggplot(data = pca_sample, aes(x = PC1, y = PC2)) +
    theme_classic()+
    theme(axis.line = element_line(colour="black"),
          axis.text = element_text(hjust=0.5, vjust=0.5,size = 6,colour="black"),
          axis.title= element_text(size=8),
          legend.text= element_text(size=6,color = "black"),
         
          legend.title =  element_text(size=8,color = "black"),
          legend.key.size = unit(10, "pt"),
          plot.title = element_text(color="BLACK", size=10,hjust = 0.5)) +
    labs(x = 'PC1(37.46%)', y = 'PC2(10.75%)')+
    geom_point(aes(color =  Group,shape=Sample),size=1) + 
    scale_shape_manual(values = c(3,4,7,8,15:18))+
   
    scale_color_manual(values = c4a("dark2", 3))
  
)

ggsave("02-result/Fig1-PCA.pdf",width = 89,height = 59,units = "mm")

