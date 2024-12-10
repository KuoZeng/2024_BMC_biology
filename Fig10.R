library(clusterProfiler)
library(ggthemes)
# library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(tibble)
library(cols4all)
library(tidyr)

# ---- fig10a ----
rm(list= ls())
path <- "PPI/MCC-top10-enrichment/"
kegg <- rio::import(paste0(path,"MCC-top10-AvsB-KEGG.txt"))
kegg_top5 <- kegg[1:5,]

kegg_top5$Description <- factor(kegg_top5$Description,levels = kegg_top5$Description,ordered = T)
colnames(kegg_top5)[11] <- "-log10.P.value"

ggplot(kegg_top5, aes(x = `-log10.P.value`, y = rev(Description), fill = `-log10.P.value`))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x=0.1,y=rev(Description),label = Description),size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
  )+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#65B0C6",high ="#CB5640")+
  geom_text(data = kegg_top5,
            aes(x = 0.1, y = rev(Description), label = geneID,color=`-log10.P.value`),
            size = 3,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_gradient(low = "#65B0C6",high ="#CB5640")+
  labs(title = "Enrichment of genes",
       y="Description")

ggsave(paste0(path,"fig10a.pdf") ,
       width = 110,height = 70,units = "mm")

# ---- fig10b ----
rm(list= ls())
path <- "PPI/MCC-top10-enrichment/"
kegg <- rio::import(paste0(path,"MCC-top10-AvsC-KEGG.txt"))
kegg_top10 <- kegg[1:10,]

kegg_top10$Description <- factor(kegg_top10$Description,levels = kegg_top10$Description,ordered = T)
colnames(kegg_top10)[11] <- "-log10.P.value"

ggplot(kegg_top10, aes(x = `-log10.P.value`, y = rev(Description), fill = `-log10.P.value`))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x=0.1,y=rev(Description),label = Description),size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
  )+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#65B0C6",high ="#CB5640")+
  geom_text(data = kegg_top10,
            aes(x = 0.1, y = rev(Description), label = geneID,color=`-log10.P.value`),
            size = 3,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_gradient(low = "#65B0C6",high ="#CB5640")+
  labs(title = "Enrichment of genes",
       y="Description")

ggsave(paste0(path,"fig10b.pdf") ,
       width = 110,height = 70,units = "mm")

