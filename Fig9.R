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

# ---- fig9a ----
rm(list= ls())
path <- "PPI/MCC-top10-enrichment/MCC-TOP10-GO/"

bp <- rio::import(paste0(path,"AvsB-BP-10.txt") )
cc <- rio::import(paste0(path,"AvsB-CC-10.txt") )
mf <- rio::import(paste0(path,"AvsB-MF-10.txt") )

table(colnames(bp) == colnames(mf))

dat <- rbind(bp,cc)
dat <- rbind(dat,mf)

dat$`-log10(p)` <- -log10(dat$pvalue)
dat$Description <- factor(dat$Description,levels = dat$Description,ordered = T)

color <- c4a("set2", 3)

dat <- dat[complete.cases(dat[, c(1)]), ]

ggplot(dat) +
  aes(x = Description, y = Count, fill = ONTOLOGY) +
  geom_bar(stat = "identity",colour="black",width = 0.7) +
  scale_fill_manual(values =color)+
  theme(
    axis.title=element_text(size=10,color="black"),
    axis.text = element_text(size=8,color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.direction = "horizontal",
    legend.position = c(0.8,0.9),
    panel.background = element_rect(fill = "transparent",colour = "black"),
    plot.margin=unit(c(1,1,1,5),'lines')
  )
ggsave(paste0(path,"fig9a.pdf"),
       width = 160,height = 110,units = "mm")


# ---- fig9a ----
rm(list= ls())
path <- "PPI/MCC-top10-enrichment/MCC-TOP10-GO/"

bp <- rio::import(paste0(path,"AvsB-BP-10.txt") )
cc <- rio::import(paste0(path,"AvsB-CC-10.txt") )
mf <- rio::import(paste0(path,"AvsB-MF-10.txt") )

table(colnames(bp) == colnames(mf))

dat <- rbind(bp,cc)
dat <- rbind(dat,mf)

dat$`-log10(p)` <- -log10(dat$pvalue)
dat$Description <- factor(dat$Description,levels = dat$Description,ordered = T)

color <- c4a("set2", 3)

dat <- dat[complete.cases(dat[, c(1)]), ]

ggplot(dat) +
  aes(x = Description, y = Count, fill = ONTOLOGY) +
  geom_bar(stat = "identity",colour="black",width = 0.7) +
  scale_fill_manual(values =color)+
  theme(
    axis.title=element_text(size=10,color="black"),
    axis.text = element_text(size=8,color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.direction = "horizontal",
    legend.position = c(0.8,0.9),
    panel.background = element_rect(fill = "transparent",colour = "black"),
    plot.margin=unit(c(1,1,1,5),'lines')
  )
ggsave(paste0(path,"fig9a.pdf"),
       width = 160,height = 110,units = "mm")


# ---- fig9b ----
rm(list= ls())
path <- "PPI/MCC-top10-enrichment/MCC-TOP10-GO/"

bp <- rio::import(paste0(path,"AvsC-BP-10.txt") )
cc <- rio::import(paste0(path,"AvsC-CC-10.txt") )
mf <- rio::import(paste0(path,"AvsC-MF-10.txt") )

table(colnames(bp) == colnames(mf))

dat <- rbind(bp,cc)
dat <- rbind(dat,mf)

dat$`-log10(p)` <- -log10(dat$pvalue)
dat$Description <- factor(dat$Description,levels = dat$Description,ordered = T)

color <- c4a("set2", 3)

dat <- dat[complete.cases(dat[, c(1)]), ]

ggplot(dat) +
  aes(x = Description, y = Count, fill = ONTOLOGY) +
  geom_bar(stat = "identity",colour="black",width = 0.7) +
  scale_fill_manual(values =color)+
  theme(
    axis.title=element_text(size=10,color="black"),
    axis.text = element_text(size=8,color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.direction = "horizontal",
    legend.position = c(0.8,0.9),
    panel.background = element_rect(fill = "transparent",colour = "black"),
    plot.margin=unit(c(1,1,1,5),'lines')
  )
ggsave(paste0(path,"fig9b.pdf"),
       width = 160,height = 110,units = "mm")