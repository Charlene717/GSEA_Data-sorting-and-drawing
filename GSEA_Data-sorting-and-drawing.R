
## Clear variables
rm(list = ls())


setwd(getwd()) ## Set current working directory
PathName <- getwd() ## Set output directroy


## Load data
FolderName <- c("#_RhoA_BRCA_20210510")

GSEA_Result_Ori <- read.csv(paste0(PathName,"/",FolderName,"_merge.csv"))
title <- ("TCGA BRCA RhoA")

## Sorting data

# Filtered by NES
GSEA_Result_1 <- GSEA_Result_Ori[abs(GSEA_Result_Ori$NES) >= 1,]

# Filtered by p-value
GSEA_Result_2 <- GSEA_Result_1[abs(GSEA_Result_1$NOM.p.val) <= 0.05,]

# Filtered by FDR
GSEA_Result_3 <- GSEA_Result_2[abs(GSEA_Result_2$FDR.q.val) <= 0.25,]
GSEA_Result_4 <- GSEA_Result_3[,c(2,6:8)]
colnames(GSEA_Result_4)[1] <- c("Pathway")

# GSEA_Result_5 <- GSEA_Result_4[order(GSEA_Result_4$NES,-GSEA_Result_4$NOM.p.val),]
GSEA_Result_5 <- GSEA_Result_4[order(GSEA_Result_4$NES),]
## Export the Sorting data

GSEA_Result_6 <- GSEA_Result_5[1:length(GSEA_Result_5[,1]),]
#GSEA_Result_6 <- GSEA_Result_5[1:30,]

write.table(GSEA_Result_5,file=paste0(PathName,"/",FolderName,"_GSEA_Sorting.csv"),quote = FALSE,row.names = FALSE,col.names = FALSE, na = "",sep = ',')


## Plot
library(ggplot2)

ggplot(GSEA_Result_6, aes(Pathway, NES, fill = FDR.q.val)) + geom_bar(stat = 'identity') +
       coord_flip()

library(dplyr)
library(forcats)

GSEA_Result_6 %>%
  mutate(Pathway = fct_reorder(Pathway, NES)) %>%
  ggplot(aes(x=Pathway, y=NES)) +
  geom_bar(stat="identity", fill="#a74cd9", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  ggtitle(title)

# GSEA_Result_5 %>%
#   mutate(Pathway = fct_reorder(Pathway, NES)) %>%
#   ggplot(aes(x=Pathway, y=NES, fill = NOM.p.val)) +
#   geom_bar(stat="identity", alpha=.6, width=.4) +
#   coord_flip() +
#   xlab("") +
#   theme_bw()

## Save the plot

pdf(paste0(PathName,"/",FolderName,"_GSEA_Sorting.pdf"))
GSEA_Result_6 %>%
  mutate(Pathway = fct_reorder(Pathway, NES)) %>%
  ggplot(aes(x=Pathway, y=NES)) +
  geom_bar(stat="identity", fill="#a74cd9", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  ggtitle(title)
dev.off()


png(paste0(PathName,"/",FolderName,"_GSEA_Sorting.png"))
GSEA_Result_6 %>%
  mutate(Pathway = fct_reorder(Pathway, NES)) %>%
  ggplot(aes(x=Pathway, y=NES)) +
  geom_bar(stat="identity", fill="#a74cd9", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  ggtitle(title)
dev.off()
