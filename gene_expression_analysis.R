#### OPEN FILE AND PREPARE IT ####

data_mrna_seq_v2_rsem <- 
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/Raw_data/luad_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")
View(data_mrna_seq_v2_rsem)

###libraries
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("dslabs")
library (dslabs)

#take columns you need
TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

data <- subset(data_mrna_seq_v2_rsem, Hugo_Symbol %in% TLS_signature)
View(data)

#put correct names and eliminate the columns which you don't need, transform
rownames(data) <- data[,1]
dat <- 
  data[, -which(colnames(data) == "Hugo_Symbol")] 
dat <-
  dat[, -which(colnames(dat)== "Entrez_Gene_Id")]
dat <- log2(dat + 1)

