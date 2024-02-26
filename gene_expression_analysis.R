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
