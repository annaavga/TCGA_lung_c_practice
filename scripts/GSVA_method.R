#### GSVA METHOD ####
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

# load data

data_mrna <- as.matrix(read.delim(
  "C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", 
  header = T))

rownames(data_mrna) <- data_mrna[,1]
data_mrna <- 
  data_mrna[, -which(colnames(data_mrna) == "Hugo_Symbol")] 
data_mrna <-
  data_mrna[, -which(colnames(data_mrna)== "Entrez_Gene_Id")]
#data_mrna <- log2(data_mrna + 1)


class(data_mrna) <- "numeric"
data_mrna <- data_mrna[complete.cases(data_mrna), ]

?complete.cases