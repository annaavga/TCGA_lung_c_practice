#### GSVA METHOD ####
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

# load data
data_mrna <- read.delim(
  "C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", 
  header = T)
#choose data from gene set
TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

data_gene_set <- subset(data_mrna, Hugo_Symbol %in% TLS_signature)
rownames(data_gene_set) <- data_gene_set[,1]
data_gene_set <- 
  data_gene_set[, -which(colnames(data_gene_set) == "Hugo_Symbol")] 
data_gene_set <-
  data_gene_set[, -which(colnames(data_gene_set)== "Entrez_Gene_Id")]
data_gene_set <- log2(data_gene_set + 1)

#take sample of genes, without those included in the dataset
data_no_set <- data_mrna[!(data_mrna$Hugo_Symbol %in% TLS_signature),]
data_sampled <- sample_n(data_no_set, 100)

rownames(data_sampled) <- data_sampled[,1]
data_sampled <- 
  data_sampled[, -which(colnames(data_sampled) == "Hugo_Symbol")] 
data_sampled <-
  data_sampled[, -which(colnames(data_sampled)== "Entrez_Gene_Id")]
data_sampled <- log2(data_sampled + 1)

#join the two dataframes
finaldat_GSVA <- data_sampled %>% bind_rows(data_gene_set)
str(finaldat_GSVA)

#prepare the df to do the analysis
finaldat_GSVA <- as.matrix(finaldat_GSVA)
class(finaldat_GSVA) <- "numeric"

#choose the gene signature
TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

TLS_signature <- as.list(TLS_signature)
class(TLS_signature)

#PERFORM GSVA
gsvaPar <- gsvaParam(finaldat_GSVA , TLS_signature)
gsvaPar

gsva.es <- gsva(gsvaPar, verbose=F)
#dim(gsva.es)

gsva.es <- as.data.frame(gsva.es)
class(gsva.es)

#add the names of the genes into the dataframe
names <- as.data.frame(data_mrna$Hugo_Symbol)
names <- subset(names, data_mrna$Hugo_Symbol %in% TLS_signature)

data_gsva <- cbind(gsva.es, names$`data_mrna$Hugo_Symbol`)

data_gsva <- data_gsva %>% 
  rename("GENE_NAME" = "names$`data_mrna$Hugo_Symbol`")

