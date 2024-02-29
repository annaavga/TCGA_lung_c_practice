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





gene_names <- read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", 
)
gene_names <- select(gene_names$Hugo_Symbol)
name <- gene_names %>% 
  select_if(is.character) 
class(name) <- "character"

TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

TLS_signature <- as.list(TLS_signature)
class(TLS_signature)
TLS_NAMES <- as.character(TLS_signature)

gsvaPar <- gsvaParam(data_mrna, TLS_signature)
gsvaPar

gsva.es <- gsva(gsvaPar, verbose=F, annotation = name)
dim(gsva.es)


?gsva


dup <- data.frame(as.numeric(duplicated(data$X))) #creates df with binary var for duplicated rows 
colnames(dup) <- c("dup") #renames column for simplicity 
df2 <- cbind(data, dup) #bind to original df 
df3 <- subset(df2, dup == 1) #subsets df using binary var for duplicated`