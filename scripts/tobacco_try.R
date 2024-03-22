####tobacco
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

data_mrna_seq_v2_rsem <- 
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", header = T)

#take columns you need
TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

mrna_expression <- subset(data_mrna_seq_v2_rsem, Hugo_Symbol %in% TLS_signature)
View(mrna_expression)

#put correct names and eliminate the columns which you don't need, transform
rownames(mrna_expression) <- mrna_expression[,1]
expression_tcga <- 
  mrna_expression[, -which(colnames(mrna_expression) == "Hugo_Symbol")] 
expression_tcga <-
  expression_tcga[, -which(colnames(expression_tcga)== "Entrez_Gene_Id")]
expression_tcga <- log2(expression_tcga + 1)

#median
mean_expression <- expression_tcga %>% 
  as_tibble() %>%
  colMeans() %>%
  as.data.frame()

median_expression <- expression_tcga %>% 
  as_tibble() %>%
  colMeans() %>% 
  median()

# separate intp tls high and low
groups_TLS <- cut(mean_expression$.,
                     breaks= c(-Inf, median_expression, Inf),
                     labels=c('TLS_LOW', 'TLS_HIGH'))
table(groups_TLS)

#join w clinical data
clinical_info <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/processed_data/data_clinical_patient.txt", 
             comment.char="#", header = T)

#putting age as numeric
str(clinical_info)
clinical_info$AGE <- as.numeric(clinical_info$AGE)
str(clinical_info)

#putting the same names to the patients
data_clinical <- clinical_info %>% 
  as.data.frame(str_replace_all(clinical_info$PATIENT_ID,"-",".")) %>%
  rownames_to_column(var = "PATIENT.ID")

class(data_clinical$PATIENT.ID)

data_clinical$PATIENT.ID <- paste0(data_clinical$PATIENT.ID, ".01")


gene_expr <- mean_expression %>% as.data.frame() %>% rownames_to_column(var = "PATIENT.ID")



#Joining the dataframes of gene expression and clinical data
merged <- gene_expr %>% left_join(data_clinical, by= 'PATIENT.ID')
str(merged)

by(merged, merged$., merged$)