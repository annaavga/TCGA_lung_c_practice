#### ANALYSIS OF BOTH DATA FRAMES JOINED ####
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

data_clinical_patient <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/processed_data/data_clinical_patient.txt", 
             comment.char="#", header = T)
str(data_clinical_patient)
data_clinical_patient$AGE <- as.numeric(data_clinical_patient$AGE)
str(data_clinical_patient)


data_gene_expression <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/processed_data/data_gene_expression.txt", 
             header = T)

data_gene_expr <- data_gene_expression %>% t() %>% as.data.frame() %>% rownames_to_column(var = "PATIENT.ID")





