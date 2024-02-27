#### ANALYSIS OF BOTH DATA FRAMES JOINED ####
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

##CLINICAL DATA
data_clinical_patient <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/processed_data/data_clinical_patient.txt", 
             comment.char="#", header = T)

#putting age as numeric
str(data_clinical_patient)
data_clinical_patient$AGE <- as.numeric(data_clinical_patient$AGE)
str(data_clinical_patient)

#putting the same names to the patients
data_clinical <- data_clinical_patient %>% 
  as.data.frame(str_replace_all(data_clinical_patient$PATIENT_ID,"-",".")) %>%
  rownames_to_column(var = "PATIENT.ID")

class(data_clinical$PATIENT.ID)

data_clinical$PATIENT.ID <- paste0(data_clinical$PATIENT.ID, ".01")

#GENE EXPRESSION DATA
data_gene_expression <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/processed_data/data_gene_expression.txt", 
             header = T)

data_gene_expr <- data_gene_expression %>% t() %>% as.data.frame() %>% rownames_to_column(var = "PATIENT.ID")



#Joining the dataframes of gene expression and clinical data
merged_data <- data_gene_expr %>% left_join(data_clinical, by= 'PATIENT.ID')
str(merged_data)

#transform the overall survival into 0 and 1
merged_data <- merged_data %>% mutate(OS_STATUS =  recode(OS_STATUS, "0:LIVING" = "0"),  OS_STATUS = recode(OS_STATUS, "1:DECEASED" = "1"))

#separate cohort at median


TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

TLS_HIGH <- data_gene_expression %>%
  as_tibble() %>%
  group_by(colMeans(data_gene_expression)> med)


TLS_HIGH <- merged_data %>% 
  as_tibble() %>%
  colMeans(colnames = TLS_signature) %>% as.data.frame() %>%

? group_by
