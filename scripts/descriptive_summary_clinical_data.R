#### DESCRIPTIVE TABLE OF CLINICAL DATA ####
#open data (patient clinical data)
data_clinical_patient <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_clinical_patient.txt", comment.char="#", header = T)

set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

View(data_clinical_patient)
str(data_clinical_patient)

data_clinical_patient$AGE <- as.numeric(data_clinical_patient$AGE)

stats <- c("SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE", 
           "AJCC_STAGING_EDITION", "HISTORY_NEOADJUVANT_TRTYN",
           "ICD_10", "ICD_O_3_HISTOLOGY", "ICD_O_3_SITE", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
           "PATH_M_STAGE", "PATH_N_STAGE", "PATH_T_STAGE", "PERSON_NEOPLASM_CANCER_STATUS", 
           "PRIOR_DX", "RACE", "RADIATION_THERAPY", "OS_STATUS", "DSS_STATUS",
           "DFS_STATUS", "PFS_STATUS", "GENETIC_ANCESTRY_LABEL")

clinical_data_descriptive_tab <-
  data_clinical_patient %>%
  select(all_of(stats)) %>%
  tbl_summary()
clinical_data_descriptive_tab

clinical_data_descriptive_tab %>%
  gtsummary::as_tibble() %>% 
  writexl::write_xlsx("./results/clin_desc_table.xlsx")

write.table(data_clinical_patient,"data_clinical_patient.txt",sep="\t")
