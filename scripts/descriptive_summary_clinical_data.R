#### DESCRIPTIVE TABLE OF CLINICAL DATA ####
#open data (patient clinical data)
data_clinical_patient <-
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_clinical_patient.txt", comment.char="#", header = T)

View(data_clinical_patient)
str(data_clinical_patient)


stats <- c("SEX", "AJCC_PATHOLOGIC_TUMOR_STAGE", 
           "AJCC_STAGING_EDITION", "HISTORY_NEOADJUVANT_TRTYN",
           "ICD_10", "ICD_O_3_HISTOLOGY", "ICD_O_3_SITE", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
           "PATH_M_STAGE", "PATH_N_STAGE", "PATH_T_STAGE", "PERSON_NEOPLASM_CANCER_STATUS", 
           "PRIOR_DX", "RACE", "RADIATION_THERAPY", "OS_STATUS", "DSS_STATUS",
           "DFS_STATUS", "PFS_STATUS", "GENETIC_ANCESTRY_LABEL")

table <-
  data_clinical_patient %>%
  select(all_of(stats)) %>%
  tbl_summary() %>%
  save.image()

