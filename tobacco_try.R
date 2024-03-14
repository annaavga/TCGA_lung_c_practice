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
median_expression <- expression_tcga %>% 
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
TLS_groups <- mutate(mean_expression = cut(mean_expression$., 
                        breaks = c(-Inf, median_expression, Inf), 
                        labels = c("LOW", "HIGH"),
                        include.lowest = F))

