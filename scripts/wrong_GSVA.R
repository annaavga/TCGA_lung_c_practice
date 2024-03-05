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

rownames(data_gsva) <- data_gsva$GENE_NAME
data_gsva <- 
  data_gsva[, -which(colnames(data_gsva) == "GENE_NAME")] 

#### DENSITY PLOT ####
density_plot_gsva <- data_gsva %>% 
  as_tibble() %>%
  colMeans() %>% as.data.frame() %>%
  # # dplyr::rename("mean_expression" = ".") %>%
  ggplot2::ggplot(aes(x = .)) +
  geom_density(fill = "lightblue2", alpha = 0.3) + 
  xlab("GSVA score") +
  ylab("Density") +
  my_theme +
  theme(legend.position = "none")
density_plot_gsva
ggsave("GSVA score.pdf")


####COMPARISON WITH MEDIAN ####
#gsva data
GSVA_COL <- gsva.es %>%
  colMeans() %>%
  data.frame() %>%
  rename("GSVA_MEAN" = '.') %>%
  mutate(GSVA_CAT =cut(GSVA_COL$GSVA_MEAN, 
      breaks = c(-Inf, 0, Inf), 
      labels = c("LOW", "HIGH"),
      include.lowest = T))
  
str(GSVA_COL)
class(GSVA_COL)

#median data
MEDIAN <- median(apply(data_gene_set, 2, median))

class(MEDIAN)

MEDIAN_COL <- data_gene_set %>%
  colMeans() %>%
  data.frame() %>%
  rename("MEDIAN_MEAN" = '.') %>%
  mutate(MEDIAN_CAT = cut(MEDIAN_COL$MEDIAN_MEAN, 
                       breaks = c(-Inf, MEDIAN, Inf), 
                       labels = c("LOW", "HIGH"),
                       include.lowest = F))


#merge dataframes and adding column of coincidence
GSVA_and_MEDIAN <- GSVA_COL %>% 
  cbind(MEDIAN_COL) %>%
  mutate(COINCIDENCE =  ifelse(GSVA_and_MEDIAN$GSVA_CAT == GSVA_and_MEDIAN$MEDIAN_CAT, "YES", "NO")) 

GSVA_and_MEDIAN %>%
  rownames_to_column() %>%
  write_xlsx("./results/GSVA_MEDIAN_comparison_table.xlsx")

#find n of coincidence
COINCIDENCE_GSVA_MED <- GSVA_and_MEDIAN %>%
  group_by(COINCIDENCE)%>%
  dplyr::summarise(n()) 
  
#PLOT of GSVA and MEDIAN
GSVA_and_MEDIAN %>%
  group_by(COINCIDENCE) %>%
  ggplot(aes(MEDIAN_MEAN, GSVA_MEAN)) +
  geom_point() +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") 

GSVA_MED_REGRESSION <- lm(GSVA_MEAN ~ MEDIAN_MEAN, data = GSVA_and_MEDIAN)
summary(GSVA_MED_REGRESSION)



