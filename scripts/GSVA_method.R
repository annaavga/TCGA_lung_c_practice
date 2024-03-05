#### GSVA METHOD ####
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

# load data
data_mrna <- read.delim(
  "C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", 
  header = T)

GSVA_data <- data_mrna %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  `rownames<-`(.[,1]) %>%
  select(-Hugo_Symbol) %>% 
  select(-Entrez_Gene_Id)

GSVA_data <- log2(GSVA_data + 1)

str(GSVA_data)

#choose gene signature
TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

#prepare the df and gene signature to do the analysis
GSVA_data <- as.matrix(GSVA_data)
class(GSVA_data) <- "numeric"

TLS_signature <- as.list(TLS_signature)
class(TLS_signature)

##PERFORM GSVA
gsvaPar <- gsvaParam(GSVA_data , TLS_signature)
gsvaPar

gsva.es <- gsva(gsvaPar, verbose=F)
#dim(gsva.es)

gsva.es <- as.data.frame(gsva.es)
class(gsva.es)
###################fins aqui, mes no carrega

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
  #geom_vline(aes(xintercept = 0.2, col = "red")) +
  #geom_vline(aes(xintercept = -0.2, col = "red")) +
  #geom_text(aes(x = 0.15, label = "TLS median signature"), y= 0.45) +
  #geom_text(aes(x = 7, label = "TLS LOW"), y= 0.35) +
  #geom_text(aes(x = 11, label = "TLS HIGH"), y= 0.35) +
  xlab("GSVA score") +
  ylab("Density") +
  my_theme +
  theme(legend.position = "none")
density_plot_gsva
ggsave("GSVA score.pdf")

