#### OPEN FILE AND PREPARE IT ####
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")


data_mrna_seq_v2_rsem <- 
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", header = T)
View(data_mrna_seq_v2_rsem)


#take columns you need
TLS_signature <- 
  c("CD79A", "FCRL5", "SSR4", "XBP1", 
    "IL7R", "CXCL12", "LUM", "C1QA", "C7",
    "CD52", "APOE", "PTGDS", "PIM2", "DERL3",
    "CCL19", "CCL21", "CXCL13", "CCL17", "CCL22", 
    "IL16", "ICAM2", "ICAM3", "VCAM1", "MADCAM1", 
    "ITGAL", "ITGA4", "ITGAD", "LTB",  "CD37")

data <- subset(data_mrna_seq_v2_rsem, Hugo_Symbol %in% TLS_signature)
View(data)

#put correct names and eliminate the columns which you don't need, transform
rownames(data) <- data[,1]
dat <- 
  data[, -which(colnames(data) == "Hugo_Symbol")] 
dat <-
  dat[, -which(colnames(dat)== "Entrez_Gene_Id")]
dat <- log2(dat + 1)

save(dat,file="data_gene_expression.txt")

#### DENSITY PLOT ####

#fer la mediana
med <- median(apply(dat, 2, median))

dat %>% 
  as_tibble() %>%
  colMeans() %>% as.data.frame() %>%
  # # dplyr::rename("mean_expression" = ".") %>%
  ggplot2::ggplot(aes(x = .)) +
  geom_density(fill = "pink", alpha = 0.3) + 
  geom_vline(aes(xintercept = med)) +
  geom_text(aes(x = 10, label = "TLS median signature"), y= 0.45) +
  geom_text(aes(x = 7, label = "TLS LOW"), y= 0.35) +
  geom_text(aes(x = 11, label = "TLS HIGH"), y= 0.35) +
  xlab("TLS signature mean expression") +
  ylab("Density") +
  my_theme

ggsave("TLS signature mean expression.pdf")

#### HEATMAP OF EXPRESSION ####
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)


dat %>% 
  as.matrix() %>%
  pheatmap( 
    cluster_rows = F, cluster_cols = T,
    scale = "row",
    #treeheight_col = 0,
    main = "Gene Expression Heatmap",
    clustering_method = 'ward.D',
    clustering_distance_cols = 'euclidean',
    clustering_distance_rows = 'euclidean',
    show_colnames = F, show_rownames = T,
    fontsize_row = 7,
    cutree_cols  = 2,
    color = colorRampPalette(c("blue", "white","red"))(50)) 


