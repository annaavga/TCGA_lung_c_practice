#### OPEN FILE AND PREPARE IT ####

data_mrna_seq_v2_rsem <- 
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", header = T)
View(data_mrna_seq_v2_rsem)

###libraries
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("dslabs")
library (dslabs)

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

#### CREATION OF MY THEME ####
my_theme <- theme(
  plot.title = element_text(size = rel(2)),
  panel.grid.major.y = element_line(color = 'gray90', linewidth = 0.3),
  panel.grid.minor.y = element_line(color = 'gray90', linewidth = 0.1),
  panel.grid.major.x = element_line(color = 'gray90', linewidth = 0.3),
  panel.grid.minor.x = element_line(color = 'gray90', linewidth = 0.1),
  panel.border = element_rect(color = 'gray50', fill = NA, linewidth = 1),
  plot.background = element_rect(fill = NULL),
  panel.background = element_rect(fill = 'gray99'),
  axis.line = element_line(color = 'gray50'),
  axis.text = element_text(color = 'gray40', face = 'bold'),
  axis.text.x = element_text(size = rel(1.2)),
  axis.text.y = element_text(size = rel(1.2)),
  axis.title = element_text(size = rel(1), face = 'bold'),
  axis.ticks = element_line(color = 'gray50'),
  legend.position = "right",
  legend.margin = margin(5, 5, 5, 5),
  legend.title = element_text(face = 'bold'),
  legend.background = element_blank(),
  legend.key = element_rect(fill = 'gray98'),
  legend.box.background = element_rect(fill = 'gray98',color = "gray50")
)

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


