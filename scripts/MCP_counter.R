#MCP-COUNTER
set.seed(123)
source("./scripts/Packages.R")
source("./scripts/Functions.R")

gene_signature <- read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/genes_MCP_counter.txt")
gene_signature %>% write_xlsx("./results/TCP_counter_genes.xlsx")
str(gene_signature)

#install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
library(devtools)
library(curl)
#install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)

##obrir dataframe
data_mrna_seq_v2_rsem <- 
  read.delim("C:/Users/FixeI/OneDrive/Bureau/Anna/R/TCGA_lung_c_practice/raw_data/data_mrna_seq_v2_rsem.txt", header = T)

expression <- data_mrna_seq_v2_rsem %>%
  dplyr::distinct(Hugo_Symbol, .keep_all = T) %>%
  column_to_rownames("Hugo_Symbol") %>%
  .[, -which(colnames(expression)== "Entrez_Gene_Id")]
  
str(expression)
expression <- log2(expression + 1)

#MCP_counter
MCPcounter <- MCPcounter.estimate(
  expression,featuresType=c("HUGO_symbols")[1],
  probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
  genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
)

?MCPcounter::MCPcounter.estimate
library(pheatmap)
MCPcounter %>%
  as.matrix()%>%
  pheatmap(
    cluster_rows = T, cluster_cols = T,
    scale = "row",
    #treeheight_col = 0,
    main = "MCP counter heatmap",
    clustering_method = 'ward.D2',
    clustering_distance_cols = 'euclidean',
    clustering_distance_rows = 'euclidean',
    show_colnames = F, show_rownames = T,
    fontsize_row = 10,
    #fontsize_col = 1,
    cutree_cols  = 2,
    color = colorRampPalette(c("blue", "white","red"))(100)
    ) 
#other version of heatmap
heatmap(as.matrix(MCPcounter),col=colorRampPalette(c("blue","white","red"))(100))
