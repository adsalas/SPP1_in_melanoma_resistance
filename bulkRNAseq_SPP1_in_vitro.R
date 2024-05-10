#####  Analysis of bulkRNAseq in vitro experiment for SPP1 effect on Melanoma cell lines  #######

# Load and preprocess the datasets ----
# Load the libraries
library(tidyverse)
library(data.table)
library(tools)
library(pheatmap)
# Define the path to the count tables for each sample
path <- "~/Documents/Projects/Jelena/Revision/RNAseq_in_vitro/mapping_outputs_43"
# Store the names of the paths:
files <- dir(path = path, pattern = "*ReadsPerGene.out.tab", full.names = T, recursive = T) 
# Create a combined table skipping the first 4 lines which contain non-required information 
counttablefull <- files %>%
  map(read_tsv,  
      skip = 4, 
      col_names = FALSE ) %>%
  reduce(cbind)

# Check the data
counttablefull[1:5, 1:10]
colnames(counttablefull)

# Extract and simplify the samples names
datasets <-
  files %>%
  stringr::str_replace("/Users/asalas/Documents/Projects/Jelena/Revision/RNAseq_in_vitro/mapping_outputs_43/", "") %>% 
  stringr::str_replace("_ReadsPerGene.out.tab", "")

datasets
datasets <- sapply(strsplit(datasets,"/"), function(x) x[2])
datasets

# Format the colum names
columnnames <- c()
for (i in datasets) {
  columnnames <- c(columnnames,
                   paste0("gene_", i),
                   paste0(i, "-unstranded"),
                   paste0(i, "-forwardstrand"),
                   paste0(i, "-reversestrand")
  )
}

names(counttablefull) <- columnnames
rm(columnnames, datasets)
head(counttablefull)
colnames(counttablefull)

# Make a column for gene identifier, and remove duplicated gene id information:
counttablefull <- counttablefull %>%
  mutate(ensembl_gene = gene_43_CD44_1_S94) %>%
  select(-starts_with("gene"))

counttablefull %>% head()

# Check for indication of strandness
counttablefull %>% # remove gene name column
  select(-ensembl_gene) %>% # find column sum
  summarise_each(funs(sum)) %>%  # gather wide data (columns) into a long table
  gather(library, counts) %>% # split library column into dataset and protocol 
  separate(library, into = c("sample", "stranding"), sep="-") %>% # split library into 2 columns
  spread(stranding,counts) %>% # make columns for each of the strandings
  mutate(propF = forwardstrand/unstranded, propR = reversestrand/unstranded) # assess rev/unst and for/unst

# Keep the reverse stranded counts only
counttablefull.FTD <- counttablefull %>% select(ends_with("reversestrand"), ends_with("gene"))
dim(counttablefull.FTD)
head(counttablefull.FTD)

# Convert the Ensembl IDs to common gene names with bioMart
# Load the libraries
library(biomaRt)
# Define the mart to use
mart <- useMart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl")

mygenesEnsem <- counttablefull.FTD$ensembl_gene
mygenesEnsem
length(mygenesEnsem)
length(unique(mygenesEnsem))

# Make the query of the required information
ext.gene.names <- getBM(attributes = c('ensembl_gene_id',
                                       'external_gene_name'),
                        filters = 'ensembl_gene_id', 
                        values = mygenesEnsem,
                        mart = mart)

dim(ext.gene.names)
length(unique(ext.gene.names$external_gene_name))

# Add to the dataframe the information of Ensembl IDs and External Gene Names, keeping the Ensembl IDs that don't have a match
counttablefull.FTD <- counttablefull.FTD[order(counttablefull.FTD$ensembl_gene, decreasing = TRUE), ]
counttablefull.FTD$ensembl_gene_id <- counttablefull.FTD$ensembl_gene

ext.gene.names <- ext.gene.names[order(ext.gene.names$ensembl_gene_id, decreasing = TRUE), ]
head(ext.gene.names)
tail(ext.gene.names)

head(counttablefull.FTD)
head(ext.gene.names)

# Join the count table with the external gene name dataframe by "ensembl_gene_id"
tail(left_join(counttablefull.FTD, ext.gene.names, by= "ensembl_gene_id"))
counttablefull.FTD <- left_join(counttablefull.FTD, ext.gene.names, by= "ensembl_gene_id")
head(counttablefull.FTD)

# Check uniqueness of the ensembl IDs and transfer them as row names 
dim(counttablefull.FTD)
length(unique(counttablefull.FTD$ensembl_gene))
rownames(counttablefull.FTD) <- counttablefull.FTD$ensembl_gene
counttablefull.FTD[1:10, 1:4]

# Store the ensembl IDs and matching gene symbols in another table
colnames(counttablefull.FTD)
genes.counttablefull.FTD <- counttablefull.FTD %>% dplyr::select(contains("gene"))
head(genes.counttablefull.FTD)
counttablefull.FTD <- counttablefull.FTD %>% dplyr::select(-contains("gene"))
head(counttablefull.FTD)

# Make the column names of the count table match the names of the samples in the metadata table
colnames(counttablefull.FTD) <- sapply(strsplit(colnames(counttablefull.FTD),"-"), function(x) x[1])
head(counttablefull.FTD)
dim(counttablefull.FTD)

# Perform DESeq2 analysis for DGE ----
library(DESeq2)
library(RColorBrewer)

# Create the coldata table
coldata <- as.data.frame(colnames(counttablefull.FTD))
colnames(coldata) <- "Sample"
rownames(coldata) <- coldata$Sample
# Create a new columns for the group and replicates
coldata$Sample <- str_replace(coldata$Sample, pattern = "1_C", "1.C") 
coldata$Group <- sapply(strsplit(coldata$Sample,"_"), function(x) x[2])
coldata$replicate <- sapply(strsplit(coldata$Sample,"_"), function(x) x[3])
coldata$Sample_new <- paste0(coldata$Group, "_", coldata$replicate)
coldata
# Convert to factors the different entries
coldata$Group <- as.factor(coldata$Group)
coldata$replicate <- as.factor(coldata$replicate)
# Change the column names of the count table
all(rownames(coldata) == colnames(counttablefull.FTD))
rownames(coldata) <- coldata$Sample_new
colnames(counttablefull.FTD) <- coldata$Sample_new
all(rownames(coldata) == colnames(counttablefull.FTD))

# Create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = counttablefull.FTD,
                                 colData = coldata,
                                 design = ~ Group)

dds
mcols(dds)

# Pre-filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# Check the levels name
unique(dds$Group)
# Indicate the reference level
dds$Group <- relevel(dds$Group, ref = "MEK")

# Perform the differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

plotMA(res, ylim=c(-12,12))

# Data transformation and visualization ----
vsd.dds <- vst(dds, blind=FALSE)
head(assay(vsd.dds), 3)

# Check the distribution of the samples in the PCA representation
plotPCA(vsd.dds, intgroup="Group") + ggtitle("Cell line M150543")

# Check the heatmap of sample to sample distance
sampleDists <- dist(t(assay(vsd.dds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.dds$Group, vsd.dds$replicate, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Make the plot
pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors)

# Define the contrasts of interest ---- 
# Load the libraries
library(ComplexHeatmap)
library(circlize)

# Retrieve the model matrix
mod.mat <- model.matrix(design(dds), colData(dds))
# Calculate coefficient vectors for each group
unique(dds$Group)
CD44 <- colMeans(mod.mat[dds$Group == "CD44", ])
MEK <- colMeans(mod.mat[dds$Group == "MEK", ])
SPP1 <- colMeans(mod.mat[dds$Group == "SPP1", ])
SPP1.CD44 <- colMeans(mod.mat[dds$Group == "SPP1.CD44", ])
# Obtain results for each pairwise contrast of interest
res_CD44_MEK <- results(dds, contrast = CD44 - MEK) # Effect of CD44
res_SPP1_MEK <- results(dds, contrast = SPP1 - MEK) # Effect of SPP1
res_SPP1.CD44_MEK <- results(dds, contrast = SPP1.CD44 - MEK) # Effect of SPP1.CD44

# Filter the results based on a padj threshold
res_CD44_MEK_05padj <- filter(as.data.frame(res_CD44_MEK), padj < 0.05)
dim(res_CD44_MEK_05padj)
res_SPP1_MEK_05padj <- filter(as.data.frame(res_SPP1_MEK), padj < 0.05)
dim(res_SPP1_MEK_05padj)
res_SPP1.CD44_MEK_05padj <- filter(as.data.frame(res_SPP1.CD44_MEK), padj < 0.05)
dim(res_SPP1.CD44_MEK_05padj)

# Comparison of the upregulated DEGs sets for each condition vs MEK as upset plots ----
# Define the up-regulated gene sets from each treatment
genes.UP.CD44 <- res_CD44_MEK_05padj %>% filter(log2FoldChange > 0 & padj < 0.05) %>% rownames() # No filtering of genes based on log2FC magnitude
genes.UP.SPP1 <- res_SPP1_MEK_05padj %>% filter(log2FoldChange > 0 & padj < 0.05) %>% rownames() # No filtering of genes based on log2FC magnitude
genes.UP.SPP1.CD44 <- res_SPP1.CD44_MEK_05padj %>% filter(log2FoldChange > 0 & padj < 0.05) %>% rownames() # No filtering of genes based on log2FC magnitude
# Combine the sets in a list
list.genes <- list("Up CD44 vs MEKi" = genes.UP.CD44,
                        "Up SPP1 vs MEKi" = genes.UP.SPP1,
                        "Up SPP1.CD44 vs MEKi" = genes.UP.SPP1.CD44
                        )

# Create a binary matrix from the list
matrix.genes <- list_to_matrix(list.genes)
dim(matrix.genes)
m <- make_comb_mat(matrix.genes)

# Define the colors
my.col <- brewer.pal(12, "Set3")
my.col <- my.col[c(1,3,4,5,6,7,8)]
# Make tue UpSet plot
UpSet(m = m, 
      set_order = c("Up SPP1 vs MEKi", 
                    "Up CD44 vs MEKi",
                    "Up SPP1.CD44 vs MEKi"),
      comb_col = my.col,
      column_title = "DEGs (padj < 0.05 & Log2FC > 0)"
)


# Check GO terms associated with the DEGs uniquely Up in SPP1 treatment ----
# Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Define the Universe
Allgenes <- rownames(as.data.frame(res_SPP1_MEK))
names(Allgenes) <- Allgenes
length(Allgenes)

# Explore the matrix
m
# Retrieve and define the gene set I want to look for enrichment
my.genes.up <- extract_comb(m, "010")
names(my.genes.up) <- my.genes.up
length(my.genes.up)

# Perform over-representation test
ego.genes.UP.SPP1.unique <- enrichGO(gene = names(my.genes.up), 
                                     universe = names(Allgenes),
                                     OrgDb = org.Hs.eg.db, 
                                     ont = "BP",
                                     pAdjustMethod = "BH", 
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05, 
                                     minGSSize = 10,
                                     maxGSSize = 500, 
                                     keyType = "ENSEMBL")


# Network-like representation for Apoptosis-linked GO terms ----
# Load the libraries
library(igraph)

# Check for GO terms associated with apoptosis
ego.genes.UP.SPP1.unique@result[grep("apopt", ego.genes.UP.SPP1.unique@result$Description, value = F), ]
# Create a dataframe with all the GO terms of interest, the Ensembl gene IDs and common gene names
data1 <- ego.genes.UP.SPP1.unique@result[grep("apopt", ego.genes.UP.SPP1.unique@result$Description, value = F), ]
# Filtering based on defined threshold
dim(data1)
data1 <- data1 %>% filter(p.adjust < 0.05)
dim(data1)

# Retrieve the Ensembl IDs for each individual GO term inside the categories of interest
my.genes <- data1$geneID
length(my.genes)
my.genes.symbol.df <- list()
# Populate the dataframes with a new column containing the gene name
for (i in 1:length(my.genes)) {
  
  genes.Ensembl <- strsplit(my.genes[i], "/") %>% unlist()
  
  my.genes.symbol.df[[i]] <- genes.counttablefull.FTD[genes.Ensembl, ]
  my.genes.symbol.df[[i]]$GO_term <- data1$Description[i]
  my.genes.symbol.df[[i]]$ID <- data1$ID[i]
  
}

# Combine the dataframes inside the list
my.df <- do.call("rbind", my.genes.symbol.df)

# Filter by specific processes of interest
my.df <- my.df[my.df$GO_term %in% c(my.df$GO_term[1], unique(grep("negative", my.df$GO_term, value = T))) , ]

# Format the data for the network representation
#### LINKS   #####
data.links <- my.df[ , c("ID", "GO_term", "external_gene_name")]
colnames(data.links)
colnames(data.links) <- c("Process_ID", "Process", "Gene")
data.links$from <- data.links$Process
data.links$to <- data.links$Gene
colnames(data.links)
data.links <- data.links[ , c(4,5,1,2,3)]
data.links

# Add the Processes as nodes in the table
unique(data.links$Process)
data.temp <- as.data.frame(matrix(data= NA,
                                  nrow = length(unique(data.links$Process)), 
                                  ncol = 5))
colnames(data.temp) <- colnames(data.links)
data.temp$from <- unique(data.links$Process)
data.temp$to <- data.temp$from[c(2:4, 1)] # with padj < 0.05 and including the general "apoptotic signaling pathway" term and terms with "negative" in the name
data.temp

data.links <- rbind(data.links, data.temp)
rownames(data.links) <- NULL
# Add a weight column so I can use it for "cutting the links" for the graphical representation
data.links$weight <- 4
data.links[data.links$to %in% data.temp$to, ]$weight <- 1
# Define the color codes
my.col <- c("steelblue1",
            "sienna2",
            "salmon1",
            "orange")

# Vertex colors
data.links$color <- my.col[1]
unique(data.links$from)
for (i in 1:length(data.temp$from)) {
  data.links[data.links$from == data.temp$from[i], ]$color <- my.col[i]  
}

data.links
tail(data.links)

# Edges colors
data.links[data.links$to %in% data.temp$to, ]$color <- NULL
# Vertex label colors
data.links$vertex.label.color <- "gray31"

# Vertex font format
data.links$vertex.label.font <- 1

# Vertex label cex
data.links$vertex.label.cex <- 0.5

data.links[data.links$to  %in% data.temp$to, ]$vertex.label.color <- "black" 

#### NODES   #####
data.nodes <- as.data.frame(unique(data.links$to))
colnames(data.nodes) <- "id"

# Define the colors for the nodes
data.nodes$color <- "gainsboro"
unique(data.nodes$id)
for (i in 1:length(data.temp$from)) {
  data.nodes[data.nodes$id == data.temp$from[i], ]$color <- my.col[i]  
}
# Defining which are going to be the labels
data.nodes$vertex.label.2 <- data.nodes$id
data.nodes[data.nodes$id %in% data.temp$to, ]$vertex.label.2 <- NA
# Define the sizes of the nodes
# Based on p.adj (more significant bigger sizes)
data.nodes$vertex.size.2 <- 2
data.nodes[data.nodes$id == data1$Description[1], ]$vertex.size.2 <- -log2(data1$p.adjust[1])
data.nodes[data.nodes$id == data1$Description[3], ]$vertex.size.2 <- -log2(data1$p.adjust[3])
data.nodes[data.nodes$id == data1$Description[7], ]$vertex.size.2 <- -log2(data1$p.adjust[7])
data.nodes[data.nodes$id == data1$Description[10], ]$vertex.size.2 <- -log2(data1$p.adjust[10])

# Defining the position of the labels
data.nodes$vertex.label.dist <- 0.7
data.nodes[data.nodes$id %in% data.temp$to, ]$vertex.label.dist <- 0

# Create the net with the formatted data
net2 <- graph_from_data_frame(d=data.links, 
                              vertices=data.nodes, 
                              directed=T) 
net2

# Remove the connection between the GO terms nodes
net2 <- delete_edges(net2, E(net2)[weight < 2])

# Make the representation
set.seed(22)
par(mar=c(0,0,1.4,0))
plot(net2, 
     layout= layout_with_dh, 
     rescale=T, 
     vertex.color= V(net2)$color,
     vertex.frame.color = "white",
     vertex.label.family= "Helvetica",
     vertex.label=  V(net2)$vertex.label.2,
     vertex.label.font = V(net2)$vertex.label.font,
     edge.arrow.width = 1.5,
     vertex.label.degree = pi/2,
     vertex.label.dist= V(net2)$vertex.label.dist,
     vertex.label.color= "black",
     edge.lty= 6, 
     vertex.size = V(net2)$vertex.size.2,
     edge.color= E(net2)$color,
     edge.arrow.mode=0,
     edge.width = 1,
     asp = 2/3.5, # y/x
     vertex.label.cex= 0.8, 
     edge.arrow.size=.25)

# Plot the legends
# Legend GO terms
legend(x=-1, 
       y= -0.8,
       legend = unique(data.temp$from),
       title = "GO term",
       pch=21, 
       pt.bg= my.col[c(1:10)],
       pt.cex= 2, 
       y.intersp = 1.1,
       cex=.8, 
       bty="n", 
       ncol=1)

# Legend significance level
legend(x= -1, 
       y= -0.4,
       legend = c("     low","     high"), 
       pch=21,
       text.col = "black",
       title.col = "black",
       title = "Significance",
       pt.bg= "gray",
       pt.cex= -log2(data1$p.adjust[c(10,1)]) *2.5/4,
       cex=.8,
       xjust = 0,
       y.intersp = c(1.2, 1.8),
       bty="n",
       ncol=1)


# Heatmap with the genes belonging to the GO terms of interest ----
# Extract the information of the genes that are related with the GO term
# negative regulation of apoptotic signaling pathway
neg.reg.apop.sign.path <- strsplit(filter(ego.genes.UP.SPP1.unique@result, Description == "negative regulation of apoptotic signaling pathway") %>% pull(geneID),"/")
neg.reg.apop.sign.path <- unlist(neg.reg.apop.sign.path)
length(neg.reg.apop.sign.path)
neg.reg.apop.sign.path.df <- genes.counttablefull.FTD[genes.counttablefull.FTD$ensembl_gene %in% neg.reg.apop.sign.path, ]
length(unique(neg.reg.apop.sign.path.df$external_gene_name))
neg.reg.apop.sign.path.Symbol <- neg.reg.apop.sign.path.df$external_gene_name
length(neg.reg.apop.sign.path.Symbol)

# Define the dataset
df_norm <- counts(dds, normalized= T)[, ]
all(colnames(df_norm) == rownames(coldata))
dim(df_norm)
data <- df_norm
data[1:5 , ]

# Define the annotations
group <- coldata$Group
group <- factor(group, levels = c("MEK","CD44","SPP1","SPP1.CD44")) 

# Annotations
column_ha <- HeatmapAnnotation(group = group,
                               col = list(group = c("MEK" = my.col[1],
                                                    "CD44" = my.col[2],
                                                    "SPP1" = my.col[3],
                                                    "SPP1.CD44" = my.col[4])
                               ))

# Filter the genes of interest
data.zsc <- data[neg.reg.apop.sign.path, ]
dim(data.zsc)
# Perform z-score transformation
data.zsc <- t(as.data.frame((scale(t(data.zsc)))))
# Create a new variable to split the columns
split.column.order <- group

# Plot the heatmap
set.seed(1234)
heatm.neg.reg.apop.sign.path <- Heatmap(data.zsc, 
                                   show_row_names = T, 
                                   show_column_names = F, 
                                   row_labels = neg.reg.apop.sign.path.Symbol, 
                                   column_title = NULL,
                                   name = "zscore", 
                                   cluster_rows = T,
                                   show_row_dend = F,
                                   row_names_gp = gpar(fontsize = 9),
                                   row_title = "GO term: negative regulation of apoptotic signaling pathway",
                                   row_title_gp = gpar(fontsize = 11),
                                   cluster_columns = F, 
                                   column_split = split.column.order, 
                                   top_annotation = column_ha, 
                                   border = T,
                                   use_raster = T
) 

heatm.neg.reg.apop.sign.path

# Retrieve the session information
sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: OS X  13.3
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      parallel  stats4    tools     stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] igraph_1.2.6                org.Hs.eg.db_3.11.4         AnnotationDbi_1.50.3        clusterProfiler_3.16.1      circlize_0.4.13             ComplexHeatmap_2.4.3       
# [7] RColorBrewer_1.1-2          DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.58.0          Biobase_2.48.0             
# [13] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0         biomaRt_2.44.4             
# [19] pheatmap_1.0.12             data.table_1.14.0           forcats_0.5.1               stringr_1.4.0               dplyr_1.0.5                 purrr_0.3.4                
# [25] readr_1.4.0                 tidyr_1.1.3                 tibble_3.1.1                ggplot2_3.3.3               tidyverse_1.3.1            
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1           backports_1.2.1        fastmatch_1.1-3        BiocFileCache_1.12.1   plyr_1.8.6             splines_4.0.2          BiocParallel_1.22.0   
# [8] urltools_1.7.3         digest_0.6.27          yulab.utils_0.0.4      GOSemSim_2.14.2        viridis_0.6.0          GO.db_3.11.4           fansi_0.4.2           
# [15] magrittr_2.0.1         memoise_2.0.0          cluster_2.1.2          annotate_1.66.0        graphlayouts_0.7.1     modelr_0.1.8           askpass_1.1           
# [22] enrichplot_1.8.1       prettyunits_1.1.1      colorspace_2.0-0       blob_1.2.1             rvest_1.0.0            rappdirs_0.3.3         ggrepel_0.9.1         
# [29] WriteXLS_6.4.0         haven_2.4.1            crayon_1.4.1           RCurl_1.98-1.3         jsonlite_1.7.2         scatterpie_0.1.7       genefilter_1.70.0     
# [36] survival_3.2-11        glue_1.4.2             polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.34.0        XVector_0.28.0         GetoptLong_1.0.5      
# [43] shape_1.4.6            scales_1.1.1           DOSE_3.14.0            DBI_1.1.1              Rcpp_1.0.6             viridisLite_0.4.0      xtable_1.8-4          
# [50] progress_1.2.2         clue_0.3-60            gridGraphics_0.5-1     bit_4.0.4              europepmc_0.4.1        httr_1.4.2             fgsea_1.14.0          
# [57] ellipsis_0.3.1         pkgconfig_2.0.3        XML_3.99-0.6           farver_2.1.0           dbplyr_2.1.1           locfit_1.5-9.4         utf8_1.2.1            
# [64] labeling_0.4.2         ggplotify_0.1.0        tidyselect_1.1.0       rlang_0.4.10           reshape2_1.4.4         munsell_0.5.0          cellranger_1.1.0      
# [71] cachem_1.0.4           downloader_0.4         cli_3.6.2              generics_0.1.0         RSQLite_2.2.7          broom_0.7.6            ggridges_0.5.3        
# [78] fastmap_1.1.0          bit64_4.0.5            fs_1.5.0               tidygraph_1.2.0        ggraph_2.0.5           DO.db_2.9              xml2_1.3.2            
# [85] compiler_4.0.2         rstudioapi_0.13        png_0.1-7              curl_4.3               reprex_2.0.0           tweenr_1.0.2           geneplotter_1.66.0    
# [92] stringi_1.5.3          lattice_0.20-41        Matrix_1.3-2           vctrs_0.3.7            pillar_1.6.0           lifecycle_1.0.0        BiocManager_1.30.12   
# [99] triebeard_0.3.0        GlobalOptions_0.1.2    cowplot_1.1.1          bitops_1.0-7           qvalue_2.20.0          R6_2.5.0               gridExtra_2.3         
# [106] MASS_7.3-53.1          assertthat_0.2.1       rjson_0.2.20           openssl_1.4.3          withr_2.4.2            GenomeInfoDbData_1.2.3 hms_1.0.0             
# [113] ggfun_0.0.9            rvcheck_0.2.1          ggforce_0.3.3          lubridate_1.7.10

