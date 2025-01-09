rm(list=ls())

# Load libraries
library(tidyverse)
library(Seurat)
library(gridExtra)
library(ggpubr)
library(SingleR)
library(celldex)
library(AUCell)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the dataset and UMI cutoffs
base_dir = "/scratch/tu.tra/FinalProject-BINF6430/Data/"
samples <- c("P1", "P2", "P3", "P4", "P5", "P6", "H1", "H2", "H3")
file_path <- paste0(base_dir, samples, "_before")
file_path[7:9] <- gsub("_before", "", file_path[7:9]) #Remove "_before" for H1, H2, H3

low_cutoff <- c(1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000)
high_cutoff <- c(60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000)


data.seurat.obj <- lapply(seq_along(samples), function(i) {
  # Read data
  sample_data <- Read10X(data.dir = file_path[i])
  # Create Seurat objects 
  sample_data <- CreateSeuratObject(counts = sample_data,
                     min.cells = 3, 
                     min.features = 200,
                     assay = "RNA")
  # Quality Control
  sample_data[["percent.mt"]] <- PercentageFeatureSet(sample_data, pattern = "^MT-")
  # Filtering based on QC metrics (remove unwanted cells from the data)
  sample_data <- subset(sample_data,
                            subset = nFeature_RNA > 200 
                            & nFeature_RNA < 2500 
                            & percent.mt < 5
                            & nCount_RNA > low_cutoff[i]
                            & nCount_RNA < high_cutoff[i]
                              )
})

# Normalize the data
data.seurat.obj <- lapply(data.seurat.obj, function(sample_data) {
  print(dim(sample_data))
  sample_data <- NormalizeData(sample_data, normalization.method = "LogNormalize", scale.factor = 10000)
  sample_data <- FindVariableFeatures(sample_data, selection.method = "vst", nfeatures = 2000)
})

# Identify high variable features (only select features that show high cell-to-cell variation)
top10_list <- lapply(data.seurat.obj, function(sample_data) {
  top10 <- head(VariableFeatures(sample_data), 10)
  })

# Name the top 10 list to match samples 
names(top10_list) <- samples


# Integrate data
# Integrated data has a cell-level meta data 
anchors <- FindIntegrationAnchors(object.list = data.seurat.obj, dims = 1:20)
comb_data <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(comb_data) <- "integrated"


# Scale the data
comb_data <- ScaleData(comb_data, verbose = TRUE)

# Perform linear dimensionality reduction
comb_data <- RunPCA(comb_data)
comb_data <- RunUMAP(comb_data, reduction = "pca", dims = 1:30)

# Get the number of cells from each original sample
cell_counts <- sapply(data.seurat.obj, ncol)

# Create a vector of sample names repeated for their respective cell counts
sample_ids <- rep(samples, cell_counts)

# Add this information to the combined object
comb_data$sample <- sample_ids
comb_data$condition <- ifelse(comb_data$sample %in% c("H1", "H2", "H3"), "Healthy", "KD_Acute")

# CLustering (cluster cells that have similar patterns)
comb_data <- FindNeighbors(comb_data, reduction = "pca", dims = 1:30)
comb_data <- FindClusters(comb_data, resolution = 0.8)

# Look at cluster IDs of the first 5 cells
head(Idents(comb_data), 5)

# Visualize clusters
plot1 <- DimPlot(comb_data, reduction = "umap", label = TRUE)

# Visualize clusters splitting by condition
plot2 <- DimPlot(comb_data, reduction = "umap", split.by = "condition")

# Visualize clusters grouping by sample type
plot3 <- DimPlot(comb_data, reduction = "umap", group.by = "sample") + ggtitle("")

# Visualize clusters grouping by sample
plot4 <- DimPlot(comb_data, reduction = "umap", group.by = "condition") + ggtitle("")

# Cell annotations by SingleR
# extract log-count expression matrix
comb_data[["RNA"]] <- JoinLayers(comb_data[["RNA"]])
exp <- comb_data[["RNA"]]$data

# Load internal reference of human immune cells for SingleR
refs<-list(ref1 = HumanPrimaryCellAtlasData(), 
           ref2 = BlueprintEncodeData(),
           ref3 = DatabaseImmuneCellExpressionData(),
           ref4 = NovershternHematopoieticData(), 
           ref5 = MonacoImmuneData())

# Scale the data (just to make sure)
comb_data <- ScaleData(comb_data)

# Find markers that define clusters via differential expression (identifies positive and negative markers of a single cluster
# and compare with all other cells)
comb_data_markers <- FindAllMarkers(comb_data, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)

# select top 10 markers per cluster 
top2_marker <- comb_data_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Because all(top10_marker$gene %in% rownames(comb_data[["RNA"]]@layers$scale.data)) returns False, we have to filter the top10_marker
#top10_marker_filtered <- top10_marker[top10_marker$gene %in% rownames(comb_data[["RNA"]]@layers$scale.data), ]

ref_names <- list() 
# Perform annotation + generating heatmaps using each reference, respectively
for (ref_idx in names(refs)) {
  ref <- refs[[ref_idx]]
  ref_name <- paste("SingleR", ref_idx, sep = ".")
  ref_names[[ref_idx]] <- ref_name
  
  rst <- SingleR(test = exp, 
                 ref = ref,
                 clusters = comb_data[["seurat_clusters"]]$seurat_clusters, 
                 labels = ref$label.main)
  
  comb_data[[ref_name]] <- rst$pruned.labels[match(comb_data[["seurat_clusters"]]$seurat_clusters, rownames(rst))]
}

ref_name_1 <- ref_names[["ref1"]]

# One UMAP to visualize 
plot5 <- DimPlot(comb_data, 
                 reduction = "umap", 
                 group.by = ref_name_1, 
                 label = FALSE) + ggtitle("")

plot6 <- DimPlot(comb_data, 
                 reduction = "umap", 
                 group.by = ref_name_1,
                 split.by = "condition",
                 label = FALSE) + ggtitle("")

# Present the first 2 rows plot 1-4 together
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2, nrow = 3)

# Heatmap(figure 1b, page 5)
DoHeatmap(subset(comb_data,
                 downsample = 1000),
          features = top2_marker$gene,
          slot = "scale.data",
          group.by = ref_name_1,
          label = FALSE,
          size = 0.5)

comb_data$ref_name_1_labels <- comb_data[[ref_name_1]]

comb_data$combined_labels <- apply(comb_data@meta.data[,grep("SingleR", colnames(comb_data@meta.data))], 1, function(x) paste(unique(x), collapse = "/"))

# Stacked bar graph (Figure 1c, page 6)
# Calculate cell type percentages for each sample
cell_type_percentages <- comb_data@meta.data %>%
  group_by(sample, combined_labels) %>%
  summarise(count = n()) %>%
  mutate(percentage = count/sum(count) * 100) %>%
  ungroup()

# Create the stacked bar graph
ggplot(cell_type_percentages, aes(x = sample, y = percentage, fill = combined_labels)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Proportion", fill = "")


# Boxplot (figure 1d, page 7)
# Calculate cell type percentages for each condition
per_condition_percentages <- comb_data@meta.data %>%
  group_by(sample, combined_labels) %>%  # Group by sample instead of condition
  summarise(count = n()) %>%
  mutate(condition_percentage = count/sum(count) * 100) %>% 
  ungroup() %>%
  mutate(condition = case_when(
    sample %in% c("H1", "H2", "H3") ~ "Healthy",
    sample %in% c("P1", "P2", "P3", "P4", "P5", "P6") ~ "KD_Acute",
    TRUE ~ as.character(sample)
  ))

# Create box plot
ggplot(per_condition_percentages, aes(x = combined_labels, y = condition_percentage, fill = condition)) +
         geom_boxplot() +
         scale_fill_manual(values = c("Healthy" = "blue", "KD_Acute" = "red")) +
         labs(x = "Cluster", y = "Percentage", fill = "Group") +
  theme(
    axis.text.x = element_text(size = 4), # Smaller x-axis text
    axis.text.y = element_text(size = 5), # Smaller y-axis text
  ) +
  stat_compare_means(aes(group = condition), label = "p.format", 
                     label.y = max(per_condition_percentages$condition_percentage) * 1.1)
  
# Create Figure 2a. DEGs in cell clusters from samples (Expression levels of DEGs across cell clusters)
top2_marker <- comb_data_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

selected_gene <- unique(top2_marker$gene)

gene_expression <- GetAssayData(comb_data, assay = "RNA", layer = "data")[selected_gene, ]
average_expression <- rowMeans(gene_expression)
average_expression <- colMeans(gene_expression)
comb_data[["average_expression"]] <- average_expression

FeaturePlot(subset(comb_data,
                   downsample = 1000),
            features = "average_expression",
            reduction = "umap",
            min.cutoff = "q10",
            max.cutoff = "q90") + ggtitle("")

# Create Figure 2b AUC histogram
IRGs <- read.csv("/scratch/tu.tra/FinalProject-BINF6430/IRGs.csv")
IRGs <- as.vector(IRGs[,2])
IRGSet <- GeneSet(IRGs, setName="IRGs")
exprMatrix <- GetAssayData(comb_data, assay = "RNA", layer = "data")
cells_AUC <- AUCell_run(exprMatrix, geneSets=IRGSet)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

# Create Figure 2c. t-SNE plots of IRG scores in all clusters
# Calculate AUC scores (calculates the signature enrichment scores on each cell)
cells_AUC <- AUCell_run(exprMatrix, geneSets = IRGSet)
comb_data$IRG_score <- cells_AUC@assays@data$AUC["IRGs", ]
comb_data <- RunTSNE(comb_data, reduction = "pca", dims = 1:30)

# Store the labels from ref_name_1 as metadata
comb_data$ref_name_1_labels <- comb_data[[ref_name_1]]

# Extract UMAP coordinates
umap_data <- as.data.frame(comb_data[["umap"]]@cell.embeddings)
umap_data$cell_labels <- comb_data$ref_name_1_labels

# Aggregate UMAP data to get name of clusters once 
umap_data_agg <- umap_data %>%
  group_by(cell_labels) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

FeaturePlot(object = comb_data,
            features = "IRG_score",
            reduction = "umap",
            label = FALSE,
            label.size = 4) +
  geom_label(data = umap_data_agg, aes(x = umap_1, y = umap_2, label = cell_labels), check_overlap = TRUE) +
  labs(color = "AUC")


# GO Enrichment and KEGG Pathway enrichment
# Filter DEGs adjusted P < 0.05
degs_GO <- comb_data_markers %>%
  filter(p_val_adj < 0.05)

degs_gene_list <- degs_GO$gene

enrichcut <- 0.05

# Convert gene symbols to Entrez IDs 
map <- bitr(degs_gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

id <- map$ENTREZID


go_mf_results <- enrichGO(gene = id,
                          OrgDb = org.Hs.eg.db, 
                          ont = "MF", 
                          pvalueCutoff = enrichcut)

go_bp_results <- enrichGO(gene = id,
                          OrgDb = org.Hs.eg.db, 
                          ont = "BP", 
                          pvalueCutoff = enrichcut)

go_cc_results <- enrichGO(gene = id,
                          OrgDb = org.Hs.eg.db, 
                          ont = "CC", 
                          pvalueCutoff = enrichcut)


go_mf_results <- simplify(go_mf_results) # Remove redundancy
go_bp_results <- simplify(go_bp_results)
go_cc_results <- simplify(go_cc_results)

go_mf_results<-setReadable(go_mf_results, 'org.Hs.eg.db', 'ENTREZID')
go_bp_results<-setReadable(go_bp_results, 'org.Hs.eg.db', 'ENTREZID')
go_cc_results<-setReadable(go_cc_results, 'org.Hs.eg.db', 'ENTREZID') # convert ID to symbol

bp_df <- as.data.frame(go_bp_results)
mf_df <- as.data.frame(go_mf_results)
cc_df <- as.data.frame(go_cc_results) # Convert results to dataframe

bp_df$Ontology <- "Biological Process"
mf_df$Ontology <- "Molecular Function"
cc_df$Ontology <- "Cellular Component" # Add ontology type

combined_counts <- bind_rows(bp_df, mf_df, cc_df)

# Create a summary table for plotting
summary_table <- combined_counts %>% 
  group_by(Description, Ontology) %>%
  summarise(total_count = sum(Count), .groups = 'drop')

summary_table_subset <- head(summary_table, 30)

desired_order <- c("Molecular Function", "Cellular Component", "Biological Process")

# reorder the y-axis Description based on ontology
summary_table_subset$Description <- factor(summary_table_subset$Description,
                                           levels = summary_table_subset$Description[order(match(summary_table_subset$Ontology, desired_order))])


ggplot(summary_table_subset, aes(x = total_count, y = Description, fill = Ontology)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "The Most Enriched GO Terms",
    x = "Gene number",
    y = "Pathway",
    fill = "Type") +
  theme(axis.text.y = element_text(size = 9)) +
  scale_fill_manual(values = c("Biological Process" = "#66CC99",
                               "Molecular Function" = "#6098c8", 
                               "Cellular Component" = "#FF5733"))

# KEGG enrichment (Figure 2e)
kegg <- enrichKEGG(id, organism = 'hsa', pvalueCutoff = enrichcut)
kegg <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
kegg <- as.data.frame(kegg)

kegg_sorted <- kegg[order(kegg$p.adjust), ]

kegg_top20 <- head(kegg_sorted, 20)

ggplot(kegg_top20, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), 
                       fill = -log10(p.adjust))) +
  geom_bar(stat = "identity") +
  labs(x = "-log10(adj.P.value)", 
       y = "Pathway",
       title = "KEGG Enrichment Barplot") +
  scale_fill_gradient(low = "white", high = "darkred")


  









