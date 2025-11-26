# ============================================
# ncRNA Differential Expression Analysis
# LUAD Primary Tumors vs Brain Metastasis
# ============================================


library(limma)              # For differential expression analysis
library(GEOquery)           # For downloading GEO datasets
library(illuminaHumanv4.db) # For probe annotation mapping
library(biomaRt)            # For gene annotation and biotype info
library(ggplot2)            # For creating visualizations
library(ggrepel)            # For non-overlapping text labels
library(dplyr)              # For data manipulation
library(pheatmap)           # For creating heatmaps
library(RColorBrewer)       # For color palettes

# ============================================
#  Brain Metastasis Data
# ============================================
#  (E-MTAB-8659)
# This dataset contains 63 brain metastasis samples from LUAD patients
bm_data <- read.table(
  "~/Downloads/omics/E-MTAB-8659/processed_all.txt",
  header = TRUE,           # First row contains column names
  sep = "\t",              # Tab-separated file
  stringsAsFactors = FALSE, # Keep strings as characters
  check.names = FALSE      # Don't modify column names
)
View(bm_data) #so here each row is a probe or gene and each column is sample in our case tumor.


#how many genes are there in data set
length(bm_data$GeneSymbol)# 47259

#how many unique genes 
length(unique(bm_data$GeneSymbol))#34659


# Checking for duplicate Gene Symbol entries
dup_genes <- bm_data$GeneSymbol[duplicated(bm_data$GeneSymbol)] 
length(dup_genes)  # how many duplicates exist #12600
head(dup_genes)    # view the first few duplicate gene symbols
length(unique(dup_genes))#8454


######### WHY WE NEED TO REMOVE DUPLICATES ########################
#In micro array data it is common that multiple probes measure the same gene . 
#in our data set for example gene A1BG comes two times but the values are different for each row , if we don't clean this then the limma package will treat it as 2 different gene .
#so we use function avergae replicate probes-avereps() to average gene expression value of same gene across probes , this gives one clean row per gene symbol.
####################################################3


# Average duplicate probe measurements for the same gene
# This ensures one expression value per gene symbol
bm_data_clean <- avereps(bm_data[, -1], ID = bm_data$GeneSymbol)
View(bm_data_clean)
dim(bm_data_clean) #34,659 x 63


# ============================================
# Primary LUAD Tumor Data
# ============================================
# Download GSE60645 dataset from GEO database
# This contains primary lung tumor samples
gse60645 <- getGEO("GSE60645", GSEMatrix = TRUE)

# Extracting expression matrix (probe-level data)
primary_data_raw <- exprs(gse60645[[1]])

# Extracting  sample metadata (phenotype information)
pheno <- pData(gse60645[[1]])

# ============================================
# Filter for Adenocarcinoma Only
# ============================================
#GSE60645 contains multiple lung cancer subtypes
# We need only AC (adenocarcinoma = LUAD) to match the BM samples
pheno_ac <- pheno[pheno$`histological subtype:ch1` == "AC", ]
ac_samples <- rownames(pheno_ac)

cat("Total GSE60645 samples:", nrow(pheno), "\n")#117
cat("AC (LUAD) samples:", length(ac_samples), "\n")#77

# Keeping  only expression data for LUAD samples
primary_data_ac_only <- primary_data_raw[, colnames(primary_data_raw) %in% ac_samples]

# ============================================
#  Map Illumina Probes to Gene Symbols
# ============================================

#so in the microarray data we see rownames as gene names
#but in the control data here we see rownames as probeID so we have to change probe Id to gene name.
# We need to convert probe IDs to gene symbols for comparison
probes <- rownames(primary_data_ac_only)

# Map probes to HGNC gene symbols using Illumina annotation
gene_symbols <- mapIds(
  illuminaHumanv4.db,
  keys = probes,           # Input: probe IDs
  column = "SYMBOL",       # Output: gene symbols
  keytype = "PROBEID",     # Key type is probe ID
  multiVals = "first"      # If multiple genes per probe, take first
)


#okay so we see in gene symbols there are probes that could not me mapped to any known gene symbol so they return NA.
# we need to remove those also remove gene_symbols vector to remove same unmapped entries.

# Remove probes that didn't map to any gene
primary_data_ac_only <- primary_data_ac_only[!is.na(gene_symbols), ]
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

# Average expression for probes mapping to the same gene
# This gives one row per gene symbol
primary_data_avg <- avereps(primary_data_ac_only, ID = gene_symbols)

cat("Primary LUAD samples after filtering:", ncol(primary_data_avg), "\n")#77

# ============================================
# Retrieve ncRNA Annotations
# ============================================
# Using Ensembl BioMart to get gene biotype information
# This tells us which genes are ncRNAs vs protein-coding
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#so now we have to collect all the unique genes , from primary data set and LAUD data set and combine them into one vector and remove duplicates.
#to Filter only ncRNAs using gene_biotype value
# all unique genes from both datasets
all_genes <- unique(c(rownames(bm_data_clean), rownames(primary_data_avg)))

# Query Ensembl for gene information
gene_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", "description"),
  filters = "hgnc_symbol",
  values = all_genes,
  mart = ensembl
)

# Filter to keep only non-coding RNAs
# Exclude protein_coding genes to focus on regulatory RNAs
ncrna_info <- gene_info %>%
  filter(!gene_biotype %in% c("protein_coding"))
head(ncrna_info)


# ============================================
#Common ncRNAs Between Datasets
# ============================================
# Now we have to Identify ncRNAs present in brain metastasis data
ncrna_in_bm <- intersect(ncrna_info$hgnc_symbol, rownames(bm_data_clean))

# Identify ncRNAs present in primary tumor data
ncrna_in_primary <- intersect(ncrna_info$hgnc_symbol, rownames(primary_data_avg))

# Find ncRNAs present in BOTH datasets 
common_ncrna <- intersect(ncrna_in_bm, ncrna_in_primary)


cat("ncRNA in BM:", length(ncrna_in_bm), "\n")#1216
cat("ncRNA in Primary LUAD:", length(ncrna_in_primary), "\n")#2423
cat("Common ncRNAs:", length(common_ncrna), "\n")#1109

# ============================================
# Subset to Common ncRNAs
# ============================================
# Extracting  only common ncRNAs from brain metastasis data
bm_ncrna <- bm_data_clean[rownames(bm_data_clean) %in% common_ncrna, ]

# Extracting only common ncRNAs from primary tumor data
primary_ncrna <- primary_data_avg[rownames(primary_data_avg) %in% common_ncrna, ]

# Ensureing  exact same genes in same order in both datasets
common_genes <- intersect(rownames(bm_ncrna), rownames(primary_ncrna))
bm_ncrna <- bm_ncrna[common_genes, ]
primary_ncrna <- primary_ncrna[common_genes, ]

# ============================================
#  Combine and Normalize Data
# ============================================
# Merging  primary and BM data into single matrix
combined_expr <- cbind(primary_ncrna, bm_ncrna)

# Creating group labels for each sample
group <- c(rep("Primary", ncol(primary_ncrna)),
           rep("BM", ncol(bm_ncrna)))

# Applying quantile normalization to remove technical variation
# This makes expression distributions comparable between datasets
combined_expr_norm <- normalizeBetweenArrays(as.matrix(combined_expr), method = "quantile")

combined_expr_norm

# ============================================
# PCA to Check for Batch Effects
# ============================================
# Performing  principal component analysis on normalized data
# This visualizes dataset separation and potential batch effects
pca_result <- prcomp(t(combined_expr_norm))

# Create dataframe for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Dataset = group
)

# Plot PCA
ggplot(pca_df, aes(PC1, PC2, color = Dataset)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA: LUAD Primary vs Brain Metastasis",
       subtitle = paste("Primary LUAD n=", ncol(primary_ncrna), 
                        ", BM n=", ncol(bm_ncrna)))

# ============================================
#  Differential Expression Analysis
# ============================================
# Creating design matrix (no intercept model)
# This specifies which samples belong to which group
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- c("BM", "Primary")

# Fiting linear model to normalized expression data
# This estimates mean expression for each group
fit <- lmFit(combined_expr_norm, design)

# Defining contrast: BM vs Primary
# Positive logFC = higher in BM, Negative = higher in Primary
contrast.matrix <- makeContrasts(BM_vs_Primary = BM - Primary, levels = design)

# Applying contrast to fitted model
fit2 <- contrasts.fit(fit, contrast.matrix)

# Computing moderated t-statistics using empirical Bayes
# This shrinks variances toward common value for more stable inference
fit2 <- eBayes(fit2)

# Extracting results for all genes
# adj.P.Val is Benjamini-Hochberg adjusted p-value (FDR)
results <- topTable(fit2, coef = "BM_vs_Primary", number = Inf, adjust.method = "BH")

# ============================================
# Annotatation and Classify Results
# ============================================
# Classifying genes as significant based on thresholds:
# - adj.P.Val < 0.05 (5% FDR)
# - |logFC| > 1 (2-fold change)
results$significant <- ifelse(
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
  ifelse(results$logFC > 0, "Upregulated", "Downregulated"),
  "Not significant"
)

# Add gene symbols as column for merging
results$hgnc_symbol <- rownames(results)

# Merge with biotype annotations from Ensembl
results_annot <- merge(results, gene_info, by = "hgnc_symbol", all.x = TRUE)

# ============================================
#Summary
# ============================================

cat("Total ncRNAs tested:", nrow(results_annot), "\n") #1178

# Filter to significant genes only
sig_results <- subset(results_annot, adj.P.Val < 0.05 & abs(logFC) > 1)

cat("Significant ncRNAs:", nrow(sig_results), "\n")#706
cat("Upregulated in BM:", sum(sig_results$logFC > 0), "\n")#634
cat("Downregulated in BM:", sum(sig_results$logFC < 0), "\n\n")#72
cat("Biotype breakdown:\n")
print(table(sig_results$gene_biotype))

#lncRNA                              miRNA                           misc_RNA               processed_pseudogene                     protein_coding 
#81                                388                                  3                                 18                                  8 
#scaRNA                             snoRNA                              snRNA   transcribed_processed_pseudogene     transcribed_unitary_pseudogene 
#8                                116                                  3                                 13                                 14 
#transcribed_unprocessed_pseudogene             unprocessed_pseudogene                          vault_RNA 
#45                                  8                                  1 


write.csv(results_annot, "ncRNA_DEG_LUAD_FINAL.csv", row.names = FALSE)

# ============================================
# Volcano Plot Data
# ============================================
# volcano plot dataframe with transformed p-values
volcano_data <- results_annot %>%
  mutate(
    logP = -log10(adj.P.Val),  # Transform p-values for volcano plot
    Regulation = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated in BM",
      adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated in BM",
      TRUE ~ "Not significant"
    )
  )

# ============================================
# Volcano Plot - Upregulated Only
# ============================================
# Selecting top 20 upregulated ncRNAs for labeling
# These have lowest adjusted p-values among upregulated genes
top_genes_up <- volcano_data %>%
  filter(Regulation == "Upregulated in BM") %>%
  arrange(adj.P.Val) %>%
  head(20)

# volcano plot
p_volcano_up <- ggplot(volcano_data, aes(x = logFC, y = logP)) +
  # Plot all points colored by regulation status
  geom_point(aes(color = Regulation), alpha = 0.6, size = 1.8) +
  
  # Define colors: red for upregulated, blue for downregulated, grey for NS
  scale_color_manual(
    values = c(
      "Upregulated in BM" = "#E64B35FF",
      "Downregulated in BM" = "#4DBBD5FF",
      "Not significant" = "grey80"
    )
  ) +
  
  # Add vertical lines at logFC = ±1 (2-fold change thresholds)
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", 
             color = "black", linewidth = 0.4) +
  
  # Add horizontal line at adj.P.Val = 0.05
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "black", linewidth = 0.4) +
  
  # Add labels ONLY for top upregulated genes
  geom_text_repel(
    data = top_genes_up,
    aes(label = hgnc_symbol),
    size = 3.5,
    max.overlaps = 100,      # Allow many labels
    segment.size = 0.2,      # Thin connecting lines
    box.padding = 0.4,       # Space around labels
    point.padding = 0.3,     # Space around points
    show.legend = FALSE,     # Don't show these in legend
    force = 2,               # Repulsion strength
    min.segment.length = 0   # Always show connecting lines
  ) +
  
  # Add titles and axis labels
  labs(
    title = "Upregulated ncRNAs in LUAD Brain Metastasis",
    subtitle = paste("Top", nrow(top_genes_up), "labeled ncRNAs (adj.P < 0.05, log2FC > 1)"),
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P~value),
    color = ""
  ) +
  
  # Apply clean theme
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Display the plot
print(p_volcano_up)

ggsave("Volcano_Upregulated_ncRNA_LUAD_BM.png", p_volcano_up, 
       width = 10, height = 8, dpi = 600)

# ============================================
#Heatmap - Top 30 Upregulated ncRNAs
# ============================================
# Select top 30 upregulated ncRNAs (lowest adj.P.Val, logFC > 1)
top_up <- sig_results %>%
  filter(logFC > 1) %>%           # Only upregulated genes
  arrange(adj.P.Val) %>%          # Sort by significance
  head(30) %>%                    # Take top 30
  pull(hgnc_symbol)               # Extract gene names

# Extract normalized expression for these genes
heatmap_data <- combined_expr_norm[top_up, ]
rownames(heatmap_data) <- top_up

# Create column annotation showing Primary vs BM groups
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heatmap_data)

# Define colors for annotation bar
ann_colors <- list(Group = c(Primary = "#4DBBD5FF", BM = "#E64B35FF"))

# Plot heatmap with row scaling (z-score per gene)
# This shows relative expression patterns across samples
pheatmap(
  heatmap_data,
  annotation_col = annotation_col,      # Color bar for groups
  annotation_colors = ann_colors,       # Custom group colors
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  scale = "row",                        # Z-score normalization per gene
  show_rownames = TRUE,                 # Show gene names
  show_colnames = FALSE,                # Hide sample names (too many)
  fontsize_row = 8,                     # Gene name font size
  clustering_method = "complete",       # Hierarchical clustering method
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top 30 Upregulated ncRNAs in LUAD Brain Metastasis",
  border_color = NA                     # No borders around cells
)

png("Heatmap_Top30_Upregulated_ncRNA_LUAD_BM.png", 
    width = 2400, height = 2000, res = 300)


###########


#additional plots

# ============================================
# Bar Plot: ncRNA Biotype Distribution
# ============================================

# ---- 1. All ncRNAs tested ----
biotype_all <- results_annot %>%
  filter(!is.na(gene_biotype)) %>%
  count(gene_biotype, name = "Count") %>%
  arrange(desc(Count))

# Bar plot for all ncRNAs
p_all <- ggplot(biotype_all, aes(x = reorder(gene_biotype, Count), y = Count, fill = gene_biotype)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Distribution of ncRNA Biotypes ",
    x = "ncRNA Biotype",
    y = "Number of Genes"
  ) +
  theme(legend.position = "none")

ggsave("Barplot_ncRNA_Biotype_Distribution_All.png", p_all, width = 8, height = 6, dpi = 600)
print(p_all)


# ---- 2. Upregulated ncRNAs in Brain Metastasis ----
biotype_up <- sig_results %>%
  filter(logFC > 1) %>%
  filter(!is.na(gene_biotype)) %>%
  count(gene_biotype, name = "Count") %>%
  arrange(desc(Count))

# Bar plot for upregulated ncRNAs
p_up <- ggplot(biotype_up, aes(x = reorder(gene_biotype, Count), y = Count, fill = gene_biotype)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Upregulated ncRNA Biotypes in LUAD Brain Metastasis",
    subtitle = paste("Total =", sum(biotype_up$Count)),
    x = "ncRNA Biotype",
    y = "Number of Upregulated ncRNAs"
  ) +
  theme(legend.position = "none")

ggsave("Barplot_Upregulated_ncRNA_Biotypes_BM.png", p_up, width = 8, height = 6, dpi = 600)
print(p_up)




library(ggplot2)
library(dplyr)
library(patchwork)

# =======================================================
# PIE 1: Overall RNA composition
# =======================================================
rna_summary <- results_annot %>%
  mutate(Class = ifelse(gene_biotype == "protein_coding", "Protein-coding", "Non-coding")) %>%
  count(Class, name = "Count") %>%
  mutate(Percent = round(Count / sum(Count) * 100, 1))

p1 <- ggplot(rna_summary, aes(x = "", y = Count, fill = Class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Class, "\n", Percent, "%")), 
            position = position_stack(vjust = 0.5), size = 5, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  theme_void(base_size = 14) +
  labs(title = "Overall RNA Composition")

# =======================================================
# PIE 2: ncRNA Upregulation in BM (using your values)
# =======================================================

noncoding_summary <- data.frame(
  Category = c("Upregulated ncRNAs", "Downregulated ncRNAs"),
  Count = c(634, 72)
) %>%
  mutate(Percent = round(Count / sum(Count) * 100, 1))

p2 <- ggplot(noncoding_summary, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Category, "\n", Count, " (", Percent, "%)")),
            position = position_stack(vjust = 0.5), size = 4.8, color = "black") +
  scale_fill_brewer(palette = "Pastel1") +
  theme_void(base_size = 14) +
  labs(title = "Significant Non-coding RNAs in BM")

# =======================================================
# Combine the two
# =======================================================

combined_plot <- p1 + p2 + 
  plot_annotation(
    title = "RNA Landscape in LUAD Brain Metastasis",
    subtitle = "Distribution of coding vs non-coding RNAs and expression pattern of significant ncRNAs",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave("Combined_RNA_PieCharts_LUAD_BM.png", combined_plot, width = 12, height = 6, dpi = 600)
print(combined_plot)
