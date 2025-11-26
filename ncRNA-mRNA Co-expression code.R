install.packages(c("Hmisc", "dplyr", "pheatmap", "RColorBrewer", "ggplot2"))
library(Hmisc)        # rcorr() for fast correlation + p-values
library(dplyr)        # For dataframe manipulation
library(pheatmap)     # For heatmaps
library(RColorBrewer) # For color palettes
library(ggplot2)      # For additional plots


# ============================================
# PREPARE EXPRESSION MATRICES
# ============================================

# ncRNA matrix (706 significant ncRNAs from Part 1)
ncRNA_expr <- combined_expr_norm[sig_results$hgnc_symbol, ]

# mRNA matrix (102 DEGs - THE CORRECT MATRIX)
mRNA_expr <- mRNA_expr  # Your 102×140 matrix

# Find common samples between both datasets
common_samples <- intersect(colnames(ncRNA_expr), colnames(mRNA_expr))
ncRNA_expr <- ncRNA_expr[, common_samples]
mRNA_expr <- mRNA_expr[, common_samples]

cat("=== DATASET SUMMARY ===\n")
cat("ncRNAs:", nrow(ncRNA_expr), "\n")
cat("mRNAs:", nrow(mRNA_expr), "\n")
cat("Samples:", length(common_samples), "\n")
cat("Total possible ncRNA-mRNA pairs:", nrow(ncRNA_expr) * nrow(mRNA_expr), "\n\n")

# ============================================
# COMBINE MATRICES FOR HMISC
# ============================================

# Combine ncRNAs and mRNAs into one matrix
combined_matrix <- rbind(ncRNA_expr, mRNA_expr)

# Transpose because rcorr() expects variables in columns, samples in rows
combined_matrix_t <- t(combined_matrix)


# Compute correlations and p-values in one step
rc <- rcorr(combined_matrix_t, type = "spearman")

# Extract full correlation matrix
cor_mat_full <- rc$r

# Extract full p-value matrix  
p_mat_full <- rc$P

# ============================================
# EXTRACT ONLY ncRNA–mRNA PAIRS
# ============================================

# Get gene names
ncRNA_ids <- rownames(ncRNA_expr)
mRNA_ids <- rownames(mRNA_expr)

# Extract submatrix: ncRNAs (rows) vs mRNAs (columns)
cor_matrix <- cor_mat_full[ncRNA_ids, mRNA_ids]
pval_matrix <- p_mat_full[ncRNA_ids, mRNA_ids]

cat("Correlation matrix dimensions:", dim(cor_matrix), "\n")

# ============================================
# ADJUST P-VALUES FOR MULTIPLE TESTING
# ============================================

cat("Applying FDR correction...\n")

# Apply Benjamini-Hochberg correction per mRNA (column-wise)
fdr_matrix <- apply(pval_matrix, 2, p.adjust, method = "BH")

# ============================================
# IDENTIFY SIGNIFICANT CORRELATIONS
# ============================================

# Set significance thresholds
cor_threshold <- 0.6    # Minimum correlation strength
fdr_threshold <- 0.05   # Maximum false discovery rate

# Find indices of significant correlations
sig_idx <- which(abs(cor_matrix) > cor_threshold & 
                   fdr_matrix < fdr_threshold, arr.ind = TRUE)

cat("Significant correlations found:", nrow(sig_idx), "\n")

# ============================================
# CREATE RESULTS DATAFRAME
# ============================================

if(nrow(sig_idx) > 0) {
  
  cor_results <- data.frame(
    ncRNA = rownames(cor_matrix)[sig_idx[, 1]],
    mRNA = colnames(cor_matrix)[sig_idx[, 2]],
    Correlation = cor_matrix[sig_idx],
    Pvalue = pval_matrix[sig_idx],
    FDR = fdr_matrix[sig_idx],
    Regulation = ifelse(cor_matrix[sig_idx] > 0, "Positive", "Negative")
  ) %>% arrange(desc(abs(Correlation)))
  
  cat("=== RESULTS SUMMARY ===\n")
  cat("Total significant correlations:", nrow(cor_results), "\n")
  cat("Positive correlations:", sum(cor_results$Regulation == "Positive"), "\n")
  cat("Negative correlations:", sum(cor_results$Regulation == "Negative"), "\n")
  cat("Average |correlation|:", round(mean(abs(cor_results$Correlation)), 3), "\n\n")
  
  # Show top 10 strongest correlations
  print(head(cor_results, 10))
  
}

# ============================================
# FOCUS ON HUB GENES
# ============================================

if(exists("cor_results") && nrow(cor_results) > 0) {
  
  # Your 20 hub genes
  hub_genes <- c(
    "CD69", "CCR7", "CD27", "CD2", "CCL5", "CD247",
    "GZMA", "CD3D", "GZMK", "IL2RB", "CXCL9", "CCL19",
    "CXCL13", "CXCL10", "CD48", "VCAM1", "CD79B", "SLAMF6", "CD79A", "SH2D1A"
  )
  
  hub_correlations <- cor_results %>% 
    filter(mRNA %in% hub_genes) %>%
    arrange(mRNA, desc(abs(Correlation)))
  
  cat("Hub gene correlations:", nrow(hub_correlations), "\n")
  
  # Hub gene summary
  hub_summary <- hub_correlations %>%
    group_by(mRNA) %>%
    summarise(
      n_correlations = n(),
      avg_correlation = mean(abs(Correlation)),
      max_correlation = max(abs(Correlation))
    ) %>%
    arrange(desc(n_correlations))
  
  print(hub_summary)
}

# ============================================
# HEATMAP VISUALIZATIONS
# ============================================

if(exists("cor_results") && nrow(cor_results) > 0) {
  cat("Creating heatmaps...\n")
  
  # HEATMAP 1: Top Correlations - WITH DATA CLEANING
  top_pairs <- head(cor_results, 50)
  
  # Create a proper matrix without warnings
  unique_ncRNAs <- unique(top_pairs$ncRNA)
  unique_mRNAs <- unique(top_pairs$mRNA)
  
  if(length(unique_ncRNAs) > 1 && length(unique_mRNAs) > 1) {
    # Initialize matrix with NAs
    heatmap_data <- matrix(NA, 
                           nrow = length(unique_ncRNAs),
                           ncol = length(unique_mRNAs),
                           dimnames = list(unique_ncRNAs, unique_mRNAs))
    
    # Fill matrix with correlation values
    for(i in 1:nrow(top_pairs)) {
      heatmap_data[top_pairs$ncRNA[i], top_pairs$mRNA[i]] <- top_pairs$Correlation[i]
    }
    
    # Remove rows/columns with all NAs AND clean data
    heatmap_data <- heatmap_data[rowSums(!is.na(heatmap_data)) > 0, 
                                 colSums(!is.na(heatmap_data)) > 0, drop = FALSE]
    
    # CRITICAL FIX: Remove any NA/NaN/Inf values that break hclust
    heatmap_data_clean <- heatmap_data
    heatmap_data_clean[is.na(heatmap_data_clean) | is.nan(heatmap_data_clean) | is.infinite(heatmap_data_clean)] <- 0
    
    if(nrow(heatmap_data_clean) > 1 && ncol(heatmap_data_clean) > 1) {
      pheatmap(heatmap_data_clean,
               color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               main = "Top 50 Strongest ncRNA-mRNA Correlations",
               na_col = "white",
               filename = "heatmap_top_correlations.png")
      cat("✓ Heatmap 1 saved: heatmap_top_correlations.png\n")
    }
  }
  
  # HEATMAP 2: Hub Genes Focus - WITH DATA CLEANING
  if(exists("hub_correlations") && nrow(hub_correlations) > 0) {
    hub_ncRNAs <- unique(hub_correlations$ncRNA)
    hub_mRNAs <- unique(hub_correlations$mRNA)
    
    if(length(hub_ncRNAs) > 1 && length(hub_mRNAs) > 1) {
      hub_matrix <- cor_matrix[hub_ncRNAs, hub_mRNAs, drop = FALSE]
      hub_matrix[abs(hub_matrix) < cor_threshold] <- 0
      hub_matrix <- hub_matrix[rowSums(abs(hub_matrix)) > 0,
                               colSums(abs(hub_matrix)) > 0, drop = FALSE]
      
      # CRITICAL FIX: Clean the data
      hub_matrix_clean <- hub_matrix
      hub_matrix_clean[is.na(hub_matrix_clean) | is.nan(hub_matrix_clean) | is.infinite(hub_matrix_clean)] <- 0
      
      if(nrow(hub_matrix_clean) > 1 && ncol(hub_matrix_clean) > 1) {
        pheatmap(hub_matrix_clean,
                 color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                 main = "ncRNA Correlations with Hub Genes",
                 filename = "heatmap_hub_genes.png")
        cat("✓ Heatmap 2 saved: heatmap_hub_genes.png\n")
      }
    }
  }
  
  # HEATMAP 3: Immune Hub Genes - WITH DATA CLEANING
  immune_genes <- c("CD69", "GZMA", "CCL5", "CXCL9", "CXCL10", "CCR7", "CD3D")
  available_immune <- intersect(colnames(cor_matrix), immune_genes)
  
  if(length(available_immune) > 0) {
    immune_matrix <- cor_matrix[, available_immune, drop = FALSE]
    immune_matrix[abs(immune_matrix) < cor_threshold] <- 0
    immune_matrix <- immune_matrix[rowSums(abs(immune_matrix)) > 0, , drop = FALSE]
    
    # CRITICAL FIX: Clean the data
    immune_matrix_clean <- immune_matrix
    immune_matrix_clean[is.na(immune_matrix_clean) | is.nan(immune_matrix_clean) | is.infinite(immune_matrix_clean)] <- 0
    
    if(nrow(immune_matrix_clean) > 1 && ncol(immune_matrix_clean) > 1) {
      pheatmap(immune_matrix_clean,
               color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               main = "ncRNA Correlations with Immune Hub Genes",
               filename = "heatmap_immune_hub_genes.png")
      cat("✓ Heatmap 3 saved: heatmap_immune_hub_genes.png\n")
    } else {
      cat("Heatmap 3: Not enough data after cleaning (", 
          nrow(immune_matrix_clean), "x", ncol(immune_matrix_clean), ")\n")
    }
  }
}

# ============================================
# BETTER VISUALIZATION FOR 1×7 DATA
# ============================================

if(length(available_immune) > 0 && nrow(immune_matrix) == 1) {
  cat("Creating specialized visualization for single ncRNA...\n")
  
  # Extract the data
  single_ncRNA_name <- rownames(immune_matrix)[1]
  correlation_data <- data.frame(
    Immune_Gene = colnames(immune_matrix),
    Correlation = as.numeric(immune_matrix[1, ]),
    Significant = abs(as.numeric(immune_matrix[1, ])) > 0
  ) %>% arrange(desc(abs(Correlation)))
  
  # Remove zero correlations
  correlation_data <- correlation_data %>% filter(Significant == TRUE)
  
  cat("Significant correlations for", single_ncRNA_name, ":\n")
  print(correlation_data)
  
  if(nrow(correlation_data) > 0) {
    # Create a lollipop plot
    p <- ggplot(correlation_data, aes(x = reorder(Immune_Gene, abs(Correlation)), 
                                      y = abs(Correlation), 
                                      color = Correlation > 0)) +
      geom_segment(aes(xend = Immune_Gene, yend = 0), linewidth = 1) +
      geom_point(size = 4) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                         labels = c("Positive", "Negative"),
                         name = "Regulation") +
      coord_flip() +
      theme_minimal() +
      labs(title = paste("Single Master ncRNA:", single_ncRNA_name),
           subtitle = "Correlations with Immune Hub Genes",
           x = "Immune Gene",
           y = "|Correlation Coefficient|") +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(color = "gray40"))
    
    print(p)
    ggsave("single_ncRNA_immune_correlations.png", p, width = 10, height = 6)
    cat("✓ Single ncRNA visualization saved: single_ncRNA_immune_correlations.png\n")
    
    # Also create a simple table
    write.csv(correlation_data, "single_ncRNA_immune_correlations.csv", row.names = FALSE)
    cat("✓ Data saved: single_ncRNA_immune_correlations.csv\n")
  }
}



# ============================================
# TRY RELAXED THRESHOLDS FOR IMMUNE GENES
# ============================================

cat("=== TRYING RELAXED THRESHOLDS ===\n")

immune_relaxed_threshold <- 0.5  # Lower from 0.6
fdr_relaxed_threshold <- 0.1     # Relax from 0.05

immune_matrix_relaxed <- cor_matrix[, available_immune, drop = FALSE]

# Apply relaxed significance filter
for(i in 1:length(available_immune)) {
  gene <- available_immune[i]
  immune_matrix_relaxed[abs(immune_matrix_relaxed[, i]) < immune_relaxed_threshold | 
                          fdr_matrix[, gene] >= fdr_relaxed_threshold, i] <- 0
}

immune_matrix_relaxed <- immune_matrix_relaxed[rowSums(abs(immune_matrix_relaxed)) > 0, , drop = FALSE]
cat("With relaxed thresholds (|r| > 0.5):", dim(immune_matrix_relaxed), "\n")

if(nrow(immune_matrix_relaxed) > 1) {
  pheatmap(immune_matrix_relaxed,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = "Immune Hub Genes (Relaxed: |r| > 0.5, FDR < 0.1)",
           filename = "heatmap_immune_relaxed.png")
  cat("✓ Relaxed threshold heatmap saved: heatmap_immune_relaxed.png\n")
}
# ============================================
# SAVE RESULTS
# ============================================

if(exists("cor_results") && nrow(cor_results) > 0) {
  write.csv(cor_results, "ncRNA_mRNA_coexpression_results.csv", row.names = FALSE)
  cat("✓ Full results saved: ncRNA_mRNA_coexpression_results.csv\n")
}

if(exists("hub_correlations") && nrow(hub_correlations) > 0) {
  write.csv(hub_correlations, "hub_gene_correlations.csv", row.names = FALSE)
  cat("✓ Hub gene results saved: hub_gene_correlations.csv\n")
}


