

#################################################################
#              ENHANCED PIPELINE: IMMUNE + WGCNA + HUB ncRNAs      
#################################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

  library(WGCNA) #
  library(Hmisc)
  library(dplyr)
  library(tibble) 
  library(tidyr) 
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer) 
  library(reshape2) 
  library(corrplot)


options(stringsAsFactors = FALSE)
allowWGCNAThreads()

#################################################################
### 1. PREPARE ncRNA EXPRESSION MATRIX (samples × ncRNAs)
#################################################################

sig_ncRNAs <- sig_results$hgnc_symbol
sig_ncRNAs <- intersect(sig_ncRNAs, rownames(master_norm_matrix))

ncRNA_matrix <- master_norm_matrix[sig_ncRNAs, ]
matrix_B_ncRNA <- t(ncRNA_matrix)          # transpose for WGCNA (samples × genes)

cat("ncRNA matrix size:", dim(matrix_B_ncRNA), "\n")
cat("Number of significant ncRNAs:", ncol(matrix_B_ncRNA), "\n")

#################################################################
### 2. LOAD + FILTER IMMUNE ESTIMATION MATRIX (CIBERSORTx)
#################################################################

immune_raw <- read.delim("CIBERSORTx_Job2_Results.txt", sep="\t")

immune_filtered <- subset(immune_raw, P.value < 0.05)
rownames(immune_filtered) <- immune_filtered$Mixture

matrix_A_immune <- immune_filtered[, 2:23]

cat("Immune matrix size:", dim(matrix_A_immune), "\n")
cat("Immune cell types:", paste(colnames(matrix_A_immune), collapse=", "), "\n")

#################################################################
### 3. ALIGN SAMPLES BETWEEN ncRNA MATRIX & IMMUNE MATRIX
#################################################################

common_samples <- intersect(rownames(matrix_A_immune), rownames(matrix_B_ncRNA))

matrix_A_immune <- matrix_A_immune[common_samples, ]
matrix_B_ncRNA  <- matrix_B_ncRNA[common_samples, ]

cat("Final shared samples:", length(common_samples), "\n")

# Create trait data for group information
traitData <- data.frame(
  BrainMetastasis = as.numeric(group_all[match(common_samples, colnames(master_norm_matrix))] == "BM"),
  row.names = common_samples
)

#################################################################
### 4. WGCNA DATA PREPARATION
#################################################################

datExpr <- as.data.frame(matrix_B_ncRNA)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("Removed", sum(!gsg$goodGenes), "genes with missing data\n")
}

sampleTree <- hclust(dist(datExpr), method = "average")

pdf("WGCNA_SampleClustering.pdf", width=12, height=6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

#################################################################
### 5. CHOOSE SOFT-THRESHOLD
#################################################################

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5, networkType="unsigned")

softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 6
  cat("Automatic power selection failed. Using power =", softPower, "\n")
} else {
  cat("Selected soft threshold:", softPower, "\n")
}

pdf("WGCNA_SoftThreshold.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.90, col="blue")

plot(sft$fitIndices[,1], sft$fitIndices[,5], type="n",
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
     main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

#################################################################
### 6. CONSTRUCT NETWORK + IDENTIFY MODULES
#################################################################

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
moduleLabels <- net$colors

module_summary <- as.data.frame(table(moduleColors))
colnames(module_summary) <- c("Module", "Number_of_ncRNAs")
print(module_summary)

pdf("WGCNA_Dendrogram_Modules.pdf", width=15, height=8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#################################################################
### 7. MODULE EIGENGENES
#################################################################

MEs <- net$MEs
MEs <- orderMEs(MEs)

# Save module eigengenes
write.csv(MEs, "WGCNA_Module_Eigengenes.csv")

#################################################################
### 8. MODULE–IMMUNE CELL CORRELATION
#################################################################

moduleTraitCor <- cor(MEs, matrix_A_immune, method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Create annotated heatmap
pdf("WGCNA_Module_Immune_Heatmap.pdf", width=12, height=10)

# Create text matrix for p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(8, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(matrix_A_immune),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = "Module–Immune Cell Correlations")

dev.off()

# Also create pheatmap version
pheatmap(moduleTraitCor,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(50),
         main = "Module–Immune Cell Correlations (Spearman)",
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         filename = "WGCNA_Module_Immune_Heatmap_pheatmap.pdf")

write.csv(moduleTraitCor, "Module_Immune_Correlation.csv")
write.csv(moduleTraitPvalue, "Module_Immune_Correlation_Pvalues.csv")

#################################################################
#  GENE IDENTIFICATION
#################################################################

cat("=== SIMPLE HUB GENE IDENTIFICATION ===\n")

# Focus only on the turquoise module (your main module)
turquoise_genes <- colnames(datExpr)[moduleColors == "turquoise"]

cat("Turquoise module has", length(turquoise_genes), "genes\n")

if(length(turquoise_genes) > 0) {
  # Calculate module membership for turquoise genes
  MM_turquoise <- as.data.frame(cor(datExpr[, turquoise_genes], MEs[, "ME1"], use = "p"))
  colnames(MM_turquoise) <- "ModuleMembership"
  
  # Calculate average immune correlation for each gene
  immune_scores <- sapply(turquoise_genes, function(gene) {
    cor_values <- sapply(colnames(matrix_A_immune), function(immune_cell) {
      abs(cor(datExpr[, gene], matrix_A_immune[, immune_cell], use = "p", method = "spearman"))
    })
    mean(cor_values, na.rm = TRUE)
  })
  
  # Create results
  hub_results <- data.frame(
    ncRNA = turquoise_genes,
    Module = "turquoise",
    ModuleMembership = abs(MM_turquoise$ModuleMembership),
    Average_Immune_Correlation = immune_scores,
    Hub_Score = abs(MM_turquoise$ModuleMembership) * immune_scores
  ) %>% arrange(desc(Hub_Score))
  
  # Save results
  write.csv(hub_results, "Turquoise_Module_Hub_ncRNAs.csv", row.names = FALSE)
  
  cat("Top 10 hub ncRNAs in turquoise module:\n")
  print(head(hub_results, 10))
  
  cat("\n✓ Hub analysis complete! Saved: Turquoise_Module_Hub_ncRNAs.csv\n")
}

#################################################################
### 10. SUMMARY STATISTICS AND FINAL OUTPUT
#################################################################

cat("\n=== WGCNA-IMMUNE PIPELINE COMPLETE ===\n")

# Count modules (excluding grey)
valid_modules <- unique(moduleColors[moduleColors != "grey"])
cat("Modules identified:", length(valid_modules), "\n")
cat("Module names:", paste(valid_modules, collapse = ", "), "\n\n")

# Hub ncRNA summary
cat("Total ncRNAs in turquoise module:", length(turquoise_genes), "\n")
cat("Hub ncRNAs identified:", nrow(hub_results), "\n")

# Show top 10 hubs again
cat("\nTop 10 Hub ncRNAs:\n")
print(head(hub_results[, c("ncRNA", "ModuleMembership", 
                           "Average_Immune_Correlation", "Hub_Score")], 10))

# Immune–module correlation summary
module_immune_summary <- data.frame(
  Module = "turquoise",
  Most_Correlated_Immune_Cell = colnames(matrix_A_immune)[which.max(abs(moduleTraitCor["ME1", ]))],
  Max_Correlation = max(abs(moduleTraitCor["ME1", ])),
  Number_of_ncRNAs = length(turquoise_genes)
)

write.csv(module_immune_summary, "Module_Immune_Summary.csv", row.names = FALSE)

cat("\nModule–Immune Association Summary:\n")
print(module_immune_summary)

cat("\n✓ SIMPLE WGCNA + IMMUNE + HUB PIPELINE COMPLETE!\n")
cat("✓ All hub ncRNAs saved to Turquoise_Module_Hub_ncRNAs.csv\n")

save.image(file = "final.RData")

