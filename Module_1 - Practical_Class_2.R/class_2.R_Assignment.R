# --------------------------
# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
# The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.



# -----------------------------------
# 1. Importing the datasets into R 
# -----------------------------------

# Making sure the folder exists
dir.create("raw_data", showWarnings = FALSE)

# files currently in the project root:
src <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")


# moving them into raw_data (renaming within same drive)
ok <- file.rename(src, file.path("raw_data", basename(src)))

# if any move failed (returns FALSE), fall back to copy+delete
if (any(!ok)) {
  to <- file.path("raw_data", basename(src))
  file.copy(src[!ok], to[!ok], overwrite = TRUE)
  file.remove(src[!ok])
}

# sanity check
list.files("raw_data")
stopifnot(all(file.exists(file.path("raw_data", c("DEGs_Data_1.csv","DEGs_Data_2.csv")))))


# -------------------------------------------------------------
# 2. classify_gene(): deciding gene status using logFC and padj
# -------------------------------------------------------------
classify_gene <- function(logFC, padj) {
  
  if ((padj < 0.05) & (logFC > 1)) {
    return("Upregulated")
  } else if ((padj < 0.05) & (logFC < -1)) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# ----------------------------------------------------
# 3. Classifying gene of two dataset within loop
# ----------------------------------------------------

#  Setting up folders 
input_dir <- "raw_data" 
output_dir <- "results"

# creating output folder if not already exists 
if(!dir.exists(output_dir)){
  dir.create(output_dir)
} 

# ----------------------------------------------------
# 4. Listing files to process 
# ----------------------------------------------------
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv") 


# -------------------------------------------------------------------------
# 5. Preparing empty list with 2 slots, named by file to store results in R 
# -------------------------------------------------------------------------
result_list <- list() 


# ---------------------------------------------------------------------------------
# 6.Applying a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
# ---------------------------------------------------------------------------------
for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path  <- file.path(input_dir, file_names)
  
  # Importing dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  # Handling missing values
  
  if("lofFC" %in% names(data)){
    missing_count <- sum(is.na(data$logFC))
    
    cat("Missing values in 'logFC':", missing_count, "\n")
    data$logFC[is.na(data$logFC)] <- mean(data$logFC, na.rm = TRUE)
  }
  
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- mean(data$padj, na.rm = TRUE)
  }
}


# --------------------------------
# 7. making sure these are numeric
# --------------------------------
data$logFC <- suppressWarnings(as.numeric(data$logFC))
data$padj  <- suppressWarnings(as.numeric(data$padj))


# --------------------------------
# 8. replacing missing padj with 1
# --------------------------------
data$padj[is.na(data$padj)] <- 1


# --------------------------------
# 9. Adding a new column 'status'
# --------------------------------
# element-wise classification (apply the scalar function row by row)
data$status <- mapply(classify_gene, data$logFC, data$padj)


# ----------------------------------------------------
# 10. Storing rsults in R
# ----------------------------------------------------
result_list[[file_names]] <- data

# Accessing results for each file
results_1 <- result_list[["DEGs_Data_1.csv"]] 
results_2 <- result_list[["DEGs_Data_2.csv"]]


# ----------------------------------------------------  
# 11. Printing summary counts of significant, upregulated, and downregulated genes
# ----------------------------------------------------
output_file_path <- file.path(output_dir, paste0("Gene_results", file_names))
write.csv(data, output_file_path, row.names = FALSE)
cat("Results saved to:", output_file_path, "\n")


# ----------------------------------------------------
# 12. Using table() for summaries
# ----------------------------------------------------
print(table(data$status))
cat("Significant (padj < 0.05):", sum(data$padj < 0.05, na.rm = TRUE), "\n")


# Saving the entire R workspace
save.image(file = "Akansha_Naithani_Class_2.R_Assignment.RData")


# Loading workspace
load("Akansha_Naithani_Class_2.R_Assignment.RData")

