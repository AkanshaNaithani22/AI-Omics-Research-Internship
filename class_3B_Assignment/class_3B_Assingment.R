# =====================================================================
# AI & Omics Research Internship (Module II • Class 3B) — AE edition
# Study: E-MTAB-11191  (ArrayExpress, Affymetrix)
# Workflow: QC (raw) → RMA normalization → QC (norm) → filtering → groups
# =====================================================================


#### 0) Packages ####
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c(
  "ArrayExpress","affy","arrayQualityMetrics","limma",
  "Biobase","AnnotationDbi","affyPLM", "hgu133plus2cdf"
), ask = FALSE, update = TRUE)

install.packages("matrixStats", dependencies = TRUE)

library(ArrayExpress)
library(affy)
library(hgu133plus2cdf)
library(arrayQualityMetrics)
library(limma)
library(Biobase)
library(AnnotationDbi)
library(matrixStats)


#### 1) Paths ####
data_dir    <- "raw_data"     # where raw files will live
results_dir <- "results"      # where outputs go
dir.create(data_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


#### 2) Downloading AE study (RAW + MAGE-TAB) ####

ae_raw  <- getAE("E-MTAB-11191", type = "raw",  path = data_dir, extract = TRUE)
ae_full <- getAE("E-MTAB-11191", type = "full", path = data_dir, extract = TRUE) # to get SDRF/IDF/ADF


#### 3) Locating CEL files + SDRF (phenodata) ####
cel_files <- list.files(data_dir, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE,
                        recursive = TRUE, full.names = TRUE)
stopifnot(length(cel_files) > 0)

sdrf_path <- ae_full$sdrf
if (is.null(sdrf_path) || !file.exists(sdrf_path)) {
  sdrf_path <- list.files(data_dir, pattern = "\\.sdrf(\\.txt)?$", 
                          recursive = TRUE, full.names = TRUE)[1]
}
stopifnot(!is.null(sdrf_path), file.exists(sdrf_path))

targets <- read.delim(sdrf_path, check.names = FALSE)

# Column with CEL filenames detection:

cel_col <- grep("array.*data.*file|hyb.*file|scan.*file|cel",
                names(targets), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(cel_col))


# Aligning SDRF rows to the CEL list

cf <- tolower(basename(cel_files))
tf <- tolower(basename(gsub("\\\\","/", targets[[cel_col]])))
tf <- sub("\\.gz$", "", tf)
tf <- ifelse(grepl("\\.cel$", tf), tf, paste0(tf, ".cel"))
idx <- match(cf, tf)
if (any(is.na(idx))) {
  cat("Unmatched CELs (disk not in SDRF):\n"); print(setdiff(cf, tf))
  cat("Unmatched SDRF entries (SDRF not on disk):\n"); print(setdiff(tf, cf))
  stop("Fix filename mismatches and rerun.")
}
targets <- targets[idx, , drop = FALSE]



#### 4) Reading raw arrays and running QC (pre-normalization) ####

raw_data <- ReadAffy(filenames = cel_files)  # AffyBatch

# (platform/annotation check)

raw_data

qc_raw_dir <- file.path(results_dir, "QC_Raw")

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "results/qc_raw_dir",
                    force = TRUE,
                    do.logtransform = TRUE)


# flagged arrays” count programmatically:

flag_file <- list.files(qc_raw_dir, pattern = "flagged.*csv$", full.names = TRUE)
if (length(flag_file)) {
  flagged_raw <- read.csv(flag_file[1])
  cat("Arrays flagged in RAW QC:", nrow(flagged_raw), "\n")
}


#### 5) RMA normalization + QC (post-normalization) ####

normalized_data <- rma(raw_data)  # ExpressionSet

qc_norm_dir <- file.path(results_dir, "QC_Normalized")

arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "results/qc_norm_dir",
                    force = TRUE)


flag_file2 <- list.files(qc_norm_dir, pattern = "flagged.*csv$", full.names = TRUE)
if (length(flag_file2)) {
  flagged_norm <- read.csv(flag_file2[1])
  cat("Arrays flagged in NORMALIZED QC:", nrow(flagged_norm), "\n")
}



#### 6) Extractung normalized expression matrix ####

expr <- as.data.frame(exprs(normalized_data))  # probes x samples
write.csv(expr, file.path(results_dir, "E-MTAB-11191_RMA_expression.csv"))


#### 7) “Soft” intensity filtering (like class) ####
# Compute per-probe medians across samples

row_median <- matrixStats::rowMedians(as.matrix(expr))

# Visualizing distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Choosing a threshold
threshold <- 3.5
keep <- row_median > threshold
filtered_expr <- expr[keep, ]
cat("Probes before filtering:", nrow(expr), 
    " | after filtering:", nrow(filtered_expr), "\n")

write.csv(filtered_expr, file.path(results_dir, "E-MTAB-11191_RMA_filtered.csv"))


#### 8) Preparing phenotype (groups) ####

# Inspecting SDRF to pick a column describing condition/groups:

View(names(targets)); head(targets[, grep("disease|group|condition|phenotype|case|control|status|characteristics", names(targets), ignore.case=TRUE), drop=FALSE])

# auto detecting sensible column:

grp_col <- grep("disease|group|condition|phenotype|status|case|control",
                names(targets), ignore.case = TRUE, value = TRUE)[1]
if (is.na(grp_col)) grp_col <- names(targets)[1]  # fallback: first column

groups_raw <- as.character(targets[[grp_col]])

# Cleaning labels: fold to lower-case and simplify common MI/control names
clean <- function(x) {
  x <- tolower(trimws(x))
  x <- sub("^characteristics\\[.*\\]:\\s*", "", x)
  x
}
groups <- factor(vapply(groups_raw, clean, ""))
table(groups)


# Relabeling to exactly two classes

levels(groups) <- c("normal","systemic_lupus_erythematosus")  # adjusting to SDRF labels


# Saving phenotype with chosen group column

phen_out <- cbind(sample = colnames(expr), group = groups, targets)
write.csv(phen_out, file.path(results_dir, "E-MTAB-11191_pheno_groups.csv"), row.names = FALSE)


cat("\nSummary for submission:\n")
cat(" RAW QC flagged:", if (exists("flagged_raw")) nrow(flagged_raw) else NA, "\n")
cat(" NORM QC flagged:", if (exists("flagged_norm")) nrow(flagged_norm) else NA, "\n")
cat(" Probes before:", nrow(expr), " after:", nrow(filtered_expr), "\n")
print(table(groups))




