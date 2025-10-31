# =====================================================================
# AI & Omics Research Internship (Module II • Class 3B) — AE edition
# ---------------------------------------------------------------------
# Study: E-MTAB-11191  (ArrayExpress, Affymetrix)
# ---------------------------------------------------------------------
# Workflow: QC (raw) → RMA normalization → QC (norm) → filtering → groups
# =====================================================================


#######################################################################
#### 0. Installing and Loading Required Packages ####
#######################################################################


# Bioconductor - R packages for omics data analysis
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Installing Bioconductor packages
BiocManager::install(c(
  "ArrayExpress","affy","limma",
  "Biobase","AnnotationDbi","affyPLM", "hgu133plus2cdf"
), ask = FALSE, update = TRUE, force = TRUE)

BiocManager::install(("arrayQualityMetrics"), ask = FALSE, update = TRUE, force = TRUE)

# Installing CRAN packages for data manipulation
install.packages("matrixStats", dependencies = TRUE)

install.packages("dplyr")

# Loading Required Libraries
library(BiocManager)
library(ArrayExpress)
library(hgu133plus2cdf)
library(limma)
library(Biobase)
library(AnnotationDbi)
library(matrixStats)
library(dplyr)
library(affy)
library(arrayQualityMetrics)


# -------------------------------------
#### Series Matrix Files and convert to ExpressionSet ####
# -------------------------------------

# Unzipped raw data and SDRF files
unzip("E-MTAB-11191/Raw_data_files.zip", exdir = "E-MTAB-11191/Raw_data_files")

unzip("E-MTAB-11191/SDRF_file.zip", exdir = "E-MTAB-11191/SDRF_files")

# Reading CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "E-MTAB-11191/Raw_data_files")

raw_data

# ---------------------------------------------------
#### Quality Control (QC) Before Pre-processing ####
# ---------------------------------------------------

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_data",
                    force = TRUE,
                    do.logtransform = TRUE)


# -------------------------------------------------------
#### RMA (Robust Multi-array Average) Normalization ####
# -------------------------------------------------------

normalized_data <- rma(raw_data)

expr_mat <- exprs(normalized_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data[,c(2,9,34,48,50,52,53,57,61,65,67,92,97,101)],
                    outdir = "Results/QC_Normalized_data",
                    force = TRUE,
                    reporttitle     = "QC Normalized (selected arrays)")


# Extracting normalized expression values into a data frame
processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)


# ---------------------------------------------------------------------------
#### Filter Low-Variance Transcripts (“soft” intensity based filtering) ####
# ---------------------------------------------------------------------------

# Calculating median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))

row_median

# Visualizing distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")


# Setting a threshold to remove low variance probes 
threshold <- 3.5

abline(v = threshold, col = "black", lwd = 2)

# Selecting probes above threshold
indx <- row_median > threshold

indx

filtered_data <- processed_data[indx, ]

# Reading SDRF and align to columns (PHENOTYPE DATA)
# Finding the SDRF file that was unzipped to E-MTAB-11191/SDRF_files/

sdrf_file <- list.files("E-MTAB-11191/SDRF_files",
                        pattern = "\\.sdrf\\.txt(\\.gz)?$",
                        ignore.case = TRUE, recursive = TRUE, full.names = TRUE)[1]
stopifnot(length(sdrf_file) == 1)

pheno_all <- read.delim(sdrf_file, check.names = FALSE)

# SDRF column that contains the array file names (matching to CEL basenames)
map_cols <- c("Array Data File", "ArrayDataFile", "Array Data File Name")
key <- map_cols[map_cols %in% names(pheno_all)][1]
stopifnot(!is.null(key))

m <- match(colnames(expr_mat), basename(pheno_all[[key]]))
stopifnot(!any(is.na(m)))          

phenotype_data <- pheno_all[m, , drop = FALSE]
rownames(phenotype_data) <- colnames(expr_mat)

probe_ids <- rownames(expr_mat)

feat_annot <- AnnotationDbi::select(hgu133plus2.db,
                                    keys = probe_ids,
                                    keytype = "PROBEID",
                                    columns = c("SYMBOL","GENENAME","ENTREZID"))

# making a tidy feature table aligned to expr_mat rows
feature_data <- feat_annot[match(probe_ids, feat_annot$PROBEID), ]
rownames(feature_data) <- probe_ids

dim(expr_mat)                  # probes x samples (expression)
dim(phenotype_data)            # samples x phenotype-fields
dim(feature_data)              # probes x annotation-fields
normalized_data

# Renaming filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)

# Overwriting processed data with filtered dataset
processed_data <- filtered_data

# -----------------------------------
#### Phenotype Data Preparation ####
# -----------------------------------

class(phenotype_data$`Characteristics[disease]`)

phenotype_data$`Characteristics[disease]`

# Defining experimental groups (normal vs disease)
groups <- factor(phenotype_data$`Characteristics[disease]`,
                 levels = c("normal", "systemic lupus erythematosus"),
                 labels = c("normal", "disease"))


class(groups)
levels(groups)
