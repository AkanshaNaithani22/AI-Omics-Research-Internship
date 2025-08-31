# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Module I: Getting Started with R
# -------------------------------------------------------------------
# ===================================================================

# --------------------
# Topics in this module
# --------------------
#   1. Operators in R
#   2. Data Structures in R
#   3. User-defined Functions in R
#   4. Automating Workflows with for-Loops
# --------------------------------------------------------------------------------------------------

#### 1. Operators in R ####

# Operators are special symbols in R.
# They instruct R to perform actions on values or variables.
# You will use them for assignment, calculations, comparisons, and logical tests.

# Types of operators: 
#   1. Assignment operators
#   2. Arithmetic operators
#   3. Comparison (relational) operators
#   4. Logical operators  

# ----------------------
# Assignment Operators
# ----------------------
# Used to store values inside variables

# <-   (Most common, assigns value on the right to the name on the left)
height <- c(1.75, 1.76, 1.82, 1.67) #height is a vector which is a varibale name that has been assigned those values 

# ->   (Same as above, but leftward assignemnt operator) 
#       first we need to assign the values that we need to store in the particular variable 
c(68, 78, 85, 75) -> weight 

# =     (Also assigns, used in function arguments)
#        works similar to "<-" operator 
smoking_status = c("Yes", "No", "No", "Yes")


# ----------------------
# Arithmetic Operators
# ----------------------
# Perform basic math: +, -, *, /, ^
#   +  addition
#   -  subtraction
#   *  multiplication
#   /  division
#   ^  exponent (to the power of)

# Example: calculate BMI (Body Mass Index = weight / height^2)

BMI <- weight/(height^2) 
BMI

# Note: R applies operations element-wise when variables are vectors.
# This is called "vectorization", every weight is divided by every height squared.

#vactorization 
# R perform operation to every value in the vector

# ---------------------
# Comparison Operators
# ---------------------
# Comparison operators ask logical questions about values.
# They return output as TRUE/FALSE
# They don't calculate a number, they return a logical output.
# They compare values 

#   > greater than
BMI > 25 
#   < less than 
BMI < 18.5
#   >= greater than or equal to 
height >= 1.75
#   <= less than or euqal to
weight <= 65
#   == equal to
smoking_status == "No"
#   != not equal to 
smoking_status != "No"


# In R 
# Yes = TRUE
# No = FALSE


# -------------------
# Logical Operators
# -------------------
# Combine multiple conditions using:
#   &   AND         (both must be TRUE)

# Example: is the patient overweight AND a smoker?
# BMI cutoff = 25, above 25 is overweight
(BMI > 25) & (smoking_status == "Yes") 
# if use single = it means assigning a value 
# using two == means to compare values and obtain exact answers
# output shows in the data only one pateint has above 25 BMI who is also a smoker 

(BMI < 25) & (smoking_status == "Yes")
# output shows only one patient has below BMI plus a smoker 


#   |   OR          (at least one must be TRUE)
# Example: is the patient overweight OR smoker?
(BMI > 25 ) | (smoking_status == "Yes")
BMI
smoking_status

#   !   NOT         (reverse the condition)
# Example: is the patient NOT a smoker?
!(smoking_status == "No")
# condition = yes 
# output = FALSE 

# ----------------------
# Summary of Operators
# ----------------------
#   - Assignment stores values
#   - Arithmetic performs calculations
#   - Comparison returns TRUE/FALSE
#   - Logical combines conditions

# --------------------------------------------------------------------------------------------------

# ------------------------
# 2. Data Structures in R
# ------------------------
# Data structures are how R organizes and store information.

# We commonly use:
#   1. Vectors
#   2. Lists
#   3. Matrices
#   4. Data Frames

# ---------
# Vectors
# ---------
# Simplest data structure in R. 
# Stores a sequence of values of the SAME type

#     - Numeric vector 
num_vec <- c(1,2,3,4)
class(num_vec)

# numeric vector used to perform mathematical calculation

#     - Character vector 
chrc_vector <- c("gene1", "gene2", "gene3")
class(chrc_vector)

#     - Logical vector 
logical_vector <- c(TRUE, FALSE, TRUE)
class(logical_vector)

# Important:
# Maybe R does not throw an error if you combine mixed types.
# Instead, it coerces all values into a single type.
mix_vector <- c("gene1", 1, "gene2", 2)
mean(mix_vector) #all values become character
class(mix_vector)

# Indexing: extract elements with []
num_vec[2] # second element 
num_vec[2:4] # extract value from 2 to 4, : indicates sequence

# you can only combine vectors of equal sequence

# We can treat vectors as columns or as rows

# by column
vec_col <- cbind(num_vec, chrc_vector) 
#in the environment "data", gene1 was repeated in vec_col because of unequal lenghts of both vector values

# when we combine unequal lengths of the vector it will duplicate the value and maybe throw an error 
# when we combine vector in form of list all vector remain in there unchnaged/original state



# ---------
# Lists
# ---------
# Unlike vectors, lists can store different types together.
# numbers, texts, logical even other data frames. 

all_vectors <- list(num_vec, chrc_vector, logical_vector)

# save your raw_data
# save processed data 
# results 

# we access elements with [[]]
all_vectors[[2]]



# ---------
# Matrices
# ---------
# A matrix is a 2D structure (rows * columns)
# A 2D table where all values must be the same type (usually numeric).

#Example: gene expression matrix where rows are genes and
# columns are samples.

my_matrix <- matrix(1:9, nrow = 3, ncol = 3)
my_matrix 

# By default R fills the matrix column wise 
# we chnage using byrow = TRUE
my_matrix <- matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE)


# we access elements with [row, columns]
my_matrix[2,3]
my_matrix[2,] # R takes values from all columns because there is no mention of column number
my_matrix[#rows, #columns]
  
# ---------
# Data Frames
# ---------
# Most important structure for biological datasets.
# Columns can have different types (numeric, character, factor).
  
data <- data.frame(
    patient_id =c("P1", "P2", "P3"),
    age = c(65, 78, NA), # NA is not available values
    diagnosis = c("cancer", "diabetes", "cancer")
)
  
print(data)
  
  
# --------------------
# Dataset Assessment
# --------------------
  
# Functions like str(), head(), dim(), and names() help us explore
# the dataset before analysis
  
# Before analyzing, inspect your dataset to understand its structure.
str(data)       # structure (column names, data types, preview)
head(data)      # first 6 rows
head(data, n=2) # first 2 rows
tail(data)      # last 6 rows
tail(data, 2)   # last 2 rows
dim(data)       # number of rows and columns
names(data)     # column names
  
# Data frame are indexed like matrix with more flexibility 
# access a column directly  
  
# Access data using:
data$patient_id     # extract single column
data [c(1,3), c(2,3)]
  
data[1:2, c(1,3)]   # extract specific rows and columns
  
# Create new columns:
data$new_column <- c(1, 2, 3)
  
  
  
# ----------------
# Missing Values
# ----------------
# Real data often contains missing values (NA).
# You must check and handle them before analysis.
  
is.na(data)                # identify missing values
sum(is.na(data))           # total missing values
colSums(is.na(data))       # missing values per column
rowSums(is.na(data))       # missing values per row
  
  
# Ways to handle NA:
  
# remove rows with NA
clean_data_1 <- na.omit(data)   
clean_data_1
  
# remove columns with NA
clean_data_2 <- data[,  colSums(is.na(data))==0]
clean_data_2
  
# replace NA value with 0
clean_data_3 <- data
clean_data_3[is.na(clean_data_3)] <- 0
clean_data_3
  
# replace NA value with mean
clean_data_4 <- data
clean_data_4[is.na(clean_data_4)] <- mean(data$age, na.rm = TRUE)
clean_data_4
  
  
  
# ------------------------------
# Summary of Data Structures:
# ------------------------------
  
#   - Vectors: simple sequences of same data type
#   - Lists: mix of different data types
#   - Matrices: numeric tables
#   - Data Frames: mixed-type tables 
  
  
#--------------------
# 3. Functions in R
#--------------------
# Functions let us wrap code into reusable blocks.
  
# function is  reusable block of code 
# Why use functions?
#   - Avoid repetition
#   - Organize and simplify code
#   - Reuse across projects (save it for later use)
#   - Share with others
  
# A function in R has 4 key parts:
#   1. Name         -> the name you give to the function
#   2. Arguments    -> the inputs you provide to the function
#   3. Body         -> the set of operations the function performs
#   4. Return Value -> the output the function gives back

# Example: A function to calculate Body Mass Index (BMI)
  
# 1. Function Name: calculate_BMI
# 2. Arguments: function(x) e.g  # weight (in kg), height (in meters)
# 3. Body: performs BMI calculation e.g   # Formula: BMI = {weight / (height^2)}
# 4. Return Value: the BMI value # return(BMI)
  
calculate_BMI <- function(weight, height) {
    #operation we want to perform
    bmi <- weight/(height^2)
    
    return(bmi)
  }
  
  
# Call the function by naming arguments explicitly
calculate_BMI(weight = 60, height = 1.75)
  
# Call the function using variables as arguments
calculate_BMI(weight = weight, height = height)
  
# If a function expects two arguments, you must provide both
# This will give an error because 'height' is missing
calculate_BMI(60) 
  
# You can assign default values to function arguments
calculate_BMI <- function(weight, height = 1.75) {
    # Perform the BMI calculation
    bmi <- weight / (height ^ 2)
    
    # Return the BMI value
    return(bmi)
  }
  
# In this case, if you don’t provide height, R automatically uses 1.75 as the default.
calculate_BMI(weight = 60)
  
# ----------------------------
# Lazy evaluation in R
# ----------------------------
# If your function has three arguments, but the body only uses two,
# R does not force you to supply the third argument for the calculation.
# Example: 'age' is defined as an argument, but not used in the formula.
  
calculate_BMI <- function(weight, height, age) {
    # Perform the BMI calculation
    bmi <- weight / (height ^ 2)
    
    # Return the BMI value
    return(bmi)
  }
  
# Here we pass only 'weight' and 'height'
# Even though 'age' exists as an argument, it is ignored because it is not used
calculate_BMI(60, 1.65)
  
#---------
# Summary:
#---------
# Functions help us package logic once and apply it to different inputs.
  
  
# ----------------------------------
# 4. Automating Workflows with for-Loop
# ----------------------------------
# Suppose you have multiple datasets and you want to:
#   - import them,
#   - check missing values,
#   - clean columns,
#   - compute BMI,
#   - and save results.

# Instead of repeating steps for each file, we use loops.
  
# -----------------------
# Typical loop workflow:
# -----------------------

# Define the input folder (where raw data files are stored) and the output folder (where results will be saved).
# Specify the list of files that need to be processed.
# Prepare an empty list in R to store results for later use.
# For each file in the list:
#          Import the data into R.
#          Check and handle missing values (NA).
#          Calculate BMI using the calculate_BMI function.
#          Save the processed results both as a CSV file and in the R results list.

# ----------------------------------------------------
# Calculate BMI of two dataset within loop
  
# Define input and output folders
input_dir <- "raw_data" 
output_dir <- "results"
  
  
# create output folder if not already exist
  
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
# the Results folder will be created if already not exist
  
# List which files to process
files_to_process <- c("BMI_data_1.csv", "BMI_data_2.csv") 
# These must match exactly with the names in your working folder,
# otherwise R will not find them.
  
# Prepare empty list to store results in R 
result_list <- list()
  
# For each file with in a loop:
#       - import data
#       - handle NA values
#       - calculate BMI using calculate_BMI function
#       - save results (both as CSV and inside R list)
  
  
for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n") # n is number of files
    
  input_file_path <- file.path(input_dir, file_names)
    
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
    
  # handling missing values
    
  if("height" %in% names(data)){
  missing_count <- sum(is.na(data$height))
      
      cat("Missing values in 'height':", missing_count, "\n") # replace this n with the actual number of missing values
      data$height[is.na(data$height)] <- mean(data$height, na.rm = TRUE)
    }
    
    if("weight" %in% names(data)){
      missing_count <- sum(is.na(data$weight))
      
      cat("Missing values in 'weight':", missing_count, "\n")
      data$weight[is.na(data$weight)] <- mean(data$weight, na.rm = TRUE)
    }

    # calculate BMI
   data$bmi <- calculate_BMI(data$weight, data$height)
   cat("BMI has been calculated successfully.\n")
    
    # save results in R
    result_list[[file_names]] <- data 
    
    # save results in Results folder
    output_file_path <- file.path(output_dir, paste0("BMI_results", file_names))
    write.csv(data, output_file_path, row.names = FALSE)
    cat("Results saved to:", output_file_path, "\n")
  }
  
  # The loop repeats until all files are processed.
  # to extract the results
  results_1 <- result_list[[1]] 
  results_2 <- result_list[[2]]
  
  
  
  # --------
  # Summary:
  # --------
  # Loops automate repetitive work 
  # making your workflow faster 
  # consistent, and reproducible
  
save.image(file="Class_2.RData")

  load("Class_2.RData")  

  write.csv(data,file="clean_data/BMI_data_1.csv")  
  write.csv(data,file="clean_data/BMI_data_2.csv")
  