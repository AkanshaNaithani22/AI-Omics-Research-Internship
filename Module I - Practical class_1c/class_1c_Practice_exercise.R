# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#                Module I: Getting Started with R
# -------------------------------------------------------------------
# =================================================================== 
# -------------------------------------------------------------------
#                 Practice excercise_Class_1c
# -------------------------------------------------------------------
# ===================================================================


# 1. Check Cholesterol level (using if) 
# Write an If statement to check cholesterol level is greater than 240, 
# if true, it will prints “High Cholesterol”

cholesterol <- 230

if (cholesterol > 240){
  print("High Cholesterol")
}

# -------------------------------------------------------------------

# 2. Blood Pressure Status (using if...else)
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”

Systolic_bp <- 130

if (Systolic_bp < 120){
  print("Blood Pressure is normal")
}else{
  print("Blood Pressure is high")
}


# ----------------------------------
# A. PATIENT DATASET 
# ----------------------------------

# --------------------------------------------------------------

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and metadata.csv
# Perform the following steps separately on each dataset (patient_info.csv data and metadata.csv)
# Create a copy of the dataset to work on.
# Identify all columns that should be converted to factor type.
# Store their names in a variable (factor_cols).

raw_data1 <- read.csv(file.choose())
dataset1 <- raw_data1

str(raw_data1)
View(raw_data1)

raw_data1$patient_id <- NULL
raw_data1

# to convert columns of raw_data to factor 
# creating a copy of raw_data to clean_data 
# to save the output into clean_data

clean_data1 <- raw_data1
str(clean_data1)

# Columns that should be converted to factor
# gender: chr to factor
# diagnosis: chr to factor 
# smoker: chr to factor 

# creating single variable for these columns to use in "for loop" function
# to covert the R type together from character to factor

factor_cols1 <- c("gender", "diagnosis" ,"smoker")

for (col in factor_cols1) {
  if (is.character(clean_data1[[col]])) {
    clean_data1[[col]] <- as.factor (clean_data1[[col]])
  }
}
  
str(clean_data1)
clean_data1

# -----------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

binary_cols1 <- c("gender","smoker")   # stored column name in a vector

#column specific conversion
for (col in binary_cols1) {
  if (col == "gender"){
    clean_data1[[col]] <- ifelse (clean_data1[[col]] == "Female", 1,0)}
  else if (col == "smoker"){
    clean_data1[[col]] <- ifelse (clean_data1[[col]] == "Yes", 1,0) }
}

str(clean_data1)
print(clean_data1)

for (col in binary_cols1) {
  if (is.numeric(clean_data1[[col]])) {
    clean_data1[[col]] <- as.integer(clean_data1[[col]])
  }
}

str(clean_data1)
str(raw_data1)


#----------------------------------
# B. META DATASET
#----------------------------------

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and metadata.csv
# Perform the following steps separately on each dataset (patient_info.csv data and metadata.csv)
# Create a copy of the dataset to work on.
# Identify all columns that should be converted to factor type.
# Store their names in a variable (factor_cols).

raw_data2 <- read.csv(file.choose())
dataset2 <- raw_data2

str(raw_data2)
View(raw_data2)

raw_data2$name <- NULL
raw_data2

# creating a copy of raw_data to save the output
clean_data2 <- raw_data2
clean_data2
str(clean_data2)


# columns that need to be converted 
# height: chr to factor
# gender: chr to factor

# creating single variable 

factor_cols2 <- c("height","gender")

for (col in factor_cols2) {
  if (is.character(clean_data2[[col]])) {
    clean_data2[[col]] <- as.factor (clean_data2[[col]])
  }
}

str(clean_data2)

# ---------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

binary_cols2 <- c("height", "gender") #stored columns in vector


#column specific conversion
for (col in binary_cols2) {
  if (col == "height"){
    clean_data2[[col]] <- ifelse (clean_data2[[col]] == "Tall",0, 
                                  ifelse (clean_data2[[col]] == "Medium",1,2))}
  else if (col == "gender"){
    clean_data2[[col]] <- ifelse (clean_data2[[col]] == "Female", 1,0) }
}



for (col in binary_cols2) {
  if (is.numeric(clean_data2[[col]])) {
    clean_data2[[col]] <- as.integer(clean_data2[[col]])
  }
}

str(clean_data2)
str(raw_data2)

save.image(file = "Class_1c_Practice_exercise_full_workspace.RData")

