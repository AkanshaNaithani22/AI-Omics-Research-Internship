


# =======================================================================
#                  AI and Biotechnology/Bioinformatics
# =======================================================================
# -----------------------------------------------------------------------
#                  AI and Omics Research Internship 2025
# -----------------------------------------------------------------------
#                  Module I: Getting started with R (Class 1c)
# -----------------------------------------------------------------------
# =======================================================================


# Topic: 
# Syntax in R: 
     # 1. Variables
     # 2. Comments 
     # 3. Keywords

# Variables
# stores values in R
# <- assignment opperator 

# gene name "TP53"

gene <- "TP53"

# to reterive value in console 
gene

print (gene)


# 2.3, 4.6, 3.6, 7.2, 4.7
# to store these values in one variable 

expression_levels <- c(2.3, 4.6, 3.6, 7.2, 4.7)

# To import data as variable
raw_data <- read.csv(file.choose())
raw_data

# Rules 

# variable name must start with letter 
1gene <- 25 # variable name can't start with number 

gene1 <- 25 # can add number at the end

# no spaces allowed in variable name
sample id <- "s01"

# Instead of spaces use underscore (_) or dot (.)
sample.id <- "s01"
sample_id <- "s02"

# R is case sensitive 
Glucose_level <- 110

glucose_level <- 110

# R overwrite variables without any warning
glucose_level <- c(110, 90, 120)

data <- raw_data # create a copy of my raw_data

raw_data$patient_id <- NULL # this code will remove patient_id column from raw_data

raw_data

# For data cleaning and transforming create a new variable for that data

clean_data <- data [,-1]
# this code deleted the patient_id column & it assign a new variable 

#Comments 
# help to understand your code 
# comments are for our own understanding R doesn't consider it as code 
# data_2 <- 23
data_2 <- 23

# pro tip: turn comments into heading
#### Heading 1 ####
#### Heading 2 ####

# Keywords
# These are reserved words in R for specific function
# if, else, TRUE, FALSE, for, so on...

help ("reserved")
help(mean)
?median

# sort values from largest to lowest 
sorted_age <- sort (raw_data$age, decreasing = TRUE)
sorted_age
raw_data$age

# sort values smallest to largest 
sorted_age2 <- sort (raw_data$age, decreasing = FALSE)
sorted_age2

# if & else, which are used for creating logical conditions

gene_expression <- 30
 
if (gene_expression > 50){
  print("Gene is highly expressed")
}

# here "if" is the keyword that check the condition if gene_expression > 25
# condition is true in this case 
# incase if the condition is false
if (gene_expression > 50){
  print("Gene is highly expressed")
}else{
  print("Gene expression is low")
}

# you cannot use keywords as variable names
if <- 28

# for loop: used to repeat same tasks multiple times 
# lets say we want to convert data type of multiple in of our 
str(raw_data)
# gender is categorical data type
# It should be in a factor format 
# gender column from chr to factor
# diagnosis: cancer/normal it is also a categorical variable 
# smoker: chr to factor 

# Instead of manually conversion we will use this for loop function for automatic conversion 
# of all these 3 columns with one command

# To convert raw_data column into factor 
# I want to save output in clean_data

# Create a copy of raw_data with name clean_data
clean_data <- raw_data
str(clean_data)

# to convert column automatically into factor
# create for loop function
for (i in 1:ncol (clean_data)) { #can put anything instead of "i"
  if (is.character(clean_data[[i]])) {
    clean_data[[i]] <- as.factor(clean_data[[i]])
  }
}

str(clean_data)
str(raw_data)

save.image(file="Akansha_Naithani_class_1c.RData")

load("Akansha_Naithani_class_1c.RData")












