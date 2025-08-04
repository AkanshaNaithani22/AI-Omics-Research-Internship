print("Hi")
data<-1:10
data
dog=readxl::read_excel(file.choose())
colnames(dog)
row.names(dog)

dog$...3
sample_names=dog$...3
sample_names

my_data=dog[,-1]
my_data=dog[,-1:-5]
my_data_selected_columns=dog[1,1:5]
colnames(dog)[1]="sample_names"

install.packages("dplyr")

library(dplyr)

dog = dog %>% rename_with(
  ~gsub("...","NEW", .x),
  .cols = contains ("...")
)

load("Akansha_Naithani_Class_1b_Assignment.RData")
