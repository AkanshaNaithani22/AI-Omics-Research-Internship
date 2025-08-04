#Subfolders 
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")


#Import data from csv file
data<-read.csv(file.choose())
View(data)
str(data)


#Convert 'gender' column to factor
data$gender_fac<-as.factor(data$gender)
str(data)

#Convert "diagnosis' column to factor
data$diagnosis_fac<-as.factor(data$diagnosis)
class(data$diagnosis)
data$diagnosis

#Set factor level order manually
data$diagnosis<-factor(data$diagnosis,
                       levels=c("normal","cancer"))
data$diagnosis
str(data)

#Convert 'smoker' column to factor
data$smoker_fac<-as.factor(data$smoker)
str(data)

#Convert factor to numeric using ifelse statement (Yes = 1, No = 0)
data$smoker_num<-ifelse(data$smoker_fac=="Yes",1,0)
class(data$smoker_num)

#Convert numeric smoker code to factor 
data$smoker_num<-as.factor(data$smoker_num)
class(data$smoker_num) 

#Save file as csv in clean_data folder 
write.csv(data,file="clean_data/patient_info_clean.csv")

#Save the entire R workspace
save.image(file="Akansha_Naithani_Class_1b_Assignment.RData")

load("Akansha_Naithani_Class_1b_Assignment.RData")
