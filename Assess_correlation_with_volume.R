##Read the radiomic features extracted from the ROIs from a csv file
Data<- read.csv()
##Load/create the list of reproducible features
Reproducible_features<- c() 
##Or
Reproducible_features<-read.csv()

##Create correlation matrix between the reproducible features and the volume of the ROI
Cor_Mat<- cor(Data[,Reproducible_features], Data$Shape_volume, use='ev')