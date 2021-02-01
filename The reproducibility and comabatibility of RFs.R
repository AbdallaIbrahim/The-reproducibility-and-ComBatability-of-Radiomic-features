library(caret)
library(epiR)
library(sva)

#LOAD THE DATA FILE
Data<- read.csv("Pyradiomics_CCR1_OriginalFts.csv",header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors = F)

#DEFINE THE BATCH BASED ON WHICH SCANS WILL BE ASSESSED AND HARMONIZED
A<- unique(Data$Batch)
B<- unique(Data$Batch)

#CREATE THE MATRICES TO RECORD THE NUMBER OF FEATURES IN EACH PAIRWISE SCENARIO
Mat_stable<- as.data.frame(matrix(0, ncol=length(A), nrow= length(B)), row.names = as.character(B))
colnames(Mat_stable)<- c(A)
Mat_ComBat<-as.data.frame(matrix(0, ncol=length(A), nrow= length(B)), row.names = as.character(B))
colnames(Mat_ComBat)<- c(A)

#RUN THE CCC AND COMBAT ON EACH PAIR
for (i in 1:length(A)){
  for (j in 1:length(B)){
    
    if (A[i]<B[j]){
      print(paste(paste("Iteration no ", i, sep=""), paste("Row no ", j, sep="")))
      Batches<- c(A[i], B[j])
      Features_All<- colnames(Data[,6:(dim(Data)[2])])
      
      
      OCC_All<- as.data.frame(matrix(0, ncol=2, nrow= length(Features_All)))
      OCC_All$V1<- as.character(Features_All)
      
      for (l in 1:length(Features_All)){
        Data2<- as.data.frame(matrix(0, ncol=(length(Batches)+1), nrow=160))
        Data2$V1<- as.character(Features_All[l])
        for (m in 2:dim(Data2)[2]){
          k<- m-1
          Data2[,m]<- Data[Data$Batch==Batches[k],Features_All[l]]
        }
        
        
        CCC<- epi.occc(Data2[,2:(length(Batches)+1)], pairs = T)
        OCC_All[l,2]<-  CCC$pairs
        
      }
      
      
      Stable_All_OCCC<- c()
      for (n in 1:length(Features_All)){
        if(!(is.na(OCC_All[n,2]))){
          if (OCC_All[n,2]=='NaN'){
            Stable_All_OCCC<- append(Stable_All_OCCC,OCC_All[n,1])
          }
          else if (OCC_All[n,2]>0.90){
            Stable_All_OCCC<- append(Stable_All_OCCC,OCC_All[n,1])
          }
        }
        }
        
      Mat_stable[j,i]<- length(Stable_All_OCCC)
      #START THE COMBAT PROCESS
      Harmonize<-Data[Data$Batch %in% Batches, ]
      #REMOVE VARIABLES WITH (Near)ZERO VARIANCE
      badCols <- nearZeroVar(Harmonize)
      to_remove<- c()
      if(length(badCols)>0){
        for (z in 1:length(badCols)){
          if (badCols[z]>5){
            to_remove<-append(to_remove, badCols[z]) 
          }
        }
      }
      
      if (length(to_remove > 0)){
        Harmonize <- Harmonize[, -to_remove]
      }
      
      
      #REMOVE FEATURES WITH NA VALUES
      Harmonize<- Harmonize[, colMeans(is.na(Harmonize)) == 0]
      
      tryCatch({
        Harmonized<- as.data.frame(t(ComBat(as.matrix(t(Harmonize[,6:(dim(Harmonize)[2])])), batch = Harmonize$Batch, mod=NULL, par.prior = T)))
        
      }, error=function(e){})
      if (!(exists("Harmonized"))){
        Mat_ComBat[j,i]<- as.character('Error')
        
        
      }
      
      else if (exists('Harmonized') & (dim(Harmonized[,colMeans(is.na(Harmonized)) == 0])[2])>0){
        Features_All<- colnames(Harmonized)
        
        Harmonized<- cbind(Harmonize[,1:6], Harmonized)
        
        OCC_All_ComBatted<- as.data.frame(matrix(0, ncol=2, nrow= length(Features_All)))
        OCC_All_ComBatted$V1<- as.character(Features_All)
        
        for (h in 1:length(Features_All)){
          Data2<- as.data.frame(matrix(0, ncol=(length(Batches)+1), nrow=160))
          Data2$V1<- as.character(Features_All[h])
          for (o in 2:dim(Data2)[2]){
            v<- o-1
            Data2[,o]<- Harmonized[Harmonized$Batch==Batches[v],Features_All[h]]
          }
          
          
          CCC<- epi.occc(Data2[,2:(length(Batches)+1)], pairs= T)
          OCC_All_ComBatted[h,2]<-  CCC$occc
          
        }
        
        
        ComBatable_All_OCCC<- c()
        for (q in 1:length(Features_All)){
          if (OCC_All_ComBatted[q,2]=='NaN'){
            ComBatable_All_OCCC<- append(ComBatable_All_OCCC,OCC_All_ComBatted[q,1])
          }
          else if (OCC_All_ComBatted[q,2]>0.9){
            ComBatable_All_OCCC<- append(ComBatable_All_OCCC,OCC_All_ComBatted[q,1])
          }
        }
        Mat_ComBat[j,i]<- paste(as.character(length(ComBatable_All_OCCC)), as.character(length(Features_All)), sep = "/")
        rm("Harmonized")
        
      }
     
    
        
      }
      
  }
  
      
 }
    
#SAVE THE MATRICES TO CSV FILES
write.csv(Mat_stable, "stable_CCR2_pairs_RadiomiX_0.9.csv") 
write.csv(Mat_ComBat, "ComBat_CCR1_pairs_RadiomiX_0.9.csv")
