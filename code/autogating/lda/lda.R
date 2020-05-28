dt2.df<-data.frame(dt2.transform)
dt2.df<-data.frame(dt2.df,cell_type=ref_label)
write.csv(dt2.df,"D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/lda/data/data.csv")
setwd('D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/lda/')
AML.data <-read.csv('data/data.csv',header = TRUE)
head(AML.data)
#AML.data <- AML.data[AML.data$cell_type != '0',]
dim(AML.data)
library(caret)
set.seed(123)
Folds <- createDataPartition(AML.data$cell_type,2)

AML.Train <- AML.data[unlist(Folds[1],use.names = FALSE),]
AML.Test <- AML.data[unlist(Folds[2],use.names = FALSE),]

write.table(AML.Train,file = 'data train/data_train.csv',col.names = FALSE,row.names = FALSE,sep = ',')
write.table(AML.Test,file = 'data test/data_test.csv',col.names = FALSE,row.names = FALSE,sep = ',')

setwd('D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/lda/')
source('CyTOF_LDAtrain.R')
source('CyTOF_LDApredict.R')




microbenchmark(
  "immunoClust"={
    # cell type labels are in column no. 23 
    LDA.Model <- CyTOF_LDAtrain(TrainingSamplesExt = 'data train/',TrainingLabelsExt = '',mode = 'CSV',
                                             RelevantMarkers =  c(2:22),LabelIndex = 23, Transformation = 'arcsinh')
  
    Predictions <- CyTOF_LDApredict(LDA.Model,TestingSamplesExt = 'data test/', mode = 'CSV', RejectionThreshold = 0)
    },
  times=10)
#Unit: seconds
#       expr      min       lq     mean   median       uq      max neval
#immunoClust 6.541852 7.012394 7.313011 7.228609 7.684723 8.048462    10

Predictions <- unlist(Predictions)
save(Predictions,Folds, file="D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/lda/result.rdata")
Accuracy <- sum(Predictions == AML.Test$cell_type)/length(AML.Test$cell_type) * 100
Accuracy #86.39686

cluster_label<-as.numeric(Predictions)

ref_label<-ref_label[unlist(Folds[2])]
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
#ids2<-ref_label[unlist(Folds[2])]!=0
adjustedRandIndex(ref_label[ids],cluster_label[ids]) # 0.9120301

# F measure 
nclust=length(levels(factor(cluster_label)))-1 #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    f1.matrix[i,j]<-2/(1/precision +1/recall)
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.9155761


