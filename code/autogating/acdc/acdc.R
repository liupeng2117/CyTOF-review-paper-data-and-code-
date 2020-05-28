setwd("D:/research/Liza/cytof/CyTOF existing methods and challenges/data")
load("ref_label.Rdata") # load reference cell type label
setwd("D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/acdc")
acdc_result<-read.csv("result.csv")
cluster_label<-acdc_result$predicted
library(plyr)
cluster_label<-mapvalues(cluster_label, from=c('unknown','naive b cells','memory b cells','classical monocytes','transitional monocytes',
                                                 'non-classical monocytes','basophils','early nks','late nks','pdc','mdc','cd8 naive','cd8 central memory',
                                                 'cd8 effector memory','cd8 terminal effector','cd4 naive','cd4 central memory','cd4 effector memory',
                                                 'cd4 terminal memory','cd4- mait/nkt','cd4-cd8-gd t cells'), to=0:20)


#PlotStars(FlowSOM_out2[[1]],backgroundValues = as.factor(FlowSOM_out2[[2]]))
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.813

# F measure 
nclust=length(unique(cluster_label))-1 #number of clusters
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
mean(f1.scores) #0.7675359