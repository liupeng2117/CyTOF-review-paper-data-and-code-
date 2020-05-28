library(mclust) # ARI calculation
#setwd("../Review/data/")

load("~/Methodology/201912_CytofClusterEvaluation/Review/data/ref_label_downsample.Rdata") # load reference cell type label
ref_label = ref_label_downsample
SWIFT_label <- read.delim("~/Methodology/201912_CytofClusterEvaluation/20200107_v1B/downsample.fcs.Cluster_Output_20clA.txt")
SWIFT_label = SWIFT_label[,3]
SWIFT_alt = read.delim("~/Methodology/201912_CytofClusterEvaluation/20200107_v1B/downsample.fcs.Cluster_Output_20clB.txt")
SWIFT_label_alt = SWIFT_alt[,3]

# ARI
cluster_label = factor(as.numeric(factor(SWIFT_label)))
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.05481306

# F measure 
nclust=length(levels(cluster_label))
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i])
    b<-names(ref_label[cluster_label==j & ref_label!=0])
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    f1.matrix[i,j]<-2/(1/precision +1/recall)
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.2873372



