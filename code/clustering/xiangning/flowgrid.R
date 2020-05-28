library(flowCore) # for loading fcs files
#library(cytofkit2) #Phenograph and other clustering methods
#BiocManager::install("FlowSOM")
library(FlowSOM) #FlowSOM

library(mclust) # ARI calculation
datafile = "D:/CyTof_review_Xiangning/Data/"
setwd(datafile)

### load CyTOF data and manual gating reference label ----
SampleName<-"HuImmProfiling_S1_PBMC_1"
dt<-read.FCS(paste0(SampleName,"_01_CD45_CD66b-__Lymphocytes__DCs__Monocytes_.fcs"),transformation = FALSE)
load("ref_label.Rdata") # load reference cell type label

### preprocessing ----- 
#Extract markers
markers<- attributes(ref_label)$markers$name
dt2<-dt[,markers]

#asinh transformation to normalize the data
asinhTrans=arcsinhTransform(a=0, b=0.2)
translist<-transformList(markers, asinhTrans)
dt2.transform<-transform(dt2, translist)
dt2.transform<-exprs(dt2.transform)
label = read.csv("ref_label.csv", row.names = 1)
write.csv(dt2.transform, "D:/CyTof_review_Xiangning/flowgrid/input.csv", row.names = FALSE)
write.csv(label, "D:/CyTof_review_Xiangning/flowgrid/label.csv", row.names = FALSE)

label = read.csv("D:/UPitts/Research/cytof/xxn/out.csv")
#label = read.csv("D:/CyTof_review_Xiangning/flowgrid/out.csv")
label = label$X3
cluster_label = label[label!= -1]
cluster_label = factor(cluster_label)
# ref_label0 = read.csv("D:/CyTof_review_Xiangning/flowgrid/label.csv")
# ref_label0 = ref_label0$label
# ref_label0 = ref_label0[label!= -1]
# ref_label = ref_label0
load("D:/UPitts/Research/cytof/data/ref_label.Rdata") # load reference cell type label
ref_label = ref_label[-1]
#ref_label = ref_label[-length(ref_label)]
ref_label = ref_label[label!=-1]
nclust=length(levels(cluster_label)) #number of clusters
f1.matrix<-matrix(, nrow=20, ncol=nclust)
#for(i in 1:20){
#  for(j in 1:nclust){
#install.packages("mclust")
library(mclust) # ARI calculation
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.536

for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    f1.matrix[i,j]<-2/(1/precision +1/recall)
  }
}

f1.scores<-apply(f1.matrix,1,max, na.rm = TRUE)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores, na.rm = TRUE) #0.482

na.row = apply(f1.matrix, 1, function(x){mean(is.na(x))})

#delete the first
#delete the last