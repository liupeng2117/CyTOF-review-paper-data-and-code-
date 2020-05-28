library(flowCore) # for loading fcs files
#library(cytofkit2) #Phenograph and other clustering methods
library(FlowSOM) #FlowSOM
library(mclust) # ARI calculation
setwd("~/data/")

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

### FlowSOM --------
FlowSOM_out<-FlowSOM(dt2.transform, colsToUse=markers, nClus=20)
#Phenograph
#phenograph_out<-Rphenograph(dt2.transform, k = 30)

### benchmark ------------

# time
library(microbenchmark)
microbenchmark(
  "FlowSOM"={FlowSOM_out2<-FlowSOM(dt2.transform, colsToUse=markers, nClus=20)},
  #"Phenograph"={phenograph_out<-Rphenograph(dt2.transform, k = 30)},
  times=10)
#Unit: seconds
#   expr     min       lq     mean   median       uq      max neval
#FlowSOM 10.0828 10.10522 10.23237 10.15262 10.25084 10.74206    10

# ARI
cluster_label<-GetMetaclusters(FlowSOM_out, FlowSOM_out$metaclustering)
#PlotStars(FlowSOM_out2[[1]],backgroundValues = as.factor(FlowSOM_out2[[2]]))
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.620

# F measure 
nclust=length(levels(cluster_label)) #number of clusters
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
mean(f1.scores) #0.672
