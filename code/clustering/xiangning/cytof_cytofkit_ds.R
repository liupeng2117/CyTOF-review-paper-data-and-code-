#BiocManager::install("cytofkit")#not available for current R version
#1. install packages
#packagefile = "D:/CyTof_review_Xiangning/cytofkit/"
#install.packages(paste0(packagefile, "cytofkit_1.1.0.tar.gz"), type = "source", repos = NULL)
#install.packages("reshape")
library(cytofkit)
library(flowCore) # for loading fcs files
#library(cytofkit2) #Phenograph and other clustering methods
#BiocManager::install("FlowSOM")
library(FlowSOM) #FlowSOM
library(mclust) # ARI calculation
datafile = "D:/CyTof_review_Xiangning/Data/"
datafile = "D:/UPitts/Research/cytof/data/"
setwd(datafile)

### load CyTOF data and manual gating reference label ----
#SampleName<-"HuImmProfiling_S1_PBMC_1"
dt<-read.FCS("downsample.fcs",transformation = FALSE)
### load downsampled manual gating reference label ----
load("ref_label_downsample.Rdata") # load reference cell type label

### preprocessing ----- 
#Extract markers
markers<- attributes(ref_label)$markers$name
load("ref_label_downsample.Rdata") # load reference cell type label
dt2<-dt[,markers]
ref_label = ref_label_downsample

#asinh transformation to normalize the data
asinhTrans=arcsinhTransform(a=0, b=0.2)
translist<-transformList(markers, asinhTrans)
dt2.transform<-transform(dt2, translist)
dt2.transform<-exprs(dt2.transform)
#dt2.notrans = exprs(dt2)

setwd('D:/CyTof_review_Xiangning/cytofkit/DensVMds/')

(time0 = Sys.time())
dt2.transform.tsne <- cytof_dimReduction(dt2.transform, method = "tsne")
(time1 = Sys.time())
(time1-time0)
cluster_DensVM = densVM_cluster(xdata = dt2.transform, ydata = dt2.transform.tsne)
#cluster_DensVM <- cytof_cluster(xdata = dt2.transform, 
#                                ydata = dt2.transform.tsne, method = "DensVM")
(time2 = Sys.time())
time2 - time1

cluster_label<-cluster_DensVM[["clusters"]]
cluster_label = cluster_label$cluster

# ARI
cluster_label<-GetMetaclusters(FlowSOM_out, FlowSOM_out$metaclustering)
#PlotStars(FlowSOM_out2[[1]],backgroundValues = as.factor(FlowSOM_out2[[2]]))
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.7136

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
mean(f1.scores) # 0.6888862
save.image("result.rData")

#benchmark
densVMWrap = function(data = dt2.transform){
  dt2.transform.tsne <- cytof_dimReduction(data, method = "tsne")
  cluster_DensVM = densVM_cluster(xdata = data, ydata = dt2.transform.tsne)
}
library(microbenchmark)
microbenchmark(
  "densVM"={FlowSOM_out2<-densVMWrap(data = dt2.transform)},
  #"Phenograph"={phenograph_out<-Rphenograph(dt2.transform, k = 30)},
  times=10)

lab1 = cluster_DensVM$clusters$cluster
lab2 = cluster_DensVM$clusters$cluster
all(lab1==lab2)


