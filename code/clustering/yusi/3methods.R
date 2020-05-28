rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("C:/Users/acdsee8602/Documents/projects/collaboration/CYTOF_review/code/data")

library(openCyto)
library(PAC)
library(cytofkit)
library(Rclusterpp)
library(SamSPECTRAL)
library(flowCore) # for loading fcs files
#library(cytofkit2) #Phenograph and other clustering methods
library(FlowSOM) #FlowSOM
library(immunoClust)
library(mclust) # ARI calculation
library(pryr)

#####loading data
SampleName<-"HuImmProfiling_S1_PBMC_1"
dt<-read.FCS(paste0(SampleName,"_01_CD45_CD66b-__Lymphocytes__DCs__Monocytes_.fcs"),transformation = FALSE)
load("ref_label.Rdata") # load reference cell type label


### preprocessing --------------------------------------------
#Extract markers
markers<- attributes(ref_label)$markers$name
dt2<-dt[,markers]

#asinh transformation to normalize the data
asinhTrans=arcsinhTransform(a=0, b=0.2)
translist<-transformList(markers, asinhTrans)
dt2.transform<-transform(dt2, translist)

dt2.transform<-exprs(dt2.transform)

###testing immunoClust-------------------------------------------------------

immunoClustOut=cell.MajorIterationLoop(dt2.transform)
immunoClustOut=test
library(microbenchmark)
microbenchmark(
  "immunoClust"={immunoClustOut2<-cell.MajorIterationLoop(dt2.transform)},
  
  times=10)

cluster_label<-immunoClustOut@label

ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.29

# F measure 
nclust=length(levels(factor(cluster_label))) #number of clusters
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
mean(f1.scores) #0.474514

###testing PAC--------------------------------------------------
PACout=PAC(dt2.transform,K=20)
microbenchmark(
  "PAC"={PACout2<-PAC(dt2.transform,K=20)},
  times=10)
#Unit: seconds
#expr      min       lq     mean   median       uq      max neval
#PAC 27.60893 27.63004 27.93816 27.80045 28.03572 29.06105    10
cluster_label<-PACout

ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.7825

# F measure 
nclust=length(levels(factor(cluster_label))) #number of clusters
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
mean(f1.scores) #0.7425


###testing SamSPECTRAL-----------------------------------------------

SamSPECTRALout=SamSPECTRAL(dt2.transform,k.for_kmeans = 20,number.of.clusters = 20,normal.sigma = 200,separation.factor = 0.7)

microbenchmark(
  "SamSPECTRAL"={SamSPECTRALout2<-SamSPECTRAL(dt2.transform,k.for_kmeans = 20,number.of.clusters = 20,normal.sigma = 200,separation.factor = 0.7)},
  times=10)

cluster_label<-SamSPECTRALout

ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.5701382

# F measure 
nclust=length(levels(factor(cluster_label))) #number of clusters
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
mean(f1.scores) #0.33624


### determine the normal.sigma
run.live <- TRUE
full=dt2.transform

for (i.column in 1:dim(full)[2]){
   ith.col <- full[,i.column]
   full[,i.column] <- (ith.col-min(ith.col)) /(max(ith.col) - min(ith.col) )
   
}
maximum.number.of.clusters <- 30;precision <- 6
space.length <- 1
community.weakness.threshold <-1
society <-Building_Communities(full,m=3000, space.length, community.weakness.threshold)
normal.sigma <- 10

conductance <- Conductance_Calculation(full, normal.sigma, space.length,
                                        society, precision)


if (run.live){
  
   clust_result.10 <- Civilized_Spectral_Clustering(full, maximum.number.of.clusters,
                                     society, conductance,stabilizer=1)
   eigen.values.10 <- clust_result.10@eigen.space$values
   } else data("eigen.values.10")
plot(eigen.values.10[1:50])



#-----------------------#
normal.sigma <- 200
conductance <- Conductance_Calculation(full, normal.sigma, space.length,
                                       society, precision)


if (run.live){ clust_result.1000 <-Civilized_Spectral_Clustering(full, maximum.number.of.clusters,
                                   society, conductance,stabilizer=1)
   eigen.values.1000 <- clust_result.1000@eigen.space$values
   } else data("eigen.values.1000")

plot(eigen.values.1000[1:50])

#####citrus--------------------------------------------------#
library(citrus) ###for at least two samples for DE testing 
citrus.launchUI()
#####fusionclust #####too long for testing -----------------------------------#run 18h and still not finished in server speed2
library(fusionclust)
x=dt2.transform
scores<- cosci_is(x,0)
plot(features,scores,type="p",col="red")
fusionclustOut=bmt(dt2.transform)
#too slow

####cydar----------------------------------------------------#bugs
library(cydar)
library(ncdfFlow)
collected.exprs=read.ncdfFlowSet(paste0(SampleName,"_01_CD45_CD66b-__Lymphocytes__DCs__Monocytes_.fcs"))
pool.ff <- poolCells(collected.exprs)
trans <- estimateLogicle(pool.ff, colnames(pool.ff))
proc.ff <- transform(pool.ff, trans)
processed.exprs <- transform(collected.exprs, trans)
cd <- prepareCellData(processed.exprs)
#Error in get(name, envir = asNamespace(pkg), inherits = FALSE) : 
#object 'normalize_names_replacement_value' not found
#looks like this is due to the update of packages "S4Vectors"


####densitycut-------------------------------------------------#
library(densitycut)
densitycutOut=DensityCut(dt2.transform)
library(microbenchmark)
microbenchmark(
  "DensityCutOut"={densitycutOut2<-DensityCut(dt2.transform)},
  
  times=10)
cluster_label=densitycutOut$cluster
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.77837 # 26 clusters


# F measure 
nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
    f1.matrix[i,j]<-2/(1/precision +1/recall)}
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.344

####flowClust--------------------------------------------
library(flowClust)
flowClustOut=flowClust(dt2.transform,trans = 0,K=20)
library(microbenchmark)
microbenchmark(
  "flowClustOut"={flowClustOut2<-flowClust(dt2.transform,trans = 0,K=20)},
  
  times=10)
#Unit: milliseconds
#expr      min       lq   mean   median       uq      max neval
#flowClustOut 452.9468 479.8633 531.62 487.9062 501.0517 941.0919    10
cluster_label=flowClustOut@label
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.5514

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
      }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.00001135112

####flowMeans-------------------------------------------------------------------------
library(flowMeans)

flowMeansOut=flowMeans(dt2.transform,NumC = 20)
library(microbenchmark)
microbenchmark(
  "flowMeansOut"={flowMeansOut2<-flowMeans(dt2.transform,NumC = 20)},
  
  times=10)
cluster_label=flowMeansOut@Label
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.8013357

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.7126421
####flowMerge-------------------------------------------------------#bugs
library(flowMerge)
flowClust.res=flowClust(dt2,trans = 1,K=15:25)
flowClustMaxBIC=flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]]
flowClust.obj=flowObj(flowClustMaxBIC,dt2.transform)
flowClust.merge=merge(flowClust.obj,metric="entropy")
i=fitPiecewiseLinreg(flowClust.merge)
####flowPeaks--------------------------------------------------------#
library(flowPeaks)
flowPeaksOut=flowPeaks(dt2.transform)
library(microbenchmark)
microbenchmark(
  "flowPeaksOut"={flowPeaksOut2<-flowPeaks(dt2.transform)},
  
  times=10)
cluster_label=flowPeaksOut$peaks.cluster
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.6461867

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.5493605



####cytometree-------------------------------------------------------#
library(cytometree)
cytometreeOut=CytomeTree(dt2.transform)
library(microbenchmark)
microbenchmark(
  "cytometreeOut"={cytometreeOut2<-CytomeTree(dt2.transform)},
  
  times=1)

cluster_label=cytometreeOut2$labels
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.0765

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.2087


####DepecheR------------------------------------------#
library(DepecheR)
depecheOut=depeche(dt2.transform,k=20)
library(microbenchmark)
microbenchmark(
  "depecheOut"={depecheOut2=depeche(dt2.transform,k=20)},
  
  times=10)

cluster_label=depecheOut$clusterVector
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.7539315

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.5363942

####kmeans--------------------------------------------------#
kmeansOut=kmeans(dt2.transform,centers = 20)

library(microbenchmark)
microbenchmark(
  "kmeansOut"={kmeansOut2=kmeans(dt2.transform,centers = 20)},
  
  times=10)

#Unit: seconds
#expr      min       lq     mean   median       uq      max neval
#kmeansOut 1.269717 1.388109 1.636957 1.499338 1.743531 2.738904    10

cluster_label=kmeansOut$cluster
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.7014708

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.6848694

####RclusterCpp--------------------------------------------------------------#
library(Rclusterpp)
RclusterppOut=Rclusterpp.hclust(dt2.transform)

library(microbenchmark)
microbenchmark(
  "RclusterppOut"={RclusterppOut2=Rclusterpp.hclust(dt2.transform)},
  
  times=10)


cluster_label=cutree(RclusterppOut,20)
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.7064929

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.7148727


####diffcyt---------------------------------------------
library(diffcyt) #for DE testing; only use flowSOM
library(SummarizedExperiment)
flowset_dt=flowSet(dt2)   ####bugs
dt_list=as(dt2,"list")
marker_info=data.frame(marker_name=colnames(exprs(dt2)),marker_class=rep("none",21))
experiment_info=data.frame(sample_id=1)
d_se=prepareData(flowset_dt,marker_info = marker_info,experiment_info = experiment_info)
d_se <- transformData(d_se)

diffcytOut=generateClusters(d_se,meta_k = 20)

library(microbenchmark)
microbenchmark(
  "diffcytOut"={diffcytOut2=generateClusters(d_se,meta_k = 20)},
  
  times=10)
#Unit: seconds
#expr            min       lq     mean   median      uq      max neval
#diffcytOut 2.154183 2.171923 2.184548 2.181377 2.18729 2.240243    10
cluster_label=rowData(diffcytOut)@listData$cluster_id
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) ####bugs 0

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.09

######ClusterX----------------------------------------------#limited memory issue
clusterOut=ClusterX(dt2.transform)
#load("ClusterX.RData")

library(microbenchmark)
microbenchmark(
  "clusterOut"={clusterOut2=ClusterX(dt2.transform)},
  
  times=10)


cluster_label=clusterOut$cluster
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) ####0.2485708

nclust=length(levels(factor(cluster_label))) #number of clusters
f1.matrix<-matrix(,nrow=20, ncol=nclust)
for(i in 1:20){
  for(j in 1:nclust){
    a<-names(ref_label[ref_label==i]) #true
    b<-names(ref_label[cluster_label==j & ref_label!=0]) #predicted, cells labeled with 0 are not included in evaluation
    recall<-length(intersect(a,b))/length(a) 
    precision<-length(intersect(a,b))/length(b) 
    if(recall==0 | precision==0){f1.matrix[i,j]=0}else{
      f1.matrix[i,j]<-2/(1/precision +1/recall)
    }
  }
}
f1.scores<-apply(f1.matrix,1,max)
names(f1.scores)<-attributes(ref_label)$annotation$cell_type
mean(f1.scores) #0.22234



######DensVM------------------------------------------------#

densVmOut=DensVM(dt2.transform,exprs(dt2))#################?bugs
