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
write.FCS(dt2.transform, "D:/CyTof_review_Xiangning/FLOCK/selectedfcs.fcs")

dt2.transform<-exprs(dt2.transform)
dt2.notrans = exprs(dt2)

write.table(dt2.transform, "D:/CyTof_review_Xiangning/FLOCK/input.txt", row.names = FALSE, col.names = TRUE)
write.table(dt2.transform, "D:/CyTof_review_Xiangning/FLOCK/input.flowtxt", row.names = FALSE, col.names = TRUE)

#see if the result is reproducible
result1 = read.table("D:/CyTof_review_Xiangning/FLOCK/Galaxy25-[flock2_with_mfi_on_No_Transformation_selectedfcs.fcs].flowclr", header = TRUE)
result2 = read.table("D:/CyTof_review_Xiangning/FLOCK/Galaxy31-[flock2_with_mfi_on_No_Transformation_selectedfcs.fcs].flowclr", header = TRUE)
all(result1[, 22]==result2[, 22])
result3 = read.table("D:/CyTof_review_Xiangning/FLOCK/Galaxy28-[flock2_with_mfi_on_No_Transformation_selectedfcs.fcs].flowclr", header = TRUE)
all(result1[, 22]==result3[, 22])
table(result3[, 22])
result4 = read.table("D:/CyTof_review_Xiangning/FLOCK/Galaxy34-[flock2_with_mfi_on_No_Transformation_selectedfcs.fcs].flowclr", header = TRUE)
table(result4[, 22])

# ARI
cluster_label<-factor(result4[, 22])
table(cluster_label)
#PlotStars(FlowSOM_out2[[1]],backgroundValues = as.factor(FlowSOM_out2[[2]]))
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
adjustedRandIndex(ref_label[ids],cluster_label[ids]) #0.732

# F measure 
nclust=length(levels(cluster_label)) #number of clusters
f1.matrix<-matrix(, nrow=20, ncol=nclust)
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
mean(f1.scores) #0.649

