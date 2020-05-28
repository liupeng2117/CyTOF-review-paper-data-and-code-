library(flowCore) # for loading fcs files
#library(cytofkit2) #Phenograph and other clustering methods
#BiocManager::install("FlowSOM")
library(FlowSOM) #FlowSOM

library(mclust) # ARI calculation
datafile = "D:/UPitts/Research/cytof/data/"
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
dt2.notrans = exprs(dt2)

dt3.transform<- dt2.transform[1:5000, ]
dt3.notrans = dt2.notrans[1:5000, ]
#1. download the package zip from the supplementary
packagefile = "D:/UPitts/Research/cytof/ccast/"
#install prerequisite packages
#install.packages("mixtools")
#install.packages("party")
#BiocManager::install("RBGL")
#install.packages(paste0(packagefile, "Algorithm_S1.gz"),repos = NULL, type="source") ## within R
library(CCAST)
library(MASS) #needed function
outfile = "D:/UPitts/Research/cytof/ccast/output/"
dir.create(outfile)
setwd(outfile)

colid = 1:dim(dt2.transform)[2]
#use dt2.transform (asinh transformation)
Result2<-ccast_main(file=dt3.transform, transformlogic=FALSE, asinhp=1, deterministic=TRUE,
                    colid, coln=NULL, rown=NULL, npmix=FALSE, k=20, boot=NULL, ylabel="Y89Di",
                    origscale=TRUE, groups=NULL,runsilhouette=TRUE)
 
#set the result to 
time1 = Sys.time()
outfile2 = "D:/UPitts/Research/cytof/ccast/output2/"
dir.create(outfile2)
setwd(outfile2)
#use dt2.transform (asinh transformation)
Result2<-ccast_main(file=dt2.notrans, transformlogic=FALSE, asinhp=150, deterministic=TRUE,
                    colid, coln=NULL, rown=NULL, npmix=FALSE, k=20, boot=NULL, ylabel="Y89Di",
                    origscale=TRUE, groups=NULL,runsilhouette=TRUE)
time2 = Sys.time()

#set the result to 
(time1 = Sys.time())
outfile2 = "D:/UPitts/Research/cytof/ccast/output2/"
dir.create(outfile2)
setwd(outfile2)
#use dt2.transform (asinh transformation)
Result2<-ccast_main(file=dt3.transform, transformlogic=FALSE, asinhp=1, deterministic=TRUE,
                    colid, coln=NULL, rown=NULL, npmix=FALSE, k=20, boot=NULL, ylabel="Y89Di",
                    origscale=TRUE, groups=NULL,runsilhouette=TRUE)
(time2 = Sys.time())


