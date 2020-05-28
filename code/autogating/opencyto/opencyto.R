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

### Opencyto
library(openCyto)
library(microbenchmark)
library(ncdfFlow)
library(flowWorkspace)
library(flowCore)
#create gating template
library(data.table)
gtFile <- "D:/research/Liza/cytof/CyTOF existing methods and challenges/programs/autogating/cd45.csv"
dtTemplate <- fread(gtFile)
dtTemplate
gt_cd45 <- gatingTemplate(gtFile, autostart = 1L)
gt_cd45
#load the raw data
ncfs  <- read.ncdfFlowSet(paste0(SampleName,"_01_CD45_CD66b-__Lymphocytes__DCs__Monocytes_.fcs"))
fr <- ncfs[[1]]
gs <- GatingSet(ncfs)
gs
#compensation
#compMat <- gh_get_compensations(gh)
#gs <- compensate(gs, compMat)
#transformation
chnls <- c("Y89Di", "Nd143Di", "Nd144Di", "Nd145Di", "Nd146Di", "Sm147Di", "Sm149Di", "Nd150Di", 
           "Eu151Di", "Sm154Di", "Gd155Di","Gd160Di", "Dy161Di", 
           "Dy163Di", "Dy164Di", "Er166Di", "Er167Di", "Er168Di", "Er170Di", "Yb171Di", "Yb173Di") 
trans <- estimateLogicle(gs[[1]], channels = chnls)
gs <- transform(gs, trans)
#gating
gt_gating(gt_cd45, gs)
autoplot(gs[[1]])
plot(gs[[1]])
autoplot(gs[[1]])

#the time depends on the gating methods users select, the ARI/F1 score also depends on the gating template users provide.
#this method only can gate two levels +/-, and no intermetiate
#this method can not make adjustment, not as flexible as manual gating