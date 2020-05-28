library(flowCore)
library(mclust) # ARI calculation
setwd("D:/research/Liza/cytof/CyTOF existing methods and challenges/data/")
SampleName<-"HuImmProfiling_S1_PBMC_1"
total.hul<-read.FCS(paste0(SampleName,"_01_CD45_CD66b-__Lymphocytes__DCs__Monocytes_.fcs"),transformation = T)
setwd("D:/research/Liza/cytof/CyTOF existing methods and challenges/data/split by population/")
#--------- liza layer1 ---------#
c1<-read.FCS("HuImmProfiling_S1_PBMC_1_07_CD19_CD27-____Naive_B_Cells.fcs",transformation = FALSE)
c2<-read.FCS("HuImmProfiling_S1_PBMC_1_08_CD19_CD27_____Total_Memory_B_Cells.fcs",transformation = FALSE)
c3<-read.FCS("HuImmProfiling_S1_PBMC_1_14_CD38_CD14hi____Classical_Monocytes.fcs",transformation = FALSE)
c4<-read.FCS("HuImmProfiling_S1_PBMC_1_15_CD38lo_-CD14int____Transitional__Monocytes.fcs",transformation = FALSE)
c5<-read.FCS("HuImmProfiling_S1_PBMC_1_16_CD38-CD14-____Non-Classical_Monocytes.fcs",transformation = FALSE)
c6<-read.FCS("HuImmProfiling_S1_PBMC_1_58_CD123_CD294_____Basophils.fcs",transformation = FALSE)
c7<-read.FCS("HuImmProfiling_S1_PBMC_1_20_CD56_CD57-____Early_NKs.fcs",transformation = FALSE)
c8<-read.FCS("HuImmProfiling_S1_PBMC_1_21_CD56_CD57_____Late_NKs.fcs",transformation = FALSE)
c9<-read.FCS("HuImmProfiling_S1_PBMC_1_23_CD123_CD11c-____pDC.fcs",transformation = FALSE)
c10<-read.FCS("HuImmProfiling_S1_PBMC_1_25_CD11c_CD38_____mDC.fcs",transformation = FALSE)
c11<-read.FCS("HuImmProfiling_S1_PBMC_1_32_CD45RA_CD45RO-____CD8_Naive.fcs",transformation = FALSE)
c12<-read.FCS("HuImmProfiling_S1_PBMC_1_33_CD45RA-CD45RO_____CD8_Central_Memory.fcs",transformation = FALSE)
c13<-read.FCS("HuImmProfiling_S1_PBMC_1_35_CD8_CD27_____CD8_Effector_Memory.fcs",transformation = FALSE)
c14<-read.FCS("HuImmProfiling_S1_PBMC_1_36_CD8_CD27-____CD8_Terminal_Effector.fcs",transformation = FALSE)
c15<-read.FCS("HuImmProfiling_S1_PBMC_1_39_CD45RA_CD45RO-____CD4_Naive_.fcs",transformation = FALSE)
c16<-read.FCS("HuImmProfiling_S1_PBMC_1_40_CD45RA-CD45RO_____CD4_Central_Memory_.fcs",transformation = FALSE)
c17<-read.FCS("HuImmProfiling_S1_PBMC_1_43_CD45RO_CD27_____CD4_Effector_Memory.fcs",transformation = FALSE)
c18<-read.FCS("HuImmProfiling_S1_PBMC_1_44_CD45RO_CD27-____CD4_Terminal_Effector.fcs",transformation = FALSE)
#c19<-read.FCS("HuImmProfiling_S1_PBMC_1_47_CD25hiCD127lo_-____Treg.fcs",transformation = FALSE)
#c20<-read.FCS("HuImmProfiling_S1_PBMC_1_51_CXCR3_CCR6-____Th1-like.fcs",transformation = FALSE)
#c21<-read.FCS("HuImmProfiling_S1_PBMC_1_53_CXCR3-CCR6-____Th2-like.fcs",transformation = FALSE)
#c22<-read.FCS("HuImmProfiling_S1_PBMC_1_54_CXCR3-CCR6_____Th17-like.fcs",transformation = FALSE)
c19<-read.FCS("HuImmProfiling_S1_PBMC_1_56_CD28_CD161hi____CD4-_MAIT_NKT.fcs",transformation = FALSE)
c20<-read.FCS("HuImmProfiling_S1_PBMC_1_60_CD3_TCR_______CD4-CD8-____T_Cells.fcs",transformation = FALSE)
#c25<-read.FCS("HuImmProfiling_S1_PBMC_1_03_CD294_CD16-____Eosinophils.fcs",transformation = FALSE)
#c26<-read.FCS("HuImmProfiling_S1_PBMC_1_04_CD294-CD16_____Neutrophils.fcs",transformation = FALSE)

nc=20
#exprs matrix
total.hul.m<-exprs(total.hul)
for(i in 1:nc){
  assign(paste0("c",i,".m"),exprs(get(paste0("c",i))))
}

#id for each cell type
total.id<-paste0(total.hul.m[,1],total.hul.m[,2])
for(i in 1:nc){
  assign(paste0("c",i,".id"),paste0(get(paste0("c",i,".m"))[,1],get(paste0("c",i,".m"))[,2]))
}

#overlap
names_cells<-paste0("c",1:nc)
overlap.matrix<-matrix(,nrow=nc,ncol=nc+2)
rownames(overlap.matrix)<-names_cells
colnames(overlap.matrix)<-c(names_cells,"total","overlap")
for(i in 1:nc){
  c1<-get(paste0(names_cells[i],".id"))
  overlap.matrix[i,nc+1]<-length(c1)
  overlap<-c()
  for(j in 1:nc){
    c2<-get(paste0(names_cells[j],".id"))
    overlap.matrix[i,j]<-length(intersect(c1,c2))
    if(i!=j) overlap<-unique(c(overlap,intersect(c1,c2)))
    if(j==nc) {
      overlap.matrix[i,nc+2]<-length(overlap)
      overlap.matrix[i,i]<-overlap.matrix[i,i]-length(overlap)
    }
  }
}

# identify overlap cell ids
com.id<-c()
for(i in 1:nc){
  com.id<-c(com.id,get(paste0("c",i,".id")))
}
overlap.id<-com.id[duplicated(com.id)]
overlap.id<-unique(overlap.id)
### creat true label ---------------
# 1. set overlap to 0 (we think those cells are ambiguous, do not use them to evaluate performance)
# 2. set cells do not belong to any type to 0
# 3. cells with 0 label were included in clustering, but not included in evaluation.
ref_label<-rep(0,length(total.id))
for(i in 1:nc){
  ref_label[total.id %in% get(paste0("c",i,".id"))]<-i
}
ref_label[total.id %in% overlap.id]<-0

#true marker table
mdf<-pData(parameters(total.hul))
markers<-c("89Y_CD45","163Dy_CD56_NCAM","168Er_CD14","144Nd_CD19","170Er_CD3","154Sm_CD27","161Dy_CD38","171Yb_CD20",
          "147Sm_CD11c","173Yb_HLA-DR","143Nd_CD123_IL-3R","166Er_CD294","150Nd_CD45RA","155Gd_CD57",
           "164Dy_TCRgd","145Nd_CD4","146Nd_CD8a","160Gd_CD28","151Eu_CD161","167Er_CD197_CCR7","149Sm_CD45RO")
attributes(ref_label)<-list(annotation=data.frame(index=1:20,
                                                  cell_type=c(
                                                      "naive b cells","memory b cells","classical monocytes","transitional monocytes","non-classical monocytes",
                                                      "basophils","early nks","late nks","pdc","mdc",
                                                      "cd8 naive","cd8 central memory","cd8 effector memory","cd8 terminal effector","cd4 naive",
                                                      "cd4 central memory","cd4 effector memory","cd4 terminal memory","cd4- mait/nkt","cd4-cd8-gd t cells")),
                            markers=mdf[mdf$desc %in% markers,])
names(ref_label)<-total.id
ref_label
save(ref_label,file="D:/research/Liza/cytof/CyTOF existing methods and challenges/data/ref_label.Rdata")
