require(Rtsne)
library(ggplot2)
library(clue)
library(flowCore) # for loading fcs files

setwd("D:/research/Liza/cytof/CyTOF existing methods and challenges/data")
### load CyTOF data and manual gating reference label ----
SampleName<-"HuImmProfiling_S1_PBMC_1"
dt<-read.FCS(paste0(SampleName,"_01_CD45_CD66b-__Lymphocytes__DCs__Monocytes_.fcs"),transformation = FALSE)
load("ref_label.Rdata")

### preprocessing ----- 
#Extract markers
markers<- attributes(ref_label)$markers$name
dt2<-dt[,markers]

#asinh transformation to normalize the data
asinhTrans=arcsinhTransform(a=0, b=0.2)
translist<-transformList(markers, asinhTrans)
dt2.transform<-transform(dt2, translist)

dt2.transform<-exprs(dt2.transform)
ids<-ref_label!=0 # cells labeled with 0 are not included in evaluation
true_label<-factor(ref_label[ids])

#tsne
#tsne<-Rtsne(dt2.transform, dim=2,perplexity=30)

#save(tsne,file="tsne.rdata")
load("tsne.RData")
true_types<-c("Naïve B cells","Memory B cells","Classical Monocytes","Transitional Monocytes",
            "Non-classical Monocytes","Basophils","Early NK cells","Late NK cells",
            "Plasmacytoid DCs","Myeloid DCs","Naïve CD8 T cells","Central Memory CD8 T cells",
            "Effector Memory CD8 T cells","Terminal Effector CD8 T cells","Naïve CD4 T cells",
            "Central Memory CD4 T cells","Effector Memory CD4 T cells",
            "Terminal Effector CD4 T cells","NKT/MAIT cells","gd T cells")
true_ann<-c("1 Naïve B cells","2 Memory B cells","3 Classical Monocytes","4 Transitional Monocytes",
            "5 Non-classical Monocytes","6 Basophils","7 Early NK cells","8 Late NK cells",
            "9 Plasmacytoid DCs","10 Myeloid DCs","11 Naïve CD8 T cells","12 Central Memory CD8 T cells",
            "13 Effector Memory CD8 T cells","14 Terminal Effector CD8 T cells","15 Naïve CD4 T cells",
            "16 Central Memory CD4 T cells","17 Effector Memory CD4 T cells",
            "18 Terminal Effector CD4 T cells","19 NKT/MAIT cells",expression("20 "*gamma*delta*" T cells"))
dataset = data.frame(populations =true_label, 
                     tsne=tsne$Y[ids,])
my_colors=c("#CC0033","#FF0033","#FF6633","#FF9900","#FFCC00","#FFFF00","#CCFFCC","#99FF33","#66CC00","#339900",
            "#336600","#003366","#0066CC","#66CCFF","#99FFFF","#9999FF","#9933FF","#660099","#FF33FF", "#000033")

pp<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = populations), size = 0.5) +
  scale_color_manual(name = "Manually gated populations", labels=true_ann,values=my_colors)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
library(grid)
library(gridExtra) 
legend <- cowplot::get_legend(pp)
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/true.pdf",width=12,height=7)
pp
dev.off()
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/legend.pdf",width=10,height=10)
grid.draw(legend)
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/true.png",width=900,height=600)
pp
dev.off()

true_medianByclusters<-t(sapply(1:20,function(x) apply(dt2.transform[ref_label==x,],2,median)))

#FlowSOM
load("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/cluster_label_flowSOM.R")
cluster_label<-as.factor(cluster_label$cluster)
flowsom_medianByclusters<-t(sapply(1:20,function(x) apply(dt2.transform[cluster_label==x,],2,median)))
dist_matrix<-as.matrix(dist(rbind(true_medianByclusters,flowsom_medianByclusters),upper=T))
dist_matrix<-dist_matrix[1:20,21:40]
assignment<-solve_LSAP(dist_matrix, maximum = FALSE)
assignment_matrix<-cbind(true=seq_along(assignment), flowsom=assignment)
#library(plyr)
#cluster_label_new<-factor(mapvalues(as.numeric(cluster_label),from=assignment_matrix[,"flowsom"],to=assignment_matrix[,"true"]))
#tsne
#tsne<-Rtsne(dt2.transform, dim=2,perplexity=30)
library(readxl)
annotations<-read_excel("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/Copy of CyTOF review manuscript clusters.xlsx",skip=2)
ann_flowsom<-as.data.frame(annotations[1:20,7:8])
#cluster_ann<-paste(1:20, ann_flowsom$ID...8)
#cluster_ann[3]<-expression("3 "*gamma*delta*" T cells")
#colors_selected<-my_colors[order(assignment_matrix[1:20,"flowsom"])]
#Find clusters annotated with the same colors
colors_selected<-rep("a",length(unique(ann_flowsom$ID...8)))
cluster_ann<-rep("a",length(unique(ann_flowsom$ID...8)))
cluster_label<-as.numeric(cluster_label)
i=1
for(ct in true_types){
  if(tolower(ct) %in% unique(tolower(ann_flowsom$ID...8))){
    colors_selected[i]<-my_colors[which(tolower(true_types)==tolower(ct))]
    ct.clusters<-which(tolower(ann_flowsom$ID...8) == tolower(ct))
    cluster_ann[i]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    cluster_label[cluster_label %in% ct.clusters]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    print(i)
    print(sum(is.na(cluster_label)))
    i=i+1
  }
}

colors_selected[16]<-"#999999"
cluster_ann[16]<- "DNT cells - 5"
cluster_label[cluster_label %in% "5"]<-"DNT cells - 5"
cluster_label<-factor(cluster_label,levels=cluster_ann)

cluster_ann2<-cluster_ann
cluster_ann2[15]<-expression(gamma*delta*" T cells - 3")
dataset = data.frame(clusters =cluster_label[ids], 
                     tsne=tsne$Y[ids,])
pp_flowsom<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Groups", labels=cluster_ann2, values=colors_selected)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/flowsom.pdf",width=12,height=7)
pp_flowsom
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/flowsom.png",width=900,height=600)
pp_flowsom
dev.off()

#SPADE
load("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/rphenograph_and_spade.RData")
cluster_label_spade<-as.factor(cluster_label_spade)
#spade_medianByclusters<-t(sapply(1:20,function(x) apply(dt2.transform[cluster_label_spade==x,],2,median)))
#dist_matrix<-as.matrix(dist(rbind(true_medianByclusters,spade_medianByclusters),upper=T))
#dist_matrix<-dist_matrix[1:20,21:40]
#assignment<-solve_LSAP(dist_matrix, maximum = FALSE)
#assignment_matrix<-cbind(true=seq_along(assignment), spade=assignment)
#library(plyr)
#cluster_label_spade_new<-factor(mapvalues(as.numeric(cluster_label_spade),from=assignment_matrix[,"spade"],to=assignment_matrix[,"true"]))
#tsne
#tsne<-Rtsne(dt2.transform, dim=2,perplexity=30)
ann_spade<-as.data.frame(annotations[1:20,1:2])
#cluster_ann_spade<-paste(1:20, ann_spade$ID...2)
#cluster_ann_spade[1]<-expression("1 "*gamma*delta*" T cells")
#colors_selected<-my_colors[order(assignment_matrix[1:20,"spade"])]
colors_selected_spade<-rep("a",length(unique(ann_spade$ID...2)))
cluster_ann_spade<-rep("a",length(unique(ann_spade$ID...2)))
cluster_label_spade<-as.numeric(cluster_label_spade)
i=1
for(ct in true_types){
  if(tolower(ct) %in% unique(tolower(ann_spade$ID...2))){
    colors_selected_spade[i]<-my_colors[which(tolower(true_types)==tolower(ct))]
    ct.clusters<-which(tolower(ann_spade$ID...2) == tolower(ct))
    cluster_ann_spade[i]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    cluster_label_spade[cluster_label_spade %in% ct.clusters]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    print(i)
    print(sum(is.na(cluster_label_spade)))
    i=i+1
  }
}
cluster_label_spade<-factor(cluster_label_spade,levels=cluster_ann_spade)
cluster_ann_spade2<-cluster_ann_spade
cluster_ann_spade2[18]<-expression(gamma*delta*" T cells - 1")

dataset = data.frame(clusters =cluster_label_spade[ids], 
                     tsne=tsne$Y[ids,])

pp_spade<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", labels=cluster_ann_spade2, values=colors_selected_spade)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/spade.pdf",width=12,height=7)
pp_spade
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/spade.png",width=900,height=600)
pp_spade
dev.off()

#phenograph
load("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/rphenograph_and_spade.RData")
cluster_label_Rphenograph<-as.factor(cluster_label_Rphenograph)
#nclusters<-length(levels(cluster_label_Rphenograph))
#Rphenograph_medianByclusters<-t(sapply(1:nclusters,function(x) apply(dt2.transform[cluster_label_Rphenograph==x,],2,median)))
#dist_matrix<-as.matrix(dist(rbind(true_medianByclusters,Rphenograph_medianByclusters),upper=T))
#dist_matrix<-dist_matrix[1:20,21:51]
#assignment<-solve_LSAP(dist_matrix, maximum = FALSE)
#assignment_matrix<-cbind(true=seq_along(assignment), Rphenograph=assignment)
#assignment_others<-cbind(true=21:31,Rphenograph=(1:31)[-c(assignment_matrix[,"Rphenograph"])])
#assignment_matrix<-rbind(assignment_matrix,assignment_others)
#library(plyr)
#cluster_label_Rphenograph_new<-factor(mapvalues(as.numeric(cluster_label_Rphenograph),from=assignment_matrix[,"Rphenograph"],to=assignment_matrix[,"true"]))
#tsne
#tsne<-Rtsne(dt2.transform, dim=2,perplexity=30)


ann_Rphenograph<-as.data.frame(annotations[1:31,3:4])
#cluster_ann_Rphenograph<-paste(1:31, ann_Rphenograph$ID...4)
#cluster_ann_Rphenograph[c(4,28)]<-c(expression("4 "*gamma*delta*" T cells"),expression("28 "*gamma*delta*" T cells"))
#color_data_Rphenograph<-unique(ggplot_build(pp_Rphenograph)$data[[1]][,c("colour","group")])
#color_data_Rphenograph<-color_data_Rphenograph[order(color_data_Rphenograph$group),]
#colors_selected<-color_data_Rphenograph$colour
#colors_selected[assignment_matrix[1:20,"Rphenograph"]]<-color_data$colour

colors_selected_Rphenograph<-rep("a",length(unique(ann_Rphenograph$ID...4)))
cluster_ann_Rphenograph<-rep("a",length(unique(ann_Rphenograph$ID...4)))
cluster_label_Rphenograph<-as.numeric(cluster_label_Rphenograph)
i=1
for(ct in true_types){
  if(tolower(ct) %in% unique(tolower(ann_Rphenograph$ID...4))){
    colors_selected_Rphenograph[i]<-my_colors[which(tolower(true_types)==tolower(ct))]
    ct.clusters<-which(tolower(ann_Rphenograph$ID...4) == tolower(ct))
    cluster_ann_Rphenograph[i]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    cluster_label_Rphenograph[cluster_label_Rphenograph %in% ct.clusters]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    print(i)
    print(sum(is.na(cluster_label_Rphenograph)))
    i=i+1
  }
}

colors_selected_Rphenograph[17]<-"#666666"
#setdiff(my_colors, colors_selected_Rphenograph)[1]
cluster_ann_Rphenograph[17]<- "NKT cells - 14"
cluster_label_Rphenograph[cluster_label_Rphenograph %in% "14"]<-"NKT cells - 14"

colors_selected_Rphenograph[18]<-"#CCCCCC"
#setdiff(my_colors, colors_selected_Rphenograph)[1]
cluster_ann_Rphenograph[18]<- "other - 3"
cluster_label_Rphenograph[cluster_label_Rphenograph %in% "3"]<-"other - 3"
cluster_label_Rphenograph<-factor(cluster_label_Rphenograph,levels=cluster_ann_Rphenograph)

cluster_ann_Rphenograph2<-cluster_ann_Rphenograph
cluster_ann_Rphenograph2[16]<-expression(gamma*delta*" T cells - 4,28")


dataset = data.frame(clusters =cluster_label_Rphenograph[ids], 
                     tsne=tsne$Y[ids,])
pp_Rphenograph<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", labels=cluster_ann_Rphenograph2, values=colors_selected_Rphenograph)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/Rphenograph.pdf",width=12,height=7)
pp_Rphenograph
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/Rphenograph.png",width=900,height=600)
pp_Rphenograph
dev.off()



#Import the downsampled data
setwd("D:/research/Liza/cytof/CyTOF existing methods and challenges/data")
load("ref_label_downsample.Rdata")
sample_id<- sapply(names(ref_label_downsample), function(x) which(names(ref_label)==x))
dt2.transform_downsample<-dt2.transform[sample_id,]
ids2=ref_label_downsample!=0
#DensVM
cluster_label_densvm<-read.csv("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/DensVM_label.csv")
cluster_label_densvm<-as.factor(cluster_label_densvm$cluster_label)
#nclusters<-length(levels(cluster_label_densvm))
#densvm_medianByclusters<-t(sapply(1:nclusters,function(x) apply(dt2.transform_downsample[cluster_label_densvm==x,],2,median)))
#dist_matrix<-as.matrix(dist(rbind(true_medianByclusters,densvm_medianByclusters),upper=T))
#dist_matrix<-dist_matrix[21:37,1:20]
#assignment<-solve_LSAP(dist_matrix, maximum = FALSE)
#assignment_matrix<-cbind(true=assignment, densvm=seq_along(assignment))
#assignment_others<-cbind(true=(1:20)[-c(assignment_matrix[,"true"])],densvm=18:20)
#assignment_matrix<-rbind(assignment_matrix,assignment_others)
#library(plyr)
#cluster_label_densvm_new<-factor(mapvalues(as.numeric(cluster_label_densvm),from=assignment_matrix[,"densvm"],to=assignment_matrix[,"true"]))
#tsne
#tsne<-Rtsne(dt2.transform, dim=2,perplexity=30)


ann_densvm<-as.data.frame(annotations[1:17,5:6])
#cluster_ann_densvm<-paste(1:17, ann_densvm$ID...6)
#cluster_ann_densvm[15]<-expression("15 "*gamma*delta*" T cells")
#colors_selected<-color_data$colour[assignment_matrix[1:17,"true"]]
colors_selected_densvm<-rep("a",length(unique(ann_densvm$ID...6)))
cluster_ann_densvm<-rep("a",length(unique(ann_densvm$ID...6)))
cluster_label_densvm<-as.numeric(cluster_label_densvm)
i=1
for(ct in true_types){
  if(tolower(ct) %in% unique(tolower(ann_densvm$ID...6))){
    colors_selected_densvm[i]<-my_colors[which(tolower(true_types)==tolower(ct))]
    ct.clusters<-which(tolower(ann_densvm$ID...6) == tolower(ct))
    cluster_ann_densvm[i]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    cluster_label_densvm[cluster_label_densvm %in% ct.clusters]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    print(i)
    print(sum(is.na(cluster_label_densvm)))
    i=i+1
  }
}
cluster_label_densvm<-factor(cluster_label_densvm,levels=cluster_ann_densvm)
cluster_ann_densvm2<-cluster_ann_densvm
cluster_ann_densvm2[15]<-expression(gamma*delta*" T cells - 15")

dataset = data.frame(clusters =cluster_label_densvm[ids2], 
                     tsne=tsne$Y[sample_id,][ids2,])
pp_densvm<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 1.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", labels=cluster_ann_densvm2,values=colors_selected_densvm)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))

pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/densvm.pdf",width=12,height=7)
pp_densvm
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/densvm.png",width=900,height=600)
pp_densvm
dev.off()


#ACCENSE
cluster_label_accense<-read.csv("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/accense_label.csv")
cluster_label_accense<-as.factor(cluster_label_accense$cluster_label)
#nclusters<-length(levels(cluster_label_accense))
#accense_medianByclusters<-t(sapply(1:nclusters,function(x) apply(dt2.transform_downsample[cluster_label_accense==x,],2,median)))
#dist_matrix<-as.matrix(dist(rbind(true_medianByclusters,accense_medianByclusters),upper=T))
#dist_matrix<-dist_matrix[1:20,21:76]
#assignment<-solve_LSAP(dist_matrix, maximum = FALSE)
#assignment_matrix<-cbind(true=seq_along(assignment), accense=assignment)
#assignment_others<-cbind(true=21:56,accense=(1:56)[-c(assignment_matrix[,"accense"])])
#assignment_matrix<-rbind(assignment_matrix,assignment_others)
#library(plyr)
#cluster_label_accense_new<-factor(mapvalues(as.numeric(cluster_label_accense),from=assignment_matrix[,"accense"],to=assignment_matrix[,"true"]))
#tsne
#tsne<-Rtsne(dt2.transform, dim=2,perplexity=30)




ann_accense<-as.data.frame(annotations[1:56,9:10])
#cluster_ann_accense<-paste(1:56, ann_accense$ID...10)
#cluster_ann_accense[c(1,53)]<-c(expression("1 "*gamma*delta*" T cells"),expression("53 "*gamma*delta*" T cells"))
#color_data_accense<-unique(ggplot_build(pp_accense)$data[[1]][,c("colour","group")])
#color_data_accense<-color_data_accense[order(color_data_accense$group),]
#colors_selected<-color_data_accense$colour
#colors_selected[assignment_matrix[1:20,"accense"]]<-color_data$colour
colors_selected_accense<-rep("a",length(unique(ann_accense$ID...10)))
cluster_ann_accense<-rep("a",length(unique(ann_accense$ID...10)))
cluster_label_accense<-as.numeric(cluster_label_accense)
i=1
for(ct in true_types){
  if(tolower(ct) %in% unique(tolower(ann_accense$ID...10))){
    colors_selected_accense[i]<-my_colors[which(tolower(true_types)==tolower(ct))]
    ct.clusters<-which(tolower(ann_accense$ID...10) == tolower(ct))
    cluster_ann_accense[i]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    cluster_label_accense[cluster_label_accense %in% ct.clusters]<-paste(ct,paste(ct.clusters,collapse = ","),sep = " - ")
    print(i)
    print(sum(is.na(cluster_label_accense)))
    i=i+1
  }
}
cluster_label_accense<-factor(cluster_label_accense,levels=cluster_ann_accense)
cluster_ann_accense2<-cluster_ann_accense
cluster_ann_accense2[19]<-expression(gamma*delta*" T cells - 1,53")

dataset = data.frame(clusters =cluster_label_accense[ids2], 
                     tsne=tsne$Y[sample_id,][ids2,])
pp_accense<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 1.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", labels=cluster_ann_accense2, values=colors_selected_accense)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))

pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/accense.pdf",width=15,height=7)
pp_accense
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/accense.png",width=1150,height=600)
pp_accense
dev.off()


