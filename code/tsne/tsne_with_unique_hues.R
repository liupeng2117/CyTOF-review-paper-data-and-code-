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
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne/true.pdf",width=10,height=7)
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
#library(readxl)
#annotations<-read_excel("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/Copy of CyTOF review manuscript clusters.xlsx",skip=2)
#ann_flowsom<-as.data.frame(annotations[1:20,7:8])

dataset = data.frame(clusters =cluster_label[ids], 
                     tsne=tsne$Y[ids,])
pp_flowsom<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Groups", values=my_colors)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/flowsom.pdf",width=9,height=7)
pp_flowsom
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/flowsom.png",width=800,height=600)
pp_flowsom
dev.off()

#SPADE
load("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/rphenograph_and_spade.RData")
cluster_label_spade<-as.factor(cluster_label_spade)

dataset = data.frame(clusters =cluster_label_spade[ids], 
                     tsne=tsne$Y[ids,])

pp_spade<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", values=my_colors)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/spade.pdf",width=9,height=7)
pp_spade
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/spade.png",width=800,height=600)
pp_spade
dev.off()

#phenograph
load("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/rphenograph_and_spade.RData")
cluster_label_Rphenograph<-as.factor(cluster_label_Rphenograph)

dataset = data.frame(clusters =cluster_label_Rphenograph, 
                     tsne=tsne$Y)
my_colors_phenograph<-c("#CC0033","#FF0033","#FF6633","#FF9900","#FFCC00","#FFFF00","#FFFF99","#CCCC66","#666633","#999933","#CCFF99","#CCFFCC","#CCFF99","#99FF33","#66FF66","#66CC00","#339900",
                          "#336600","#003366","#0066CC","#0000FF","#66CCFF","#99FFFF","#9999FF","#9933FF","#660099","#FF33FF","#CC33FF","#9900FF", "#000033","#CCCCCC")
pp_Rphenograph<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", values=my_colors_phenograph)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/Rphenograph.pdf",width=9,height=7)
pp_Rphenograph
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/Rphenograph.png",width=800,height=600)
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

dataset = data.frame(clusters =cluster_label_densvm[ids2], 
                     tsne=tsne$Y[sample_id,][ids2,])
my_colors_densvm<-c("#CC0033","#FF6633","#FF9900","#FFCC00","#FFFF00","#CCFFCC","#99FF33","#66CC00","#339900",
                    "#336600","#003366","#0066CC","#66CCFF","#99FFFF","#9999FF","#9933FF","#FF33FF")

pp_densvm<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 1.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", values=my_colors_densvm)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))

pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/densvm.pdf",width=9,height=7)
pp_densvm
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/densvm.png",width=800,height=600)
pp_densvm
dev.off()


#ACCENSE
cluster_label_accense<-read.csv("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/tsne_label/accense_label.csv")
cluster_label_accense<-as.factor(cluster_label_accense$cluster_label)
dataset = data.frame(clusters =cluster_label_accense[ids2], 
                     tsne=tsne$Y[sample_id,][ids2,])
my_colors_accense<-c("#CC0033","#FF0033","#FF0099","#FF66CC","#FF3333","#FF0000","#FF6633","#FF9900","#FF9999","#FFCC00","#FFFF00","#FFFF99","#CCFF99","#CCCC66","#666633","#999933","#669933","#CC9933","#CCFF99","#CCFFCC","#CCFF99","#33FF33","#99FF33","#66FF66","#33CC33","#66CC00","#339900","#66CC33",
                     "#003366","#336600","#66FFFF","#00CCCC","#0099CC","#0066CC","#006699","#3366CC","#0066FF","#0000FF","#000099","#003399","#0033FF","#66CCFF","#66FFFF","#99FFFF","#9999FF","#9933FF","#993399","#663399","#660099","#FF33FF","#FF33CC","#CC33FF","#9900FF","#FF0033", "#000033","#CCCCCC")

pp_accense<-ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = clusters), size = 1.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  scale_color_manual(name="Clusters", values=my_colors_accense)+
  guides(colour = guide_legend(override.aes = list(size=2.5)))+
  theme_classic()+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))

pdf("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/accense.pdf",width=10,height=7)
pp_accense
dev.off()
png("D:/research/Liza/cytof/CyTOF existing methods and challenges/figures and tables/tsne/five methods tsne with unique hues/accense.png",width=900,height=600)
pp_accense
dev.off()
