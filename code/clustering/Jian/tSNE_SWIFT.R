load("/Users/jianzou/Methodology/201912_CytofClusterEvaluation/Review/data/ref_label_downsample.Rdata")
load("/Users/jianzou/Methodology/201912_CytofClusterEvaluation/Review/data/ref_label.Rdata")
cell_selected = rownames(data.frame(ref_label_downsample))
cell_select = which(rownames(data.frame(ref_label)) %in% cell_selected)

library(ggplot2)
source("SWIFT_downsampled.R")
load("~/Methodology/201912_CytofClusterEvaluation/20200107_v1B/tsne.RData")
tsne_complete = tsne
label = as.factor(label)
dataset = data.frame(group = cluster_label,
                     tsne=tsne_complete$Y[cell_select,])
pdf("SWIFT_tsne.pdf", width = 10, height = 10, onefile = FALSE)
ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = group), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  guides(colour = guide_legend(override.aes = list(size=2.5))) +
  theme_classic() + 
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(size=3.5)))+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text=element_text(size=25),
        axis.title=element_text(size=25))
dev.off()