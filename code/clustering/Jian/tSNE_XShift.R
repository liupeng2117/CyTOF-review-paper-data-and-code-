load("~/Methodology/201912_CytofClusterEvaluation/20200107_v1B/tsne.RData")
tsne_complete = tsne
load("~/Methodology/201912_CytofClusterEvaluation/20200107_v1B/cytof_tsne_downsample.RData")
tsne_downsample = tsne
library(ggplot2)

# X-shift
VorteX_26_Apr_2018 <- read.csv("~/Methodology/201912_CytofClusterEvaluation/20191216_v1A/VorteX.26-Apr-2018.csv")
Vortex_label = VorteX_26_Apr_2018[order(VorteX_26_Apr_2018$Index.in.File), 1]
cluster_label = factor(as.numeric(factor(Vortex_label)))
dataset = data.frame(group = cluster_label,
                     tsne=tsne_complete$Y)
pdf("X-shift_tsne.pdf", width = 10, height = 10, onefile = FALSE)
ggplot(dataset) + geom_point(aes(tsne.1, tsne.2, colour = group), size = 0.5) +
  labs(x = "tsne1",
       y = "tsne2") +
  guides(colour = guide_legend(override.aes = list(size=2.5))) +
  theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=3.5)))+
  theme(legend.text.align = 0,
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text=element_text(size=25),
        axis.title=element_text(size=25))
dev.off()