require(Yano)

obj <- readRDS("AHCA.rds")


########################################################################################################################################
##                  Supplementary Figure 10A
########################################################################################################################################
## Plotted at supp_fig09.R

########################################################################################################################################
##                  Supplementary Figure 10B
########################################################################################################################################

p <- FeaturePlot(obj, features = c("chr6:80273134G>A/+"), order = TRUE, pt.size=1, coord.fixed = TRUE) + theme_void() + ggtitle("")
ggsave(p, file = "chr6_80273134G.png", width = 5,height = 5)

p <- FeaturePlot(obj, features = c("BCKDHB"), order = TRUE, pt.size=1, coord.fixed = TRUE) + theme_void() + ggtitle("")
ggsave(p, file = "BCKDHB.png", width = 5,height = 5)

########################################################################################################################################
##                  Supplementary Figure 10C
########################################################################################################################################

p <- FeaturePlot(obj, features = c("chr22:20709921TTGTC>T/-"), order = TRUE, pt.size=1, coord.fixed = TRUE) + theme_void() + ggtitle("")
ggsave(p, file = "chr22_20709921.png", width = 5,height = 5)

p <- FeaturePlot(obj, features = c("PI4KA"), order = TRUE, pt.size=1, coord.fixed = TRUE) + theme_void() + ggtitle("")
ggsave(p, file = "PI4KA.png", width = 5,height = 5)
