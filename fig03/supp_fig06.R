require(Yano)
require(RColorBrewer)

obj.sm <- readRDS("object_sm.rds")
obj.sm1 <- readRDS("object_sm1.rds")

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf")

###################################################################################################
##   Supplementary Figure 6A
###################################################################################################

p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Dock4", cell.tag = "RG", umi.tag = NULL, cell.group = obj.sm$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(3,1), highlights=c(40799657,40799750))
ggsave(p, file = "Dock4_Smartseq2.png", width = 15, height = 5)

###################################################################################################
##   Supplementary Figure 5B
###################################################################################################
p1 <- FeaturePlot_scCustom(obj.sm, features = c("chr12:40799657-40799750/Dock4"), order = TRUE, pt.size=1, figure_plot=TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p2 <- FeaturePlot_scCustom(obj.sm, features = c("Dock4"), order = TRUE, pt.size=1, figure_plot=TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p <- cowplot::plot_grid(p1,p2)
ggsave(p, file = "Dock4_feature.pdf")

###################################################################################################
##   Supplementary Figure 5C
###################################################################################################
p1 <- FbtPlot(obj.sm1, assay = c("exon", "junction"), val = "gene_name.padj", gene = "Dock4", gtf=gtf, shape.by = "assay", col.by = "assay", pt.size=2, cols = c("blue", "red"), layout.heights = c(2,1))
p2 <- FbtPlot(obj.sm, assay = c("exon", "junction"), val = "gene_name.padj", gene = "Dock4", gtf=gtf, shape.by = "assay", col.by = "assay", pt.size=2, cols = c("blue", "red"), layout.heights = c(2,1))

p <- cowplot::plot_grid(p1,p2, ncol=2)
ggsave(p, file = "Dock4_fbt.pdf")

