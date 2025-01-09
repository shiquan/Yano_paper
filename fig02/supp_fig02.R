require(Yano)

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf")

obj.sm <- readRDS("object_sm.rds")
obj.10x <- readRDS("object_10X.rds")

####################################################################################################
##  Supplementary Figure 2A
####################################################################################################
Idents(obj.sm) <- obj.sm$subclass_label

p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Caly", cell.tag = "RG", umi.tag = NULL, cell.group = obj.sm$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(6,1), highlights=c(140082462,140082548))
ggsave(p, file = "Caly_smartseq2.png", width=21, height = 7.6)

####################################################################################################
##  Supplementary Figure 2B
####################################################################################################
Idents(obj.10x) <- obj.10x$subclass_label
sl <- split(Idents(obj.10x),obj.10x$library_id)
sl <- lapply(sl, function(x) {
  nn <- gsub("(.*)_L.*", "\\1-1", names(x))
  names(x) <- nn
  x
  })

bamfiles <- list(
  "L8TX_181211_01_A01" = "10xV3/L8TX_181211_01_A01_CellRangerV5/L8TX_181211_01_A01_CellRangerV5.bam",
  "L8TX_190430_01_C08" = "10xV3/L8TX_190430_01_C08_CellRangerV5/L8TX_190430_01_C08_CellRangerV5.bam",
  "L8TX_190430_01_G08" = "10xV3/L8TX_190430_01_G08_CellRangerV5/L8TX_190430_01_G08_CellRangerV5.bam",
  "L8TX_181211_01_B01" = "10xV3/L8TX_181211_01_B01_CellRangerV5/L8TX_181211_01_B01_CellRangerV5.bam",
  "L8TX_190430_01_D08" = "10xV3/L8TX_190430_01_D08_CellRangerV5/L8TX_190430_01_D08_CellRangerV5.bam",
  "L8TX_190430_01_H08" = "10xV3/L8TX_190430_01_H08_CellRangerV5/L8TX_190430_01_H08_CellRangerV5.bam",
  "L8TX_190430_01_A08" = "10xV3/L8TX_190430_01_A08_CellRangerV5/L8TX_190430_01_A08_CellRangerV5.bam",
  "L8TX_190430_01_E08" = "10xV3/L8TX_190430_01_E08_CellRangerV5/L8TX_190430_01_E08_CellRangerV5.bam",
  "L8TX_190430_01_B08" = "10xV3/L8TX_190430_01_B08_CellRangerV5/L8TX_190430_01_B08_CellRangerV5.bam",
  "L8TX_190430_01_F08" = "10xV3/L8TX_190430_01_F08_CellRangerV5/L8TX_190430_01_F08_CellRangerV5.bam")

p <- TrackPlot(bamfile = bamfiles, gtf = gtf, gene = "Caly", cell.group = sl, layout.heights = c(6,1), junc = TRUE,  highlights=c(140082462,140082548))
ggsave(p, file = "Caly_10xV3.png", width=21, height = 7.6)

####################################################################################################
##  Supplementary Figure 2B and 2C
####################################################################################################
p1 <- FbtPlot(obj.10x, assay = c("exon", "junction"), shape.by="assay", col.by = "assay", val = "gene_name.padj", cols = c("red", "blue"), pt.size=2, gene = "Caly", gtf=gtf, layout.heights = c(6,1))
p2 <- FbtPlot(obj.sm, assay = c("exon", "junction"), shape.by="assay", col.by = "assay", val = "gene_name.padj", cols = c("red", "blue"), pt.size=2, gene = "Caly", gtf=gtf, point.label = c("chr7:140082462-140082508/Caly", "chr7:140082462-140082548/Caly", "chr7:140082093-140082227/Caly", "chr7:140082093-140082268/Caly", "chr7:140074100-140082462/Caly", "chr7:140074100-140082093/Caly"), label.size=5, layout.heights = c(6,1))
p <- cowplot::plot_grid(p1,p2, ncol=2)
ggsave(p, file = "Caly_zoom_in.pdf")

####################################################################################################
##  Supplementary Figure 2D
####################################################################################################
p <- FeaturePlot_scCustom(obj.sm, features = c("chr7:140082462-140082508/Caly", "chr7:140082093-140082227/Caly", "chr7:140082093-140082268/Caly", "chr7:140074100-140082462/Caly", "chr7:140074100-140082093/Caly"), pt.size=1, figure_plot = TRUE, num_columns = 5) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p, file = "Caly_features.pdf")
