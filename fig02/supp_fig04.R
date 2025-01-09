require(Yano)

obj.sm <- readRDS("object_sm.rds")
obj.10x <- readRDS("object_10X.rds")

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf")

############################################################################################################
## Supplementary Figure 4A
############################################################################################################
p1 <- FbtPlot(obj.sm, assay = c("exon", "junction"), shape.by = "assay", col.by = "assay", gtf = gtf, gene = "Stk16", val = "gene_name.padj", cols = c("red", "blue"), pt.size=2)
p2 <- FbtPlot(obj.10x, assay = c("exon", "junction"), shape.by = "assay", col.by = "assay", gtf = gtf, gene = "Stk16", val = "gene_name.padj", cols = c("red", "blue"), pt.size=2)
p <- cowplot::plot_grid(p1,p2,ncol=2)
ggsave(p, file = "Stk16_fbt.png", width=18, height = 8)


############################################################################################################
## Supplementary Figure 4B
############################################################################################################

p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Stk16", cell.tag = "RG", umi.tag = NULL, cell.group = obj.sm$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(6,1),max.depth = 50)#, display.genes = c("Stk16"))
#p
ggsave(p, file = "Stk16_smartseq.png", width=20, height = 8)

############################################################################################################
## Supplementary Figure 4C
############################################################################################################

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

p <- TrackPlot(bamfile = bamfiles, gtf = gtf, gene = "Actb", cell.group = sl, layout.heights = c(6,1), junc = TRUE, highlights=c(142904067,142904248))
ggsave(p, file = "Actb_10xV3.png", width=18, height = 10)

p <- TrackPlot(bamfile = bamfiles, gtf = gtf, gene = "Stk16", cell.group = sl, layout.heights = c(6,1), junc = TRUE)
ggsave(p, file = "Stk16_10xV3.png", width=20, height = 8)
