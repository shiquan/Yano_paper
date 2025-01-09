require(Yano)
require(RColorBrewer)

obj.sm <- readRDS("object_sm.rds")
obj.sm1 <- readRDS("object_sm1.rds")

###################################################################################################
##   Supplementary Figure 7A
###################################################################################################
p <- FbtPlot(obj.sm1, assay = c("exon", "exclude"), val = "exon_name.padj", remove.chr = TRUE, sel.chrs = c(1:21,"X"), shape.by = "assay", col.by = "assay", pt.size=2, cols = c("grey40", "yellow"), point.label = "chr18:38447713-38447821/Ndfip1", label.size=8) + xlab("") + ylab("") + theme(legend.position = "none")
ggsave(p, file = "sm3_exon_junc_mode1.png", width = 15, height = 6)


###################################################################################################
##   Supplementary Figure 7B, 7C, 7D and other two genes
###################################################################################################

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf")


p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Clta", cell.tag = "RG", umi.tag = NULL, cell.group = obj$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(3,1), highlights=list(c(44030215,44030268),c(44031400,44031435)))
ggsave(p, file = "Clta_Smartseq2.png", width = 15, height = 5)


p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Gria3", cell.tag = "RG", umi.tag = NULL, cell.group = obj$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(3,1), highlights=list(c(41669501,41669615),c(41654810,41654924)))
ggsave(p, file = "Gria3_Smartseq2.png", width = 15, height = 5)


p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Uap1", cell.tag = "RG", umi.tag = NULL, cell.group = obj$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(3,1), highlights=list(c(170147996,170148046),c(170147996,170148043)))
ggsave(p, file = "Uap1_Smartseq2.png", width = 15, height = 5)

p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Fam49b", cell.tag = "RG", umi.tag = NULL, cell.group = obj$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(3,1), highlights=list(c(63956600,63956682),c(63957585,63957667)))
ggsave(p, file = "Fam49b_Smartseq2.png", width = 15, height = 5)

p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Ttc7b", cell.tag = "RG", umi.tag = NULL, cell.group = obj$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(3,1), highlights=c(100340497,100340547))
ggsave(p, file = "Ttc7b_Smartseq2.png", width = 15, height = 5)

