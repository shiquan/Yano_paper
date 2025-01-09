require(Yano)
require(RColorBrewer)

obj.sm <- readRDS("object_sm.rds")
obj.sm1 <- readRDS("object_sm1.rds")

###################################################################################################
##   Supplementary Figure 5A
###################################################################################################

p <- FbtPlot(obj.sm1, assay = c("exon", "exclude"), shape.by="assay", col.by = "assay", remove.chr = TRUE, val = "exon_name.padj")
ggsave(p, file = "sm3_exon_junc_mode1.png", width = 15, height = 6)


###################################################################################################
##   Supplementary Figure 5B
###################################################################################################
obj.sm[['exon']][[]] -> df1
obj.sm1[['exon']][[]] -> df2

nm <- intersect(rownames(df1), rownames(df2))
df1 <- df1[nm,]
df2 <- df2[nm,]
df1[['sel.padj']]<- df2$gene_name.padj
p <- ggplot() + geom_point(data=df1,aes(-log10(gene_name.padj), -log10(sel.padj)))  + geom_label_repel(data=df1[c("chr9:78478449-78478858/Eef1a1"),], aes(-log10(gene_name.padj), -log10(sel.padj), label = exon_name), nudge_y = 10, nudge_x=10, size=5) + theme_pubr() + xlab("") + ylab("")
ggsave(p, file = "sel_vs_all.png")

###################################################################################################
##   Supplementary Figure 5C
###################################################################################################

p1 <- FeaturePlot(obj.sm1, features = c("chr9:78478449-78478858/Eef1a1", "Eef1a1"), slot = "counts") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p2 <- FeaturePlot(obj.sm, features = c("chr9:78478449-78478858/Eef1a1", "Eef1a1"),  slot = "counts") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p <- cowplot::plot_grid(p1,p2,ncol=2)
ggsave(p, file = "Eef1a1_features.pdf")


###################################################################################################
##   Supplementary Figure 5D
###################################################################################################

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf")
p <- TrackPlot(bamfile = "./Smartseq2/merged.bam", gtf = gtf, gene = "Ly6e", cell.tag = "RG", umi.tag = NULL, cell.group = obj.sm1$subclass_label, junc=TRUE, strand = "ignore", layout.heights = c(6,1), highlights=c(74955071,74955717))
ggsave(p, file = "Ly6e_track.png", width=18, height = 8)
