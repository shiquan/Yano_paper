require(Yano)

obj <- readRDS("AHCA.rds")

########################################################################################################################################
##                  Supplementary Figure 8A
########################################################################################################################################

## gencode.v38.annotation.gtf.gz download from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gtf <- gtf2db("./gencode.v38.annotation.gtf.gz")
p <- FbtPlot(obj, val = "locus.padj", assay = "var", col.by = "GnomAD", shape.by = "strand", cols = c("red", "blue"), gtf = gtf, chr = "chr2", start = 88500000, end = 90500000, pt.size=2)
ggsave(p, file = "IGKVs_fbt.pdf")

########################################################################################################################################
##                  Supplementary Figure 8B
########################################################################################################################################
p <- FbtPlot(obj, val = "locus.padj", assay = "var", col.by = "GnomAD", shape.by = "strand", cols = c("red", "blue"), gtf = gtf, chr = "chr11", start = 290000, end = 330000, pt.size=2)
ggsave(p, file = "IFITMs_fbt.pdf")

########################################################################################################################################
##                  Supplementary Figure 8C
########################################################################################################################################
FeaturePlot(obj, features = c("chrX:119468424C=/+",  "chrX:119468424C>G/+"), blend = TRUE, reduction = "tsne", pt.size=1, order = TRUE)
