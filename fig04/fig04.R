require(Yano)
require(data.table)

type.values <- c("#131313","blue",RColorBrewer::brewer.pal(12, "Paired"))
values <- c("#131313","blue",RColorBrewer::brewer.pal(12, "Paired"))
names(type.values) <- c("antisense_complex", "antisense_exon", "antisense_intron", "antisense_utr3", "antisense_utr5", "exon","exonintron", "intergenic", "intron", "multiexons", "multigenes","utr3", "utr5", "whole_gene")
sel.chrs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "MT")
sel.chrs <- factor(sel.chrs, levels = sel.chrs)

##########################################################################################################
## Construct Seurat Object
##########################################################################################################
## meta data download from 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144236/suppl/GSE144236%5Fpatient%5Fmetadata%5Fnew%2Etxt%2Egz'
meta <- fread("./GSE144236_patient_metadata_new.txt.gz")

cells <- meta$V1
#unique(gsub("(.*)_[ACGTN]+","\\1", cells))
rename_counts <- function(dir = NULL, name = NULL, cells = NULL, min.features = 0) {
  cnt <- ReadPISA(dir)
  nm <- paste0(name, "_", colnames(cnt))
  nm <- gsub("(.*)-1","\\1",nm)
  colnames(cnt) <- nm
  
  if (!is.null(cells)) {
    cells0 <- intersect(cells, colnames(cnt))
    cnt <- cnt[,cells0]
  }

  if (min.features > 0) {
    cs <- colSums(cnt > 0)
    idx <- which(cs >= min.features)
    cnt <- cnt[,idx]
  }
  
  cnt
}

exp18 <- rename_counts("./P8_cSCC_scRNA_1/exp", "P8_Tumor", cells=cells, min.features=500)
exp19 <- rename_counts("./P8_cSCC_scRNA_2/exp", "P8_Tumor", cells=cells, min.features=500)
idx1 <- setdiff(colnames(exp18),colnames(exp19))
idx2 <- setdiff(colnames(exp19),colnames(exp18))
exp18 <- exp18[,idx1]
exp19 <- exp19[,idx2]
p8.cSCC.cells.1 <- colnames(exp18)
p8.cSCC.cells.2 <- colnames(exp19)

exp20 <- rename_counts("./P8_normal_scRNA_1/exp", "P8_Normal", cells=cells, min.features=500)
exp21 <- rename_counts("./P8_normal_scRNA_2/exp", "P8_Normal", cells=cells, min.features=500)
idx1 <- setdiff(colnames(exp20),colnames(exp21))
idx2 <- setdiff(colnames(exp21),colnames(exp20))
exp20 <- exp20[,idx1]
exp21 <- exp21[,idx2]
p8.normal.cells.1 <- colnames(exp20)
p8.normal.cells.2 <- colnames(exp21)

exp <- mergeMatrix(exp18, exp19, exp20, exp21)

obj <- QuickRecipe(exp, min.cells=10)
meta <- as.data.frame(meta)
rownames(meta) <- meta$V1
meta <- meta[colnames(obj),]

obj$level1_celltype <- meta$level1_celltype
obj$level2_celltype <- meta$level2_celltype
obj$level3_celltype <- meta$level3_celltype

obj$sample <- gsub("(.*)_[ACGT]+", "\\1", colnames(obj))

##############################################################################################
##  Figure 4A and 4B
##############################################################################################
values0 <- c("#D9D9D9", "blue", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928")[1:12]
names(values0) <- c("Keratinocyte", "Fibroblast", "CD1C", "Tcell", "Mac", "LC", "CLEC9A", "B Cell", "NK", "MDSC", "PDC", "Pilosebaceous")
values1 <-  c("#000000", "#FFFFFF", "#000000", "#000000", "#000000", "#000000", "#000000", "#FFFFFF", "#000000", "#000000", "#000000", "#FFFFFF")
names(values1) <- c("Keratinocyte", "Fibroblast", "CD1C", "Tcell", "Mac", "LC", "CLEC9A", "B Cell", "NK", "MDSC", "PDC", "Pilosebaceous")

p1 <- DimPlot(obj, pt.size=2, group.by = "level2_celltype", label = TRUE, label.size=8, cols = values0, label.box = TRUE, label.color = values1)+ theme(legend.position = "none")
p2 <- DimPlot(obj, pt.size=2, group.by = "sample") + theme(legend.position = "none")
p <- cowplot::plot_grid(p1,p2)
ggsave(p, file = "cSCC_cluster.pdf")
ggsave(p, file = "cluster.pdf")


##############################################################################################
##  Figure 4C
##############################################################################################
cells <- colnames(obj)
var18 <- rename_counts("./P8_cSCC_scRNA_1/var", "P8_Tumor", cells=p8.cSCC.cells.1)
var19 <- rename_counts("./P8_cSCC_scRNA_2/var", "P8_Tumor", cells=p8.cSCC.cells.2)
var20 <- rename_counts("./P8_normal_scRNA_1/var", "P8_Normal", cells=p8.normal.cells.1)
var21 <- rename_counts("./P8_normal_scRNA_2/var", "P8_Normal", cells=p8.normal.cells.2)
p8.var <- mergeMatrix(var18, var19, var20, var21)

obj[['var']] <- CreateAssayObject(p8.var[, colnames(obj)], min.cells=10)
DefaultAssay(obj) <- 'var'
obj <- NormalizeData(obj)
obj <- RunAutoCorr(obj)
obj <- SetAutoCorrFeatures(obj)

gtf <- gtf2db("Homo_sapiens.GRCh37.87.chr.gtf.gz")
obj <- annoVAR(obj, gtf = gtf)

# No binding assay is specified, so test features which have the same binding feature with be aggreated together as the binding feature
# Use mode 2 indicate that test the EAT with other EATs at the same locus
obj <- RunBlockCorr(obj, bind.name = "locus", mode =2)

obj[['var']][[]] -> df
obj[['var']][['chr0']] <- paste0("chr", df$chr)

# Download from https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.25.gz
# The orginal vcf from NCBI is not well formated, because the allele frequenscy tags collapsed together within FREQ tag, make it hard to prase.
# I manually extract allele frequency and create tags for different databases.
# The VCF header of GCF_000001405.25_Freq.vcf.gz can be found at header.txt.
obj <- annoVAR(obj, vcf = "/home/projects/ku_00009/data/TranscriptScar/hg19/var/GCF_000001405.25_Freq.vcf.gz", tags = c("GnomAD", "dbGaP_PopFreq", "KOREAN", "TOMMO"), chr="chr0")

# Make sure all the record is numeric properly.
obj[['var']][[]] -> df
obj[['var']][['GnomAD']] <- as.numeric(df$GnomAD)
obj[['var']][['dbGaP']] <- as.numeric(df$dbGaP_PopFreq)
obj[['var']][['KOREAN']] <- as.numeric(df$GnomAD)
obj[['var']][['TOMMO']] <- as.numeric(df$TOMMO)

## Clinvar database is download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2024/clinvar_20240312.vcf.gz
obj <- annoVAR(obj, vcf = "/home/projects/ku_00009/data/TranscriptScar/hg19/clinvar_20240312.vcf.gz", tags = c("CLNSIG", "RS"), check.alt.only = TRUE)

wes <- ReadPISA("WES_var/")
wes1 <- as.data.frame(wes)
wes1 %>% filter(P8_normal_WES > 0) %>% rownames -> norm.var
wes1 %>% filter(P8_cSCC_WES > 0) %>% rownames -> tumor.var
wes1 %>% filter(P8_cSCC_WES > 0 & P8_normal_WES==0) %>% rownames -> tumor.only.var

df <- obj[['var']][[]]
df$var <- gsub("(.*)/[-+]","\\1", rownames(df))

df$anno <- "others"
idx <- match(norm.var,df$var)
idx <- idx[!is.na(idx)]
df[idx,"anno"] <- "normal"
idx <- match(tumor.only.var,df$var)
idx <- idx[!is.na(idx)]
df[idx,"anno"] <- "tumor.only"
obj[['var']][['var.type']] <- df$anno

#######################################################################################################################################
##  Figure 4C
#######################################################################################################################################
demo.sel <- c("17:73894804G=/-", "17:73894804G>A/-", "MT:11152T=/+", "MT:11152T>C/+")
#demo.sel <- c("11:252818G=/+", "11:252818G>A/+", "MT:15274A=/+", "MT:15274A>C/+")

p <- FbtPlot(obj, val = "locus.pval", col.by = "dbGaP", pt.size = 2, point.label = demo.sel, label.size=5, sel.chrs = sel.chrs, shape.by="var.type", xlab="", ylab="") + theme(legend.position = "none")
ggsave(p, file = "var.png", width = 15, height = 6)
ggsave(cowplot::get_legend(p), file = "var_legend.pdf")


#######################################################################################################################################
##  Figure 4D
#######################################################################################################################################
p <- FeaturePlot(obj, features = demo.sel, pt.size = 2, order = TRUE, coord.fixe= TRUE) & scale_colour_gradientn(colours = c("grey",brewer.pal(n = 8, name = "Reds")))


#######################################################################################################################################
##  Figure 4E, 4F use IGV to visualize specific regions and snapshot mannually
#######################################################################################################################################


#######################################################################################################################################
##  Figure 4G
#######################################################################################################################################
require(ggvenn)
p <- ggvenn(list(norm=norm.var, tumor=tumor.var), stroke_color="white")
ggsave(p, file = "wes_var.pdf", width=6, height = 5)


#######################################################################################################################################
##  Figure 4H
#######################################################################################################################################

norm.sl <- .Call("parse_var_names", norm.var)
tumor.only.sl <- .Call("parse_var_names", tumor.only.var)

norm.df <- varanno(chr=paste0("chr",norm.sl[[1]]), start = norm.sl[[2]], ref = norm.sl[[3]], alt = norm.sl[[4]], vcf = "/home/projects/ku_00009/data/LiverZonation/hg19/var/GCF_000001405.25_Freq.vcf.gz", tags = c("GnomAD", "dbGaP_PopFreq", "KOREAN", "TOMMO", "RS"))

tumor.only.df <- varanno(chr=paste0("chr",tumor.only.sl[[1]]), start = tumor.only.sl[[2]], ref = tumor.only.sl[[3]], alt = tumor.only.sl[[4]], vcf = "/home/projects/ku_00009/data/LiverZonation/hg19/var/GCF_000001405.25_Freq.vcf.gz", tags = c("GnomAD", "dbGaP_PopFreq", "KOREAN", "TOMMO", "RS"))


require(ggpubr)

norm.df$lowFreq <- "lowFreq"
idx <- which(norm.df$dbGaP_PopFreq>1e-4)
norm.df[idx,]$lowFreq <- "HighFreq"
norm.df1 <- as.data.frame(table(norm.df$lowFreq))

norm.df$lowFreq <- "lowFreq"
idx <- which(norm.df$GnomAD>1e-4)
norm.df[idx,]$lowFreq <- "HighFreq"
norm.df2 <- as.data.frame(table(norm.df$lowFreq))

norm.df$lowFreq <- "lowFreq"
idx <- which(norm.df$KOREAN>1e-4)
norm.df[idx,]$lowFreq <- "HighFreq"
norm.df3 <- as.data.frame(table(norm.df$lowFreq))

norm.df$lowFreq <- "lowFreq"
idx <- which(norm.df$TOMMO>1e-4)
norm.df[idx,]$lowFreq <- "HighFreq"
norm.df4 <- as.data.frame(table(norm.df$lowFreq))


tumor.only.df$lowFreq <- "lowFreq"
idx <- which(tumor.only.df$dbGaP_PopFreq>1e-4)
tumor.only.df[idx,]$lowFreq <- "HighFreq"
tumor.df1 <- as.data.frame(table(tumor.only.df$lowFreq))

tumor.only.df$lowFreq <- "lowFreq"
idx <- which(tumor.only.df$GnomAD>1e-4)
tumor.only.df[idx,]$lowFreq <- "HighFreq"
tumor.df2 <- as.data.frame(table(tumor.only.df$lowFreq))

tumor.only.df$lowFreq <- "lowFreq"
idx <- which(tumor.only.df$KOREAN>1e-4)
tumor.only.df[idx,]$lowFreq <- "HighFreq"
tumor.df3 <- as.data.frame(table(tumor.only.df$lowFreq))

tumor.only.df$lowFreq <- "lowFreq"
idx <- which(tumor.only.df$TOMMO>1e-4)
tumor.only.df[idx,]$lowFreq <- "HighFreq"
tumor.df4 <- as.data.frame(table(tumor.only.df$lowFreq))

norm.df1$ratio <- paste0(as.integer(norm.df1$Freq/sum(norm.df1$Freq)*100), "%")
norm.df2$ratio <- paste0(as.integer(norm.df2$Freq/sum(norm.df2$Freq)*100), "%")
norm.df3$ratio <- paste0(as.integer(norm.df3$Freq/sum(norm.df3$Freq)*100), "%")
norm.df4$ratio <- paste0(as.integer(norm.df4$Freq/sum(norm.df4$Freq)*100), "%")

tumor.df1$ratio <- paste0(as.integer(tumor.df1$Freq/sum(tumor.df1$Freq)*100), "%")
tumor.df2$ratio <- paste0(as.integer(tumor.df2$Freq/sum(tumor.df2$Freq)*100), "%")
tumor.df3$ratio <- paste0(as.integer(tumor.df3$Freq/sum(tumor.df3$Freq)*100), "%")
tumor.df4$ratio <- paste0(as.integer(tumor.df4$Freq/sum(tumor.df4$Freq)*100), "%")


p1 <- ggpie(data=norm.df1, x="Freq", label="ratio", fill = "Var1", lab.pos = "in") + ggtitle("Normal, dbGaP")
p2 <- ggpie(data=norm.df2, x="Freq", label="ratio", fill = "Var1", lab.pos = "in") + ggtitle("GnomAD")
p3 <- ggpie(data=norm.df3, x="Freq", label="ratio", fill = "Var1", lab.pos = "in") + ggtitle("KOREAN")
p4 <- ggpie(data=norm.df4, x="Freq", label="ratio", fill = "Var1", lab.pos = "in")+ ggtitle("TOMMO")

p5 <- ggpie(data=tumor.df1, x="Freq", label="ratio", fill = "Var1", lab.pos = "in")
p6 <- ggpie(data=tumor.df2, x="Freq", label="ratio", fill = "Var1", lab.pos = "in")
p7 <- ggpie(data=tumor.df3, x="Freq", label="ratio", fill = "Var1", lab.pos = "in")
p8 <- ggpie(data=tumor.df4, x="Freq", label="ratio", fill = "Var1", lab.pos = "in")

p <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4)
p
ggsave(p, file = "mut_freq.pdf")



#######################################################################################################################################
##  Figure 4I
#######################################################################################################################################
p <- FeaturePlot(obj, features = c("16:69745145G>A/-", "MT:8362T>G/+"), pt.size = 2, order = TRUE, coord.fixe= TRUE) & scale_colour_gradientn(colours = c("grey",brewer.pal(n = 8, name = "Reds")))





saveRDS(obj, file = "p8.rds")
               

