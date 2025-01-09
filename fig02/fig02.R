require(Yano)
require(data.table)
require(harmony)
require(scCustomize)
require(RColorBrewer)

########################################################################
##                         Process 10X data                           ##
########################################################################

# The annotation files download from https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/
## And rename the annotation files with capped "10X_" to distinguish with Smartseq data.
meta <- as.data.frame(fread("10X_sample_metadata.csv"))
mem <- as.data.frame(fread("10X_cluster.membership.csv"))
anno <- as.data.frame(fread("10X_cluster.annotation.csv"))
rownames(anno) <- anno$cluster_id

meta <- left_join(meta, mem, by = "V1")
rownames(meta) <- meta$V1
meta <- meta[mem$V1,]

cells <- gsub("-[0-9]+","_",meta$V1)
rownames(meta) <- cells

Read1 <- function(dir = NULL, suffix = NULL, cells = NULL) {
  exp <- ReadPISA(dir, suffix = suffix)
  colnames(exp) <- gsub("-[0-9]+","_",colnames(exp))
  cells1 <- intersect(colnames(exp), cells)
  exp[,cells1]
}

# We found four libraries have no annotation information, so we filter out these libraries
exp1 <- Read1("10xV3/L8TX_181211_01_A01_CellRangerV5/exp/", suffix = "L8TX_181211_01_A01", cells=cells)
exp2 <- Read1("10xV3/L8TX_190430_01_B08_CellRangerV5/exp/", suffix = "L8TX_190430_01_B08", cells = cells)
#exp3 <- Read1("10xV3/L8TX_190430_01_E08_CellRangerV5/exp/", suffix = "L8TX_190430_01_E08", cells =cells)
#exp4 <- Read1("10xV3/L8TX_190430_01_H08_CellRangerV5/exp/", suffix = "L8TX_190430_01_H08", cells =cells)
exp5 <- Read1("10xV3/L8TX_181211_01_B01_CellRangerV5/exp/", suffix = "L8TX_181211_01_B01", cells =cells)
#exp6 <- Read1("10xV3/L8TX_190430_01_C08_CellRangerV5/exp/", suffix = "L8TX_190430_01_C08", cells =cells)
exp7 <- Read1("10xV3/L8TX_190430_01_F08_CellRangerV5/exp/", suffix = "L8TX_190430_01_F08", cells =cells)
exp8 <- Read1("10xV3/L8TX_190430_01_A08_CellRangerV5/exp/", suffix = "L8TX_190430_01_A08", cells =cells)
#exp9 <- Read1("10xV3/L8TX_190430_01_D08_CellRangerV5/exp/", suffix = "L8TX_190430_01_D08", cells =cells)
exp10 <- Read1("10xV3/L8TX_190430_01_G08_CellRangerV5/exp/", suffix = "L8TX_190430_01_G08", cells =cells)

## merge All gene counts into one big matrix
exp <- mergeMatrix(exp1, exp2, exp5, exp7, exp8, exp10)

## Perform standard Seurat pipeline
obj <- CreateAssayObject(exp, min.cells=20)

## Add meta data back to the object
meta <- meta[colnames(obj),]
obj <- AddMetaData(obj, metadata = meta)

anno <- anno[obj$x,]
rownames(anno) <- colnames(obj)
obj <- AddMetaData(obj, anno)

## Filter out low quality
obj <- subset(obj, subclass_label != "Low Quality" & subclass_label != "doublet")

# We only select neurons for analysis
obj <- subset(obj, class_label %in% c("Glutamatergic")) 

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj,verbose=FALSE)
obj <- RunHarmony(obj, "library_id")
obj <- FindNeighbors(obj, dims = 1:30, verbose=FALSE, reduction = "harmony") %>% RunUMAP(dims=1:20, verbose=FALSE, reduction = "harmony")

cells <- colnames(obj)

exon1 <- Read1("10xV3/L8TX_181211_01_A01_CellRangerV5/exon/", suffix = "L8TX_181211_01_A01", cells=cells)
exon2 <- Read1("10xV3/L8TX_190430_01_B08_CellRangerV5/exon/", suffix = "L8TX_190430_01_B08", cells = cells)
exon5 <- Read1("10xV3/L8TX_181211_01_B01_CellRangerV5/exon/", suffix = "L8TX_181211_01_B01", cells =cells)
exon7 <- Read1("10xV3/L8TX_190430_01_F08_CellRangerV5/exon/", suffix = "L8TX_190430_01_F08", cells =cells)
exon8 <- Read1("10xV3/L8TX_190430_01_A08_CellRangerV5/exon/", suffix = "L8TX_190430_01_A08", cells =cells)
exon10 <- Read1("10xV3/L8TX_190430_01_G08_CellRangerV5/exon/", suffix = "L8TX_190430_01_G08", cells =cells)

exon <- mergeMatrix(exon1, exon2, exon5, exon7, exon8, exon10)

junction1 <- Read1("10xV3/L8TX_181211_01_A01_CellRangerV5/junction/", suffix = "L8TX_181211_01_A01", cells=cells)
junction2 <- Read1("10xV3/L8TX_190430_01_B08_CellRangerV5/junction/", suffix = "L8TX_190430_01_B08", cells = cells)
junction5 <- Read1("10xV3/L8TX_181211_01_B01_CellRangerV5/junction/", suffix = "L8TX_181211_01_B01", cells =cells)
junction7 <- Read1("10xV3/L8TX_190430_01_F08_CellRangerV5/junction/", suffix = "L8TX_190430_01_F08", cells =cells)
junction8 <- Read1("10xV3/L8TX_190430_01_A08_CellRangerV5/junction/", suffix = "L8TX_190430_01_A08", cells =cells)
junction10 <- Read1("10xV3/L8TX_190430_01_G08_CellRangerV5/junction/", suffix = "L8TX_190430_01_G08", cells =cells)
junction <- mergeMatrix(junction1, junction2, junction5, junction7, junction8, junction10)

exclude1 <- Read1("10xV3/L8TX_181211_01_A01_CellRangerV5/exclude/", suffix = "L8TX_181211_01_A01", cells=cells)
exclude2 <- Read1("10xV3/L8TX_190430_01_B08_CellRangerV5/exclude/", suffix = "L8TX_190430_01_B08", cells = cells)
exclude5 <- Read1("10xV3/L8TX_181211_01_B01_CellRangerV5/exclude/", suffix = "L8TX_181211_01_B01", cells =cells)
exclude7 <- Read1("10xV3/L8TX_190430_01_F08_CellRangerV5/exclude/", suffix = "L8TX_190430_01_F08", cells =cells)
exclude8 <- Read1("10xV3/L8TX_190430_01_A08_CellRangerV5/exclude/", suffix = "L8TX_190430_01_A08", cells =cells)
exclude10 <- Read1("10xV3/L8TX_190430_01_G08_CellRangerV5/exclude/", suffix = "L8TX_190430_01_G08", cells =cells)
exclude <- mergeMatrix(exclude1, exclude2, exclude5, exclude7, exclude8, exclude10)

obj[['junction']] <- CreateAssayObject(junction[,colnames(obj)], min.cells=20)
obj[['exon']] <- CreateAssayObject(exon[,colnames(obj)], min.cells=20)
obj[['exclude']] <- CreateAssayObject(exclude[,cells],min.cells=20)

DefaultAssay(obj.10x) <- "RNA"

##################################################################################
##            Figure 2A, left
##################################################################################
p <- DimPlot(obj, pt.size=1, group.by = "subclass_label") + ggtitle("")
ggsave(p, file = "fig02_A.pdf", width = 7, height = 6)

DefaultAssay(obj) <- 'exon'
obj <- NormalizeData(obj)
obj <- RunAutoCorr(obj, reduction = "harmony")
obj <- SetAutoCorrFeatures(obj)
obj <- ParseExonName(obj)
obj <- RunBlockCorr(obj, bind.assay = "RNA", wm.name = "harmony_wm", bind.name = "gene_name") # Mode 1

DefaultAssay(obj) <- 'junction'
obj <- NormalizeData(obj)
obj <- RunAutoCorr(obj, reduction = "harmony")
obj <- SetAutoCorrFeatures(obj)
obj <- ParseExonName(obj)
obj <- RunBlockCorr(obj, bind.assay = "RNA", wm.name = "harmony_wm", bind.name = "gene_name") # Mode 1

#####################################################################################
##                   Figure 2B
#####################################################################################
p <- FbtPlot(obj, assay = c("exon", "junction"), shape.by="assay", col.by = "assay", val = "gene_name.padj", sel.chrs=c(1:21,"X"), remove.chr = TRUE, point.label = c("chr10:68136486-68136626/-/Arid5b", "chr7:140082462-140082548/-/Caly", "chr5:142904067-142904248/-/Actb"), cols = c("red", "blue"), pt.size=2, label.size = 6) + xlab("") + ylab("") + theme(legend.position="none")
ggsave(p, file = "fbtplot_10x.png", width = 15, height = 6)

cells <- colnames(obj)
cells.1 <- sample(cells, length(cells)/2)
cells.2 <- setdiff(cells, cells.1)
DefaultAssay(obj) <- 'exon'
obj <- RunBlockCorr(obj, bind.assay = "RNA", prefix="group2_rep1", cells = cells.1, reduction = "harmony")
obj <- RunBlockCorr(obj, bind.assay = "RNA", prefix="group2_rep2", cells = cells.2, reduction = "harmony")

#######################################################################################
##                   Figure 2D
#######################################################################################
p <- ggplot(data=df1) + geom_point(aes(-log10(group2_rep1.padj), -log10(group2_rep2.padj))) + theme_pubr() + xlab("") + ylab("")
ggsave(p, file = "10X_group2_replicates.png", width=4, height = 4)
subset(df1, !is.na(group2_rep1.padj) & !is.na(group2_rep2.padj)) -> df0
cor(-log10(df0$group2_rep1.padj), -log10(df0$group2_rep2.padj)) ## should be 0.97

saveRDS(obj, file = "object_10X.rds")
rm(obj)
gc()

#####################################################################################
##                      Smart-seq2                                                 ##
#####################################################################################

exp <- ReadPISA("Smartseq2/exp")
exon <- ReadPISA("Smartseq2/exon/")
junction <- ReadPISA("Smartseq2/junction/")
exclude <- ReadPISA("Smartseq2/exclude/")

obj.sm <- QuickRecipe(exp, min.cells = 20)
obj.sm[['exon']] <- CreateAssayObject(exon, min.cells=20)
obj.sm[['junction']] <- CreateAssayObject(junction, min.cells=20)
obj.sm[['exclude']] <- CreateAssayObject(exclude, min.cells=20)

## The meta data for Smartseq2 download form https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/processed/analysis/SMARTer_cells_MOp/
## And rename the annotation files with capped "Sm_" to distinguish with 10X data.

meta <- as.data.frame(fread("Sm_sample_metadata.csv"))
rownames(meta) <- meta$V1
meta <- meta[colnames(obj.sm),]
obj.sm <- AddMetaData(obj.sm, metadata = meta)

anno <- as.data.frame(fread("Sm_cluster.membership.csv"))
rownames(anno) <- anno$V1
obj.sm$cluster_id <- anno[colnames(obj.sm),"x"]

anno <- as.data.frame(fread("Sm_cluster.annotation.csv"))
rownames(anno) <- anno$cluster_id

anno <- anno[obj.sm$cluster_id,]
rownames(anno) <- colnames(obj.sm)
obj.sm <- AddMetaData(obj.sm, metadata=anno)

obj.sm <- subset(obj.sm, class_label %in% c("Glutamatergic"))

DefaultAssay(obj.sm) <- "RNA"
######################################################################################
##             Figure 2A, right
######################################################################################
p <- DimPlot(obj.sm, pt.size=1, group.by = "subclass_label")+ ggtitle("")
ggsave(p, file = "fig02_A_2.pdf", width = 7, height = 6)

DefaultAssay(obj.sm) <- 'exon'
obj.sm <- NormalizeData(obj.sm)
obj.sm <- RunAutoCorr(obj.sm)
obj.sm <- SetAutoCorrFeatures(obj.sm)
obj.sm <- ParseExonName(obj.sm)
obj.sm <- RunBlockCorr(obj.sm, bind.assay = "RNA", bind.name = "gene_name") # Mode 1

DefaultAssay(obj.sm) <- 'junction'
obj.sm <- NormalizeData(obj.sm)
obj.sm <- RunAutoCorr(obj.sm)
obj.sm <- SetAutoCorrFeatures(obj.sm)
obj.sm <- ParseExonName(obj.sm)
obj.sm <- RunBlockCorr(obj.sm, bind.assay = "RNA", bind.name = "gene_name") # Mode 1

saveRDS(obj.sm, file = "object_sm.rds")

######################################################################################
##             Figure 2C
######################################################################################
p <- FbtPlot(obj.sm, assay = c("exon", "junction"), shape.by="assay", col.by = "assay", val = "gene_name.padj", sel.chrs=c(1:21,"X"), remove.chr = TRUE, point.label = c("chr10:68136486-68136626/Arid5b", "chr7:140082462-140082548/Caly", "chr5:142904067-142904248/Actb"), cols = c("red", "blue"), pt.size=2, label.size = 6) + xlab("") + ylab("") + theme(legend.position="none")
ggsave(p, file = "fbtplot_sm.png", width = 15, height = 6)


####################################################################################
##       Figure 2E,  Downsampling analysis                                        ##
####################################################################################
obj.10x <- readRDS("object_10X.rds")
cells <- coalnames(obj.10x)
cells.500 <- sample(cells, 500)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell500", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 1000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell1000", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 2000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell2000", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 5000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell5000", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 10000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell10K", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 15000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell15K", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 20000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell20K", reduction = "harmony", dims = 1:20)
cells.500 <- sample(cells, 30000)
obj.10x <- RunBlockCorr(obj.10x, bind.assay = "RNA", cells = cells.500, prefix = "cell30K", reduction = "harmony", dims = 1:20)


DefaultAssay(obj.sm) <- "exon"
cells <- colnames(obj.sm)
cells.500 <- sample(cells, 500)
obj.sm <- RunBlockCorr(obj.sm, bind.assay = "RNA", cells = cells.500,  prefix = "cell500")
cells.500 <- sample(cells, 1000)
obj.sm <- RunBlockCorr(obj.sm, bind.assay = "RNA", cells = cells.500,  prefix = "cell1000")
cells.500 <- sample(cells, 2000)
obj.sm <- RunBlockCorr(obj.sm, bind.assay = "RNA", cells = cells.500,  prefix = "cell2000")
cells.500 <- sample(cells, 3000)
obj.sm <- RunBlockCorr(obj.sm, bind.assay = "RNA", cells = cells.500,  prefix = "cell3000")

obj.10x[['exon']][[]] -> df1
obj.sm[['exon']][[]] -> df2

df1 %>% filter(cell500.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.500
df1 %>% filter(cell1000.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.1000
df1 %>% filter(cell2000.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.2000
df1 %>% filter(cell5000.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.5000
df1 %>% filter(cell10K.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.10k
df1 %>% filter(cell15K.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.15k
df1 %>% filter(cell20K.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.20k
df1 %>% filter(cell30K.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.30k
df1 %>% filter(gene_name.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.all

df2 %>% filter(cell500.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.sm.500
df2 %>% filter(cell1000.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.sm.1000
df2 %>% filter(cell2000.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.sm.2000
df2 %>% filter(cell3000.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.sm.3000
df2 %>% filter(gene_name.padj<1e-10) %>% pull(gene_name) %>% unique %>% length -> sel.sm.all

tb <- data.frame(lib = c(rep("10Xv3",9), rep("Smartseq2",5)), cell_num = c(500,1000,2000, 5000, 10000, 15000, 20000, 30000, 33257, 500, 1000, 2000, 3000, 4135), detect = c(sel.500, sel.1000, sel.2000, sel.5000, sel.10k, sel.15k, sel.20k, sel.30k, sel.all, sel.sm.500, sel.sm.1000, sel.sm.2000, sel.sm.3000, sel.sm.all))

p <- ggplot(data=tb) + geom_point(aes(x = cell_num, y = detect, color = lib), size=3) + geom_line(aes(x = cell_num, y = detect, group = lib, color = lib), size=1.5) + theme_pubr()
p <- p + scale_x_log10()
ggsave(p, file = "simulate.pdf", width=3, height=3)

############################################################################################
##                      Figure 2F
###########################################################################################

obj.10x[['exon']][[]] -> df1
obj.sm[['exon']][[]] -> df2

df1$nm1 <- gsub("(.*)/[-+](.*)", "\\1\\2",rownames(df1))
itr <- intersect(df1$nm1, rownames(df3))
df <- subset(df1, nm1 %in% itr)
df$name <- rownames(df)
rownames(df) <- df$nm1
df2$nm1 <- rownames(df2)
df <- df[itr,]
df.1 <- df2[itr,]
df$sm.padj <- df.1$gene_name.pval

subset(df, gene_name.padj <1e-10 & sm.padj<1e-10) %>% nrow
subset(df, gene_name.padj <1e-10 & sm.padj>0.01) %>% nrow
subset(df, gene_name.padj >0.01 & sm.padj<1e-10) %>% nrow

df0 <- df[c("chr10:68136486-68136626/Arid5b", "chr3:80690404-80690518/Gria2", "chr5:142904067-142904248/Actb"),]
p <- ggplot(df) + geom_point(aes(x= -log10(gene_name.padj), y = -log10(sm.padj))) + geom_point(data=df0,aes(-log10(gene_name.padj), -log10(sm.padj)), color = "red")+ theme_pubr() + geom_label_repel(data=df0, aes(x = -log10(gene_name.padj), y = -log10(sm.padj), label = nm1), nudge_y=15, max.iter = 10000000, max.overlaps = 3) + xlab("10xV3") + ylab("Smartseq2")

ggsave(p, "10x_vs_sm2.pdf")

############################################################################################
##                     Figure 2G
############################################################################################
DefaultAssay(obj.10x) <- "exon"
p1 <- FeaturePlot_scCustom(obj.10x, features = c("chr10:68136486-68136626/-/Arid5b"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
DefaultAssay(obj.10x) <- "RNA"
p2 <- FeaturePlot_scCustom(obj.10x ,features = c("Arid5b"), order=TRUE, pt.size=1,figure_plot = TRUE, reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- FeaturePlot_scCustom(obj.sm,  features = c(  "chr10:68136486-68136626/Arid5b"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p4 <-  FeaturePlot_scCustom(obj.sm,  features = c("Arid5b"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p00 <- cowplot::plot_grid(p1, p2, p3, p4, ncol=4)
ggsave(p00, file = "Arid5b.pdf")

###########################################################################################
##     Figure 2H
##########################################################################################
DefaultAssay(obj.10x) <- "exon"
p1 <- FeaturePlot_scCustom(obj.10x, features = c("chr5:142904067-142904248/-/Actb"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
DefaultAssay(obj.10x) <- "RNA"
p2 <- FeaturePlot_scCustom(obj.10x ,features = c("Actb"), order=TRUE, pt.size=1,figure_plot = TRUE, reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- FeaturePlot_scCustom(obj.sm,  features = c("chr5:142904067-142904248/Actb"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p4 <-  FeaturePlot_scCustom(obj.sm,  features = c("Actb"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p11 <- cowplot::plot_grid(p1, p2, p3, p4, ncol=4)
ggsave(p11, file = "Actb.pdf")

########################################################################################
##     Figure 2I
#######################################################################################
p <- RatioPlot(obj.10x, features = c("chr5:142904067-142904248/-/Actb"), assay = "exon", bind.assay = "RNA", bind.name = "gene_name", order = TRUE, pt.size = 1, figure_plot = TRUE, reduction = "umap") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p, file = "Actb_exon_ratio.pdf")

##########################################################################################
##     Figure 2J
#########################################################################################
DefaultAssay(obj.10x) <- "exon"
p1 <- FeaturePlot_scCustom(obj.10x, features = c("chr7:140082462-140082548/-/Caly"), pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
DefaultAssay(obj.10x) <- "RNA"
p2 <- FeaturePlot_scCustom(obj.10x ,features = c("Caly"), order=TRUE, pt.size=1,figure_plot = TRUE, reduction = 'umap') & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- FeaturePlot_scCustom(obj.sm,  features = c("chr7:140082462-140082548/Caly"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p4 <-  FeaturePlot_scCustom(obj.sm,  features = c("Caly"), order=TRUE, pt.size=1, figure_plot = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p22 <- cowplot::plot_grid(p1, p2, p3, p4, ncol=4)
ggsave(p22, file = "Caly.pdf")

