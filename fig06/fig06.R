require(Yano)
require(data.table)

## Directory scRNA-AHCA cloned from https://github.com/bei-lab/scRNA-AHCA
meta <- fread("./scRNA-AHCA/Cell_barcode_and_corresponding_cell_types_of_AHCA/Annotation_AHCA_alltissues_meta.data_84363_cell.txt")

cells <- meta$V1
major.cluster <- meta$Cell_type_in_merged_data
names(major.cluster) <- cells

exp1 <-ReadPISA("./Bladder/Bladder/outs/exp", "Bladder_cDNA_", cells = cells)
exp2 <- ReadPISA("./Blood/Blood/outs/exp", "Blood_cDNA_", cells = cells)
exp3 <- ReadPISA("./Common_bile_duct/Common_bile_duct/outs/exp", "Common.bile.duct_cDNA_", cells = cells)
exp4 <- ReadPISA("./Esophagus/Esophagus/outs/exp", "Esophagus_cDNA_", cells = cells)
exp5 <- ReadPISA("./Heart/Heart/outs/exp", "Heart_cDNA_", cells = cells)
exp6 <- ReadPISA("./Marrow/Marrow/outs/exp", "Marrow_cDNA_", cells = cells)
exp7 <- ReadPISA("./Skin/Skin/outs/exp", "Skin_cDNA_", cells = cells)
exp8 <- ReadPISA("./Stomach/Stomach/outs/exp", "Stomach_cDNA_", cells = cells)
exp9 <- ReadPISA("./Liver/Liver/outs/exp", "Liver_cDNA_", cells = cells)
exp10 <- ReadPISA("./Muscle/Muscle/outs/exp", "Muscle_cDNA_", cells = cells)
exp11 <- ReadPISA("./Small_intestine/Small_intestine/outs/exp", "Small.intestine_cDNA_", cells = cells)
exp12 <- ReadPISA("./Trachea/Trachea/outs/exp", "Trachea_cDNA_", cells = cells)
exp13 <- ReadPISA("./Lymph_node/Lymph_node/outs/exp", "Lymph.node_cDNA_", cells = cells)
exp14 <- ReadPISA("./Rectum/Rectum/outs/exp", "Rectum_cDNA_", cells = cells)
exp15 <- ReadPISA("./Spleen/Spleen/outs/exp", "Spleen_cDNA_", cells = cells)

exp <- mergeMatrix(exp1, exp2, exp3, exp4, exp5, exp6, exp7, exp8, exp9, exp10, exp11, exp12, exp13, exp14, exp15)

obj <- CreateSeuratObject(exp, min.cells=20)
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
obj$tissue <- gsub("(.*)_cDNA.*","\\1",colnames(obj))
obj$major.cluster <- major.cluster[colnames(obj)]

df <- data.frame("tsne_1" = meta$tSNE_1,"tsne_2" = meta$tSNE_2)
df <- as.matrix(df)
rownames(df) <- meta$V1
names(obj)
obj[['tsne']] <- CreateDimReducObject(embeddings = df, assay = "RNA")
Idents(obj) <- obj$major.cluster

cells <- colnames(obj)
var1 <- ReadPISA("./Bladder/Bladder/outs/var", "Bladder_cDNA_", cells = cells)
var2 <- ReadPISA("./Blood/Blood/outs/var", "Blood_cDNA_", cells = cells)
var3 <- ReadPISA("./Common_bile_duct/Common_bile_duct/outs/var", "Common.bile.duct_cDNA_", cells = cells)
var4 <- ReadPISA("./Esophagus/Esophagus/outs/var", "Esophagus_cDNA_", cells = cells)
var5 <- ReadPISA("./Heart/Heart/outs/var", "Heart_cDNA_", cells = cells)
var6 <- ReadPISA("./Marrow/Marrow/outs/var", "Marrow_cDNA_", cells = cells)
var7 <- ReadPISA("./Skin/Skin/outs/var", "Skin_cDNA_", cells = cells)
var8 <- ReadPISA("./Stomach/Stomach/outs/var", "Stomach_cDNA_", cells = cells)
var9 <- ReadPISA("./Liver/Liver/outs/var", "Liver_cDNA_", cells = cells)
var10 <- ReadPISA("./Muscle/Muscle/outs/var", "Muscle_cDNA_", cells = cells)
var11 <- ReadPISA("./Small_intestine/Small_intestine/outs/var", "Small.intestine_cDNA_", cells = cells)
var12 <- ReadPISA("./Trachea/Trachea/outs/var", "Trachea_cDNA_", cells = cells)
var13 <- ReadPISA("./Lymph_node/Lymph_node/outs/var", "Lymph.node_cDNA_", cells = cells)
var14 <- ReadPISA("./Rectum/Rectum/outs/var", "Rectum_cDNA_", cells = cells)
var15 <- ReadPISA("./Spleen/Spleen/outs/var", "Spleen_cDNA_", cells = cells)

var <- mergeMatrix(var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, var13, var14, var15)
obj[["var"]] <- CreateAssayObject(var[,cells], min.cells = 20)

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

DefaultAssay(obj) <- "var"
obj <- NormalizeData(obj)
obj <- ParseVarName(obj)
obj <- RunAutoCorr(obj)
obj <- SetAutoCorrFeatures(obj)
obj <- RunBlockCorr(obj, bind.name = "locus", mode=2)
obj <- RunBlockCorr(obj, bind.name = "locus", mode=1, prefix="locus1") # for Figure 1D

################################################################################################################################
##        Figure 6A
################################################################################################################################
obj$cluster <- gsub("(.*Cell).*", "\\1",obj$major.cluster)
obj$cluster <- gsub("(Fibroblast).*", "\\1",obj$cluster)
obj$cluster <- gsub("(Macrophage).*", "\\1",obj$cluster)
obj$cluster <- gsub("(Keratinocyte).*", "\\1",obj$cluster)

clusters <- unique(obj$cluster)
ids <- c(1:length(clusters))
names(ids) <- clusters
cluster_1 <- paste0(ids, " ", clusters)

ids0 <- ids[obj$cluster]
names(ids0) <- colnames(obj)
obj[['ids']] <- ids0
obj$cluster_1 <- factor(paste0(obj$ids, " ", obj$cluster), levels = cluster_1)
tissues <- unique(obj$tissue)
tissue.cols <- c("blue","#D9D9D9","green","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#666666", "#FF7F00", "#CAB2D6", "#6A3D9A", "black", "#B15928")
names(tissue.cols) <- tissues
cell.types <- unique(obj$cluster)
type.cols1 <- c(rainbow(5)[2:5], brewer.pal(12, name = "Paired"),  brewer.pal(8, name = "Dark2"))
type.cols2 <- type.cols1
type.cols3 <- type.cols1
names(type.cols1) <- unique(obj$ids)
names(type.cols2) <- unique(obj$cluster_1)
names(type.cols3) <- cell.types

p1 <- DimPlot(obj, reduction = "tsne", pt.size=0.5, group.by = "ids", label = TRUE, label.size = 10, cols = type.cols1) + theme_void() + theme(legend.position="none",plot.title = element_blank())
p2 <- DimPlot(obj, reduction = "tsne", group.by = "cluster_1", cols = type.cols2)
leg <- cowplot::get_legend(p2)
ggsave(p1, file = "major_cluster.png", width = 10, height=10)
ggsave(leg, file = "major_cluster_legend.pdf", width = 10, height=9)

################################################################################################################################
##        Figure 6B
################################################################################################################################
p1 <- DimPlot(obj, reduction = "tsne", group.by = "tissue", pt.size=0.5, cols = tissue.cols) + theme_void()+theme(legend.position="none", plot.title = element_blank())
p2 <- DimPlot(obj, reduction = "tsne", group.by = "tissue", cols=tissue.cols)
leg <- cowplot::get_legend(p2)
ggsave(p1, file = "tissue.png", width = 10, height=10)
ggsave(leg, file = "tissue_legend.pdf", width = 10, height=9)

################################################################################################################################
##        Figure 6C
################################################################################################################################

gtf <- gtf2db("/home/projects/ku_00009/data/TranscriptScar/GRCh38/hg38.110.ncbiRefSeq.gtf.gz")
obj <- annoVAR(obj, gtf = gtf, fasta = "/home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-GRCh38-2020-A/fasta/genome.fa")

## GCF_000001405.40_Freq.vcf.gz is manually edited and prased from https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz
obj <- annoVAR(obj,vcf = "/home/projects/ku_00009/data/TranscriptScar/GRCh38/var/GCF_000001405.40_Freq.vcf.gz", tags = c("GnomAD", "dbGaP_PopFreq", "TOMMO", "KOREAN", "RS"))

obj[['var']][[]] -> df
obj[['var']][['chr0']] <- gsub("chr","", df$chr)

obj <- annoVAR(obj, vcf = "/home/projects/ku_00009/data/TranscriptScar/GRCh38/var/clinvar_20240317.vcf.gz", tags = c("CLNSIG"), check.alt.only = TRUE, chr = "chr0")

obj[['var']][[]] -> df

renames <- c("down/up_1K",
             "antisense", "antisense",
             "down/up_1K",
             "down/up_1K",
             "exon_loss",
             "frameshift", "frameshift",
             "intergenic",
             "intron",
             "missense", "noncoding", "noncoding", "noncoding",
             "ref", "splice_region","splice_region","splice_region","splice_sites",
             "start_loss", "synonymous",
             "stop_gained", "stop_loss", "stop_retained", "synonymous", "down/up_1K", "utr5/3","utr5/3")
names(renames) <-  c("antisense_downstream_1K", "antisense_exon", "antisense_intron", "antisense_upstream_1K", "downstream_1K", "exon_loss", "frameshift_elongation", "frameshift_truncation", "intergenic", "intron", "missense", "noncoding_exon", "noncoding_intron", "noncoding_splice_region", "ref", "splice_acceptor", "splice_donor", "splice_region", "splice_sites", "start_loss", "start_retained", "stop_gained", "stop_loss", "stop_retained", "synonymous", "upstream_1K", "utr3", "utr5")

obj[['var']][["conseq"]] <- renames[df$consequence]
df$freq <- "<0.01%"
idx <- which(df$dbGaP_PopFreq > 0.0001|df$consequence %in% "ref")
df[idx, "freq"] <- ">=0.01%"
obj[['var']][["freq"]] <- df$freq

demo.sel <- c( "chr11:320649G=/-",    "chr11:320649G>A/-",  "chrX:119468424C=/+",  "chrX:119468424C>G/+") 
p <- FbtPlot(obj, val="locus.padj", sel.chrs=c(1:22,"X","M"), remove.chr = TRUE, col.by='conseq', pt.size=2, cols=c("#1F78B4","blue","#A6CEE3","#33A02C","#FB9A99","#E31A1C","#FDBF6F","grey97","#CAB2D6","#6A3D9A","#FFFF99", "hotpink3", "seashell2", "#131313", "#B2DF8A", "palegreen"), shape.by = "freq", point.label = demo.sel, label.size = 5)
p1 <- p+theme(legend.position="none")+xlab("")+ylab("")
ggsave(p1, file = "locus_pval_1.png", width=16,height = 6)
ggsave(cowplot::get_legend(p), file = "locus_legend.pdf", width = 10,height = 10)

################################################################################################################################
##        Figure 6D, 6E
################################################################################################################################
demo.sel <- c( "chr11:320649G=/-",    "chr11:320649G>A/-",  "chr9:33447422C=/-",  "chr9:33447422C>T/-")
p <- FeaturePlot(obj, features = demo.sel, reduction = "tsne", ncol=4, pt=1, order = TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p, file = "demo_features.pdf", width = 23, height = 6)

################################################################################################################################
##        Figure 6F
################################################################################################################################

df <- obj[['var']][[]]
df  %>% filter(locus.pval < 1e-10) %>% pull(locus) -> sel0
df  %>% filter(locus %in% sel0) %>% rownames %>% sort ->sel

tc0 <- paste0(obj$tissue,"_",obj$cluster)
tc <- table(tc0)
sel.tc <- names(tc[tc>50])

idx <- which(tc0 %in% sel.tc)

obj1 <- AggregateExpression(object = obj[,idx], group.by = c("tissue", "cluster"), return.seurat = TRUE, features = sel)

dat <- GetAssayData(obj1, layer = "counts")
dat@x[which(dat@x<5)] <- 0
dat@x[which(dat@x>=5)] <- 1
dat <- as.matrix(dat)
rs <- rowSums(dat)
idx <- which(rs > 1 & rs < 103)
dat <- dat[idx,]
dat <- t(dat)
d0 <- dist(dat, method = "binary")
hc0 <- hclust(d0, method = "ward")
dend0 <- as.dendrogram(hc0)
km0 <- cutree(hc0,10)
km0 <- as.factor(km0)

d <- dist(t(dat), method = "binary")
hc <- hclust(d, method = "ward")
dend <- as.dendrogram(hc)
km <- cutree(hc,10)
km <- as.factor(km)

df[["cluster"]] <- "others"
df[names(km),]$cluster <- km
idx <- which(df$gene_name<1e-10 &df$cluster == "others")
df[idx,]$cluster <- "filtered"

obj[['var']][["cluster"]] <- factor(df[['cluster']], levels = gtools::mixedsort(unique(df$cluster)))

eats <- colnames(dat)

ha <- rowAnnotation(tissue=gsub("(.*)_(.*)", "\\1",rownames(dat)), type=gsub("(.*)_(.*)", "\\2",rownames(dat)), col = list(tissue = tissue.cols, type = type.cols3), border = TRUE)

ann <- df[colnames(dat),]$anno
names(ann) <- colnames(dat)

ann.cols <- c("white", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "black", "red")
names(ann.cols) <- factor(c("Others", "IFITM",  "IGH",   "IGL",   "IGK",   "MHC", "mito"), levels = c("Others", "IFITM",  "IGH",   "IGL",   "IGK",   "MHC", "mito"))

df0 <- df[names(km),]

df0$nm <- rownames(df0)
df0 %>%  arrange(cluster, anno) %>% pull(nm) -> sel.ord

ha1 <- HeatmapAnnotation(anno=ann[sel.ord],cluster=km[sel.ord],col = list(anno=ann.cols, cluster=cluster.cols), border = TRUE)

pdf(file = "cluster_ht.pdf", width = 26, height = 18, family="Arial")
Heatmap(dat[,sel.ord], border =TRUE, use_raster=TRUE, col = c("white", "black"), show_column_names=FALSE, right_annotation = ha, top_annotation = ha1, cluster_rows = dend0, row_split=10, cluster_columns = FALSE)
dev.off()

saveRDS(obj, file= "AHCA.rds")
