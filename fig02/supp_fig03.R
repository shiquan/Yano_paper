require(Yano)
require(RColorBrewer)
require(ComplexHeatmap)

obj.sm <- readRDS("object_sm.rds")

obj.sm[['exon']][[]] %>% filter(gene_name.padj<1e-30) %>% rownames -> exons
obj.sm[['exon']][[]] %>% filter(gene_name.padj<1e-30) %>% pull(gene_name) -> bind.genes
DefaultAssay(obj.sm) <- "RNA"
obj.sm <- ScaleData(obj.sm, features = unique(bind.genes))
DefaultAssay(obj.sm) <- "exon"
obj.sm <- ScaleData(obj.sm, features = exons)

dat1 <- GetAssayData(obj.sm, assay = 'exon', layer = 'scale.data')
dat2 <- GetAssayData(obj.sm, assay = 'RNA', layer = 'scale.data')
idents <- sort(Idents(obj.sm))
order.cells <- names(idents)

dat2 <- dat2[bind.genes,]
rownames(dat2) <- exons

dat <- cbind(dat1, dat2)


d <- dist(dat)
hc <- hclust(d)
idx <- hc$labels[hc$order]

ha <- HeatmapAnnotation(group=idents, border = TRUE)
ht1 <- Heatmap(dat1[idx, order.cells], cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, border = TRUE, use_raster = TRUE, top_annotation = ha, name = "exon", column_title = "exon")
ht2 <- Heatmap(dat2[idx, order.cells], cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, border = TRUE, use_raster = TRUE, top_annotation = ha, name = "gene", column_title = "gene", row_names_max_width = max_text_width(rownames(dat2), gp = gpar(fontsize = 12)))

ht <- ht1 + ht2
pdf(file = "ht_sm2.pdf", width=15, height = 22)
draw(ht, heatmap_legend_side = "left",  annotation_legend_side = "left")
dev.off()
