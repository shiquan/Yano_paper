require(Yano)

obj <- readRDS("p8.rds")

#############################################################################################################################
## Figure 5A
#############################################################################################################################
p8.k <- subset(obj, level2_celltype %in% "Keratinocyte")
DefaultAssay(p8.k) <- "var"
p8.k <- QuickRecipe(p8.k)
p <- DimPlot(p8.k, group.by = "sample", pt.size=2)
#p2 <- DimPlot(p8.k, group.by = "level3_celltype", pt.size=2)
#p <- cowplot::plot_grid(p1,p2)

ggsave(p, file="p8_subset.pdf")

p8.k <- NormalizeData(p8.k)
p8.k <- RunAutoCorr(p8.k)
p8.k <- SetAutoCorrFeatures(p8.k)
p8.k <- ParseVarName(p8.k)
p8.k <- RunBlockCorr(p8.k, bind.name = "locus", mode=2)

#############################################################################################################################
## Figure 5B
#############################################################################################################################

demo.sel <- c("12:52880996C=/-",  "12:52880996C>G/-")

p <- FbtPlot(p8.k, val="locus.padj", sel.chrs = sel.chrs, point.label = demo.sel, label.size = 6, xlab = "", ylab="", size=2, col.by = "dbGaP", shape="var.type") + theme(legend.position="none")
ggsave(p, file="p8_k_var.png",width = 16, height = 6) 

#############################################################################################################################
## Figure 5C
#############################################################################################################################

p <- FeaturePlot(p8.k, features = c("12:52880996C=/-",  "12:52880996C>G/-"),pt.size=3, order = TRUE) & scale_colour_gradientn(colours = c("grey",brewer.pal(n = 8, name = "Reds")))
ggsave(p, file = "fig5_features.pdf")

#############################################################################################################################
## Figure 5D
#############################################################################################################################

p8.k[['var']][[]] %>% filter (locus.padj<1e-10) %>% rownames -> sel
dat <- FetchData(p8.k, vars = sel)
dat <- as.matrix(dat)
d <- dist(dat)
hc <- hclust(d)
cells <- rownames(dat)
an <- gsub("P8_(.*)_.*", "\\1",rownames(dat))
names(an) <- rownames(dat)

hc <- hclust(dist(dat), method = 'ward.D2')
n<-5
mode.cols <- c( "blue", "#FF7F00", "#FDBF6F", "#131313", "#1F78B4" )
names(mode.cols) <- c("1","2","3","4","5")

km <- cutree(hc, k=n)
km <- as.factor(km)
df <- p8.k[['var']][[]]

row_dend = as.dendrogram(hc)

pdf("ht.pdf", width = 10, height = 6)
Heatmap(km, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, border=TRUE, col= mode.cols, cluster_rows = row_dend, width = unit(5, "mm"), row_split = 5)  +Heatmap(dat, use_raster = TRUE, col= c(brewer.pal(n = 9, name = "Reds")), border=TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE) + Heatmap(an, col = c("black","red"), show_column_names=FALSE,show_row_names = FALSE, use_raster = TRUE, border=TRUE, width = unit(5, "mm"))
dev.off()

#############################################################################################################################
## Figure 5E, 5F
#############################################################################################################################
p8.k$subclones <- km[colnames(p8.k)]
Idents(p8.k) <- p8.k$subclones
p8.k[['var']][[]] %>% filter(locus.padj<1e-10) %>% rownames -> hvf
p8.k1 <- ScaleData(p8.k, features = hvf)
p8.k1 <- RunPCA(p8.k1, features = hvf, npcs = 5)
p8.k1 <- FindNeighbors(p8.k1, dims = 1:5)
p8.k1 <- FindClusters(p8.k1)
p8.k1 <- RunUMAP(p8.k1, dims=1:5)
p8.k1$subclones <- km[colnames(p8.k1)]
Idents(p8.k1) <- p8.k1$subclones

p0 <- DimPlot(p8.k, pt.size = 2, group.by = "subclones", shape.by = "sample",cols=mode.cols) + ggtitle("") + coord_fixed()
p1 <- DimPlot(p8.k1, pt.size=2, shape.by="sample", cols = mode.cols) + coord_fixed()
p0 <- p0 + theme(legend.position="none")
p1 <- p1 + theme(legend.position="none")
p <- cowplot::plot_grid(p0,p1)
ggsave(p, file = "recluster_subclone.pdf", width = 11, height = 5)


#############################################################################################################################
## Figure 5G
#############################################################################################################################
DefaultAssay(p8.k1) <- "RNA"

markers1 <- FindMarkers(p8.k1, ident.1 = "1", ident.2 = c("2","3"))
markers1 %>% dplyr::filter(avg_log2FC > 1 & p_val_adj<1e-5) %>% head(10) %>% rownames -> sel1

markers2 <- FindMarkers(p8.k1, ident.1 = "2", ident.2 = c("1","3"))
markers2 %>% dplyr::filter(avg_log2FC > 1 & p_val_adj<1e-5) %>% head(10) %>% rownames -> sel2

markers3 <- FindMarkers(p8.k1, ident.1 = "3", ident.2 = c("1","2"))
markers3 %>% dplyr::filter(avg_log2FC > 1 & p_val_adj<1e-5) %>% head(10) %>% rownames -> sel3

markers4 <- FindMarkers(p8.k1, ident.1 = "4", ident.2 = c("1", "2", "3"))
markers4 %>% dplyr::filter(avg_log2FC > 1 & p_val_adj<1e-10) %>% head(10) %>% rownames -> sel4

markers5['KRT6A',]

markers5 <- FindMarkers(p8.k1, ident.1 = "5")
markers5 %>% dplyr::filter(avg_log2FC > 1 & p_val_adj<1e-10) %>% head(10) %>% rownames -> sel5

p8.k1 <- ScaleData(p8.k1, features = unique(c(sel1,sel2, sel3, sel4, sel5)))
p <- DoHeatmap(p8.k1, features = c(sel1, sel2,sel3, sel4, sel5), raster = TRUE, group.colors = mode.cols)
ggsave(p, file = "subclones_ht.pdf", width = 10, height = 6)

#############################################################################################################################
## Figure 5H
#############################################################################################################################

p <- VlnPlot(obj,features = c("CCL3","CCR1"), group.by = "sample" )
ggsave(p, file="CCL3_CCR1_vln.pdf")


