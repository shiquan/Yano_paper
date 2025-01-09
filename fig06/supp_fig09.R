require(Yano)
require(RColorBrewer)

obj <- readRDS("AHCA.rds")

df <- obj[['var']][[]]
p.cutoff <1e-10
df %>% filter(locus.padj<p.cutoff) %>% pull(locus) %>% table -> tb

names(which(tb > 1)) -> sel
names(which(tb > 2)) %>% length # 7

df %>% filter(locus %in% sel) %>% pull(gene_name) %>% unique -> genes

ifs <- grep("^IFITM", genes, value = TRUE)
df %>% filter(locus %in% sel & gene_name %in% ifs) %>% nrow
genes <- grep("^IG", genes, value = TRUE, invert = TRUE)
genes <- grep("^HLA", genes, value = TRUE, invert = TRUE)
genes <- grep("^TRB", genes, value = TRUE, invert = TRUE)
genes <- grep("^IFITM", genes, value = TRUE, invert = TRUE)
genes <- grep("^LOC", genes, value = TRUE, invert = TRUE)
genes <- grep("^PRAMENP", genes, value = TRUE, invert = TRUE)
genes <- grep("^AQP3", genes, value = TRUE, invert = TRUE)

df %>% filter(locus %in% sel & gene_name %in% genes) %>% rownames -> sel0
df[sel0,]%>% pull(locus) -> sel

df %>% filter(locus %in% sel) %>% rownames -> sel0
sel0 <- gtools::mixedsort(sel0)

########################################################################################################################################
##                  Supplementary Figure 9
########################################################################################################################################

p <- FeaturePlot(obj, features = sel0[c(1:42)],  reduction = "tsne", ncol=6, pt=1, order = TRUE, coord.fixed=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & theme_void() & theme(legend.position="none")
ggsave(p, file = "ASE_features_1.png", width =15, height=23)

########################################################################################################################################
##                  Supplementary Figure 10A
########################################################################################################################################

p <- FeaturePlot(obj, features = sel0[43:60],  reduction = "tsne", ncol=6, pt=1, order = TRUE, coord.fixed=TRUE)& scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & theme_void() & theme(legend.position="none")
#ggsave(p, file = "ASE_features.pdf", width =15, height=20)
ggsave(p, file = "ASE_features_2.png", width =15, height=8)
