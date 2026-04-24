library(Seurat)
library(plyr)
library(scater)
library(stringr)
library(future)
library(ggplot2)
library(cowplot)
library(reshape2)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)
library(SingleCellExperiment)
library(data.table)
library(tibble)
library(harmony)
library(scales)



mycol <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Set2"),
           brewer.pal(11,"Set3"), brewer.pal(12,"Paired"),
           brewer.pal(8,"Accent"), brewer.pal(11,"Spectral"),
           brewer.pal(11,"BrBG"), brewer.pal(11,"PiYG"),
           brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"))

source("~/additional_fuctions_for_FUSCC.R")

## 数据读取 ======

`1866-counts_DC_cohort1` <- readRDS("~/1866-counts_DC_cohort1.rds")

datafilt <- `1866-counts_DC_cohort1`

# reference: 2021_Bassez_BRCA_ICB_cohort

Bassez <- readRDS("~/ICB_2021_BRCA_Bassez_all.rds")

## data pre-processing ======

datafilt <- CreateSeuratObject(
  counts = datafilt,
  project = "2021_Bassez_BRCA_ICB_DC"
)

# clinical information

common_cell <- intersect(colnames(datafilt),colnames(Bassez)) # 2410

meta <- Bassez@meta.data
meta <- meta[common_cell,]

# merge meta.data

datafilt@meta.data <- meta
datafilt$sample <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", rownames(meta)) # BIOKEY_10_Pre

table(datafilt$patient_id) # 31
table(datafilt$sample) # 61
table(datafilt$response) # non, NR, R

## seurat pipeline ======

nfeatures = 2000
ndim = 15
neigh = 50
dist = 0.5
res = 0.2

datafilt <- NormalizeData(datafilt, scale.factor = 10000,
                          normalization.method = "LogNormalize")

datafilt <- FindVariableFeatures(datafilt, nfeatures = nfeatures, 
                                 selection.method = "vst")

datafilt <- ScaleData(datafilt, features = VariableFeatures(datafilt))

datafilt <- RunPCA(datafilt, assay = 'RNA', slot = 'scale.data')

datafilt <- RunHarmony(datafilt, group.by.vars = "sample",dims.use = 1:50,
                       assay.use = "RNA")

datafilt <- FindNeighbors(datafilt, k.param = neigh,
                          dims = 1:ndim, reduction = "harmony")

datafilt <- FindClusters(datafilt, resolution = 1, n.iter = 50)

datafilt <- RunUMAP(datafilt, dims = 1:ndim,
                    n.neighbors = neigh, min.dist = dist, 
                    reduction = "harmony", reduction.name = "umap_harmony")

datafilt <- RunTSNE(datafilt, dims = 1:ndim,
                    n.neighbors = neigh, min.dist = dist, 
                    reduction = "harmony", reduction.name = "tsne_harmony")


fig1 <- dimplot_new(data = datafilt,
                    reduction = "tsne_harmony",
                    pt.size = 0.1, label = T,
                    group.by = c("seurat_clusters"))
fig1

i = "MPEG1"

my_colors <- colorRampPalette(c("#FBF4F8","#E5E0ED","#BFC6DD","#8CADCC",
                                "#4E92BA","#1871A8","#085889","#003758"))(100)

fig2 <-  FeaturePlot(datafilt, features = i, reduction = "tsne_harmony",
                     pt.size = 0.1, raster=FALSE) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),            
    axis.ticks = element_blank(),            
    axis.text = element_blank(),             
    axis.title = element_text(colour = "black", size = 15),  
    plot.title = element_text(size = 17, hjust = 0.5),  
    panel.border = element_blank(),          
    legend.position = "right"                
  ) + 
  scale_color_gradientn(colors = my_colors) + 
  labs(x = ' ', y = ' ', title = i)  

fig2

## annotation ======

datafilt.celltypemarkers <- FindAllMarkers(datafilt, only.pos = TRUE, min.pct = 0.15)

datafilt@meta.data <- datafilt@meta.data %>%
  mutate(celltype_minor = recode(seurat_clusters,
                                 "0" = "pDC",
                                 "1" = "cDC2",
                                 "2" = "cDC2",
                                 "3" = "CCL19+DC",
                                 "4" = "cDC2",
                                 "5" = "CD1A+DC",
                                 "6" = "CCL19+DC",
                                 "7" = "cDC1",
                                 "8" = "cDC1",
                                 "9" = "pDC",
                                 "10" = "pDC",
                                 "11" = "CD1A+DC"))

fig1 <- dimplot_new(data = datafilt,
                    reduction = "tsne_harmony",
                    pt.size = 0.1, label = T,
                    group.by = c("celltype_minor"))
fig1

# CCL19+DC cluster 3 vs 7 enrichment

diff_3vs7 <- seurat_diff2(datafilt = datafilt,
                          group.by = "seurat_clusters",
                          group1 = "3",
                          group2 = "7",
                          assay = "RNA",
                          min.pct = 0.15,
                          thres.fc = 0.15)

input <- data.frame(id = diff_3vs7$gene,
                    value = diff_3vs7$logfc)

path = "~/all_genesets.rds"

geneset1 = "GProfiler_GOBP"
geneset2 = "KEGG"

result_GOBP <- gsea_analysis(input = input,
                             source = path,
                             geneset = geneset1,
                             set.min = 3,
                             set.max = 1000)

result_KEGG <- gsea_analysis(input = input,
                             source = path,
                             geneset = geneset2,
                             set.min = 3,
                             set.max = 1000)


## 可视化 ======


## Fig.1e -----

library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

celltype_order <- c("cDC1", "CD1A+DC", "CCL19+DC", "cDC2", "pDC")


gene_sets <- list(
  "Antigen" = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DPB1", "CD74"),
  "Costimulation" = c("CD40", "CD80", "CD83", "CD70", "CD86", "RELB", "TNFSF4", "ICOSLG", "TNFSF9", "TNFSF18"),
  "Immune checkpoints" = c("CD274", "PDCD1LG2", "PVR", "LGALS9", "IDO1", "FAS", "CD200", "HAVCR2", "LILRB1"),
  "Soluble factors" = c("IL10", "IL1B", "CCL2", "CCL4", "CCL5", "XCL1", "CXCL9", "CXCL10", "IFI6", "ISG15"),
  "TLRs and adaptors" = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "MYD88", "TICAM1"),
  "Migration" = c("ICAM1", "MARCKS", "MARCKSL1", "MYO1G", "CCR7")
)

genes_use <- unlist(gene_sets, use.names = FALSE)

# 取平均值

avg_expr <- AverageExpression(
  object = datafilt,
  assays = "RNA",
  features = genes_use,
  group.by = "celltype_minor",
  slot = "data"
)$RNA

# 调整顺序

avg_expr <- avg_expr[genes_use, , drop = FALSE]
keep_celltypes <- celltype_order[celltype_order %in% colnames(avg_expr)]
avg_expr <- avg_expr[, keep_celltypes, drop = FALSE]
mat <- t(avg_expr)
mat <- mat[celltype_order[celltype_order %in% rownames(mat)], , drop = FALSE]

# scale

mat_z <- scale(mat)
mat_z <- as.matrix(mat_z)
mat_z[is.na(mat_z)] <- 0

mat_z[mat_z > 2] <- 2
mat_z[mat_z < -2] <- -2

gene_group <- rep(names(gene_sets), times = sapply(gene_sets, length))
names(gene_group) <- genes_all
gene_group <- factor(gene_group[colnames(mat_z)], levels = names(gene_sets))

# heatmap

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#3969A1", "#F7F7F7", "#A8204C")
)

ht <- Heatmap(
  mat_z,
  name = "Normalized exp.",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 16),
  row_title = NULL,
  
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 11, fontface = "italic"),
  
  column_split = gene_group,
  column_gap = unit(3, "mm"),
  row_gap = unit(0, "mm"),
  
  rect_gp = gpar(col = "#D9D9D9", lwd = 0),
  border = TRUE,
  
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_title_rot = 0,
  
  heatmap_legend_param = list(
    title = "Normalized exp.",
    at = c(-2.2, 0, 2),
    labels = c("-2", "0", "2"),
    direction = "horizontal",
    legend_width = unit(4.5, "cm"),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 11)
  )
)

ht


## Fig.1b ======

datafilt$main_cluster <- dplyr::recode(
  datafilt$celltype_minor,
  "cDC1"     = "DC cluster 1",
  "CD1A+DC"  = "DC cluster 2",
  "CCL19+DC" = "DC cluster 3",
  "cDC2"     = "DC cluster 4",
  "pDC"      = "DC cluster 5"
)

cluster_levels <- c("DC cluster 1", "DC cluster 2", "DC cluster 3", "DC cluster 4", "DC cluster 5")
datafilt$main_cluster <- factor(datafilt$main_cluster, levels = cluster_levels)


cluster_cols <- c(
  "DC cluster 1" = "#36A9B3",  
  "DC cluster 2" = "#E27773",  
  "DC cluster 3" = "#BD9AC9",  
  "DC cluster 4" = "#43AF74",  
  "DC cluster 5" = "#A6A11C")

tsne <- DimPlot(
  object = datafilt,
  reduction = "tsne_harmony",
  group.by = "main_cluster",
  cols = cluster_cols,
  pt.size = 0.7,
  raster = FALSE
) +
  labs(x = "tSNE_1", y = "tSNE_2", color = "Main clusters")

tsne


## Fig.1c ======

Idents(datafilt) <- "main_cluster"

main_markers <- FindAllMarkers(
  datafilt,
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.25
)


# 取定值
top_n <- 40

heatmap_markers <- main_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[["avg_log2FC"]], n = top_n, with_ties = FALSE) %>%
  ungroup()

heatmap_genes <- unique(heatmap_markers$gene)
heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(datafilt)]

## scale
datafilt <- ScaleData(datafilt, features = heatmap_genes, verbose = FALSE)

emb <- as.data.frame(Embeddings(datafilt, "tsne_harmony"))
colnames(emb) <- c("tSNE_1", "tSNE_2")
emb$cell <- rownames(emb)
emb$main_cluster <- datafilt$main_cluster

cell_order <- emb %>%
  mutate(main_cluster = factor(main_cluster, levels = cluster_levels)) %>%
  arrange(main_cluster, tSNE_1, tSNE_2) %>%
  pull(cell)

# DoHeatmap
ht1 <- DoHeatmap(
  object = datafilt,
  features = heatmap_genes,
  cells = cell_order,
  group.by = "main_cluster",
  group.colors = cluster_cols,
  disp.min = -2,
  disp.max = 2,
  raster = FALSE,
  draw.lines = TRUE,
  group.bar = TRUE,
  label = FALSE
) +
  scale_fill_gradientn(
    colors = c("#A64E9F", "black", "#E5D51A"),
    values = rescale(c(-2, 0, 2)),
    limits = c(-2, 2),
    name = "Expression\nz-score\n(log2)"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

ht1

## Fig.1d ======

p1 <- FeaturePlot(
  object = datafilt,
  features = "LILRA4",
  reduction = "tsne_harmony",
  pt.size = 0.35,
  raster = FALSE,
  cols = c("lightgrey", "#4C5AA8")
) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "italic", hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 8)
  )

p1

