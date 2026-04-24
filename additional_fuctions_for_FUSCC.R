
# dimplot ====================

dimplot_new <- function(data = datafilt,
                        reduction = "umap",
                        pt.size = 1, label = T,
                        group.by = c("seurat_clusters")) {
  
  library(Seurat)
  
  # 设定颜色
  
  info <- data@meta.data
  number <- length(unique(info[,group.by]))
  
  if (number < 17) {
    col <- c("#A6D719", "#176EBF", "#00A8DE", "#AEE0E8",
             "#00A9A3", "#FBD324", "#F28A24", "#A52828",
             "#A37CB7", "#F2D7EE", "#CD6981", "#FBD0C0",
             "#F15E4C", "#ECB2C8", "#B2DBBF", "#CCDAD1")
  } else {
    col <- c("#B8E3EA", "#5CB3DA", "#0070B2", "#FBDD7E", "#F7AE24", "#FF7149", 
             "#F2D7EE", "#A37CB7", "#A231A1", "#ECB2C8", "#E93B8C", "#B91372", 
             "#FF9F99", "#F15E4C", "#DA1735", "#CDE391", "#8BBE53", "#679436", 
             "#98D4C6", "#00A385", "#067D69", "#B2DBBF", "#028090", "#114B5F", 
             "#FBD0C0", "#CD6981", "#A23E48", "#CCDAD1", "#9CAEA9", "#788585")
    col <- colorRampPalette(col)(number)
  }
  
  # 开始画图
  
  DimPlot(data, pt.size = pt.size, label = label, repel = T, 
          raster = TRUE, label.size = 5, reduction = reduction,
          group.by = group.by) + 
    scale_color_manual(values = c(my36colors)) + 
    
    theme_bw() +
    theme(panel.grid = element_blank(), # 删去网格线
          axis.ticks = element_blank(), # 删去刻度线
          axis.text = element_blank(), # 删去刻度标签
          axis.title = element_text(colour = "black", size = 15),
          legend.text = element_text(size = 15),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1)) +   
    
    labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
    guides(colour = guide_legend(title = NULL,
                                 override.aes = list(size = 7.5)))
}


# 两组间差异分析 ====================

seurat_diff2 <- function(datafilt = datafilt,
                         group.by = "seurat_clusters",
                         group1 = "test1",
                         group2 = "test2",
                         assay = "RNA",
                         min.pct = 0.25,
                         thres.fc = 0.25) {
  library(Seurat)
  library(future)
  
  # 设定差异分析的分组
  
  datafilt <- SetIdent(datafilt, value = group.by)
  
  # 并行运算
  
  plan(multisession, workers = 4)
  options(future.globals.maxSize = 300*1024^3)
  
  diff <- FindMarkers(datafilt,
                      test.use = 'LR',
                      ident.1 = group1,
                      ident.2 = group2,
                      assay = assay,
                      slot = 'data', 
                      min.pct = min.pct,
                      logfc.threshold = thres.fc)
  plan(sequential)
  
  # 整理结果
  
  data.frame(gene = rownames(diff),
             logfc = diff$avg_log2FC,
             pctfc = diff$pct.1 - diff$pct.2,
             pct.1 = diff$pct.1, pct.2 = diff$pct.2,
             pvalue = diff$p_val, adjp = diff$p_val_adj)
}


# GSEA分析 ====================

gsea_analysis <- function(input = input,
                          source = path,
                          geneset = NULL,
                          set.min = 5,
                          set.max = 1000) {
  
  library(clusterProfiler)
  library(enrichplot)
  library(gridExtra)
  library(msigdbr)
  options(connectionObserver = NULL)
  
  # 处理注释信息
  
  if (is.null(geneset)) {
    
    if (is.data.frame(source)) {
      
      sig_list = source
      
    } else {
      
      sig_list <- read.table(source, sep = ",", header = F,
                             row.names = 1, na.strings = "")
      
      sig_list <- lapply(rownames(sig_list), function(i){
        subdata <- t(sig_list[i,])[,1]
        subdata <- data.frame(term = i, gene = subdata)
      })
      
      sig_list <- do.call(rbind, sig_list)
      sig_list <- na.omit(sig_list)
      
    }
    
  } else {
    
    all_genesets <- readRDS(source)
    genesets <- reshape2::melt(all_genesets[[geneset]])
    sig_list <- data.frame(term = genesets[,2], gene = genesets[,1])
    
  }
  
  # 整理数据
  
  input <- input[order(input$value, decreasing = T),]
  data <- input$value
  names(data) <- input$id
  
  # 开始进行GSEA分析
  
  kk <- GSEA(data, TERM2GENE = sig_list,
             minGSSize = set.min, maxGSSize = set.max,
             pvalueCutoff = 1, nPermSimple = 10000, eps = 0)
  
  # 结果整理
  
  result <- data.frame(id = kk$Description, 
                       ES = kk$enrichmentScore,
                       NES = kk$NES, 
                       pvalue = kk$pvalue, 
                       FDR = kk$p.adjust, 
                       setsize = kk$setSize, 
                       leading_prop = kk$leading_edge,
                       leading_gene = kk$core_enrichment)
  
  result$leading_prop <- do.call(rbind, strsplit(result$leading_prop, ','))[,1]
  result$leading_prop <- do.call(rbind, strsplit(result$leading_prop, '='))[,2]
  
  # 结果输出
  
  result <- result[order(result$NES, decreasing = T),]
  result
}


