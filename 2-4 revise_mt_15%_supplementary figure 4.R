########### seurat ###########
# 1. Quality control -------------------------------------------------
rm(list = ls())
load("data/GSE132465_rmd.rda")
# 添加线粒体、核糖体、红细胞基因的比例
GSE132465 <- PercentageFeatureSet(GSE132465, pattern = "^MT-", col.name = 'percent.mt')
GSE132465 <- PercentageFeatureSet(GSE132465, pattern = "^RP[SL]", col.name = 'percent.ribo')
GSE132465 <- PercentageFeatureSet(GSE132465, pattern = "^HB[^(P)]", col.name = 'percent.hb')

# Remove cells with low quality or potential doublets
GSE132465<- subset(GSE132465, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & # Gene numbers
                     nCount_RNA >= 1000 & nCount_RNA <= 40000 & # UMI numbers
                     log10GenesPerUMI >= 0.8 & # Novelty
                     percent.mt <= 15) # 线粒体比例

# 2. Batch removing -------------------------------------------------------
GSE132465 <- SCTransform(GSE132465, vars.to.regress = "percent.mt")
GSE132465 <- RunPCA(GSE132465,features = VariableFeatures(object = GSE132465)) 
ElbowPlot(GSE132465, ndims = 50, reduction = "pca")
GSE132465_2 <- RunHarmony(GSE132465, group.by.vars = 'orig.ident')

# 3. Find cluster -------------------------------------------- 
GSE132465_2 <- FindNeighbors(GSE132465_2, dims = 1:26, verbose = F) #使用原来的dims、resolution
GSE132465_2 <- FindClusters(GSE132465_2, resolution = 1)

GSE132465_2 <- RunUMAP(GSE132465_2, reduction = "harmony", dims = 1:26, reduction.name = "umap")
GSE132465_2 <- RunTSNE(GSE132465_2, reduction = "harmony", dims = 1:26, reduction.name = "tsne")

plot <- DimPlot(GSE132465_2, reduction = "umap", label = T, pt.size = 0.5, group.by = "seurat_clusters")+
  DimPlot(GSE132465_2, reduction = "umap", label = T, pt.size = 0.5, group.by = "Cell_type")
plot  
ggsave(plot, filename = "output/revise_mt_10_15/UMAP_mt 15.pdf",width = 12,height = 6,limitsize = FALSE)

# 4. annotation -------------------------------------------- 
degs <- FindAllMarkers(
  GSE132465_2, 
  only.pos = TRUE,              # 只考虑上调基因
  logfc.threshold = 0.25,       # 设置 log fold change 的阈值
  min.pct = 0.1)

top10 <- degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(degs, "data/GSE132465_markers_mt_15.csv")

Marker <- c(
  "PTPRC", # immune
  "CD3D", "CD3E",  # T
  "CD8A", "CD8B",  "NKG7", "PRF1",# CTL
  "IL7R", # Th / Memory T
  "FOXP3", # Treg
  "CCR7", # Naive T
  "CD79A", "CD79B", "MS4A1", # B
  "MZB1", # plasma
  "CD68", # M1 macrophage
  "CD163",# M2 macrophage
  "IRF7", "HLA-DRA",  # DC
  "KIT",  # mast cell
  "CDH5", "VWF",  # endothelial cell
  "COL1A1",  "COL6A1",  "THY1", # fibroblast
  "EPCAM", "KRT18",  # # malignant epithelial 
  "SPINK4", "TFF3", # goblet
  "OLFM4", "ASCL2", # ISC
  "MKI67",  "TOP2A", # TAC
  "FABP1",  "CEACAM1" # enterocyte
)

p <- 
  DotPlot(GSE132465_2, group.by = "seurat_clusters", features = Marker) +
  theme_bw() +
  coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 13, face = "plain"), 
        legend.title = element_text(size = 13, face = "bold")) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colours = c("#1F77B4","#9EDAE5", "#FF9896","#D62728")) + 
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3))

ggsave(filename = "output/revise_mt_10_15/celltype marker_dotplot_mt 15.pdf", plot = p,
       width = 16,height = 6,limitsize = FALSE)

# add annotation
# 定义 Cluster 0–34 对应的细胞类型
celltype_map <- c(
  "0"  = "memory T",
  "1"  = "memory T",
  "2"  = "ISC",
  "3"  = "naive B",
  "4"  = "macrophage",
  "5"  = "CTL",
  "6"  = "macrophage",
  "7"  = "naive T",
  "8"  = "memory T",
  "9"  = "goblet",
  "10" = "CTL",
  "11" = "enterocyte",
  "12" = "malignant epithelial",
  "13" = "memory T",
  "14" = "memory T",
  "15" = "malignant epithelial",
  "16" = "CTL",
  "17" = "endothelial cell",
  "18" = "fibroblast",
  "19" = "goblet",
  "20" = "fibroblast",
  "21" = "malignant epithelial",
  "22" = "CTL",
  "23" = "fibroblast",
  "24" = "ISC",
  "25" = "fibroblast",
  "26" = "plasma",
  "27" = "DC",
  "28" = "enterocyte",
  "29" = "naive B",
  "30" = "fibroblast",
  "31" = "DC",
  "32" = "fibroblast",
  "33" = "mast",
  "34" = "TAC")

# 写入 metadata
GSE132465_2$celltype <- unname(celltype_map[clusters])
# 查看各细胞类型数量
table(GSE132465_2$celltype)
# 绘制 UMAP
DimPlot(GSE132465_2, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)

# 保存对象
save(GSE132465_2, file = "data/GSE132465_anno_mt 15.rda")



########### mtGene Score ###########
rm(list = ls())
### 加载R包
library(openxlsx)
library(Seurat)
library(AUCell)
library(GSVA)

## 读入mtGene
mtGene <- read.xlsx("data/mtGene/mitochondria gene_geneID.xlsx") 
mtGene <- mtGene$Gene.Name
## 读入seurat
load("data/GSE132465_anno_mt 15.rda")

## 1.AddModuleScore
gene_set <- mtGene
# 添加模块评分
GSE132465_2 <- AddModuleScore(object = GSE132465_2, 
                              features = list(gene_set), 
                              name = "mtGene_ModuleScore")
## 2.AUC评分
gene_set <- mtGene
geneSet <- list(Module = gene_set)
#表达矩阵
exprMatrix <- GSE132465_2@assays$RNA@data
#计算AUC评分
cell_rankings <- AUCell_buildRankings(exprMatrix) #构建基因排序矩阵
cells_AUC <- AUCell_calcAUC(geneSet, cell_rankings) #计算AUC评分
auc_scores <- as.data.frame(t(getAUC(cells_AUC))) #提取AUC分数 
GSE132465_2 <- AddMetaData(GSE132465_2, metadata = auc_scores[, "Module"], col.name = "mtGene_AUC") 

## 3.GSVA
gene_set <- mtGene
geneSet <- list(Module = gene_set)
exprMatrix <- GSE132465_2@assays$RNA@data
# 计算 GSVA 分值
gsva_scores <- gsva(as.matrix(exprMatrix), geneSet, method = "ssgsea")
head(gsva_scores)
# 添加 GSVA 分数到metadata
gsva_scores_transposed <- t(gsva_scores)
all.equal(rownames(gsva_scores_transposed), colnames(GSE132465_2))  # 确定行名相同

for (module_name in rownames(gsva_scores)) {
  GSE132465_2 <- AddMetaData(GSE132465_2, 
                             metadata = gsva_scores_transposed[, module_name], 
                             col.name = module_name)
}

colnames(GSE132465_2@meta.data)[(ncol(GSE132465_2@meta.data) - 2):ncol(GSE132465_2@meta.data)] <- 
  c("mtGene_ModuleScore", "mtGene_AUC" ,"mtGene_GSVA")

# 保存seurat
save(GSE132465_2, file = "data/GSE132465_mt 15.rda")



########### 绘图 ###########
rm(list = ls())

###### 0. 加载R包与数据
library(openxlsx)
library(Seurat)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(Matrix)

mtGene <- read.xlsx("data/mtGene/mitochondria gene_geneID.xlsx")
mtGene <- unique(na.omit(mtGene$Gene.Name))
load("data/GSE132465_mt 15.rda")

###### 1. 统一设置图片风格
group_colors <- c("Normal" = "#66A61E", "Tumor" = "#FF7F0E")
celltype_levels <- c("CTL", "DC", "endothelial cell", "enterocyte", "fibroblast", "goblet", "ISC", "macrophage",
                     "malignant epithelial", "mast", "memory T", "naive B", "naive T", "plasma", "TAC")

theme_fig <- theme_bw(base_size = 14) +
  theme(axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, angle = 30, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 13),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        plot.margin = margin(8, 8, 8, 8))

theme_top <- theme_fig +
  theme(axis.text.x = element_text(size = 13, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"))

theme_cell <- theme_fig +
  theme(axis.text.x = element_text(size = 11, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.position = "top")

theme_gene <- theme_fig +
  theme(axis.text.x = element_text(size = 10, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(6, 6, 6, 6))

###### 2. 数据检查
meta <- GSE132465_2@meta.data
table(meta$Class, useNA = "ifany")
min_cells <- 10
min_samples <- 3

###### A. pseudobulk score
merged_data <- meta %>%
  filter(Class %in% c("Normal", "Tumor"), !is.na(Sample)) %>%
  group_by(Sample, Class) %>%
  summarise(mtGene_ModuleScore = mean(mtGene_ModuleScore, na.rm = TRUE),
            mtGene_GSVA = mean(mtGene_GSVA, na.rm = TRUE),
            mtGene_AUC = mean(mtGene_AUC, na.rm = TRUE), .groups = "drop")

gene_expr_long <- merged_data %>%
  pivot_longer(cols = c(mtGene_ModuleScore, mtGene_GSVA, mtGene_AUC), names_to = "gene", values_to = "score") %>%
  filter(is.finite(score))
gene_expr_long$gene <- factor(gene_expr_long$gene, levels = c("mtGene_ModuleScore", "mtGene_GSVA", "mtGene_AUC"))
gene_expr_long$Class <- factor(gene_expr_long$Class, levels = c("Normal", "Tumor"))

p1 <-
  ggplot(gene_expr_long, aes(gene, score, fill = Class)) +
  geom_boxplot(width = 0.65, linewidth = 0.5) +
  scale_fill_manual(values = group_colors) +
  theme_top +
  labs(title = NULL, x = NULL, y = "Mitochondria gene score") +
  stat_compare_means(aes(group = Class), method = "wilcox.test", label = "p.signif", size = 4)
p1

###### B. pseudobulk mtGene expression
cellname <- rownames(GSE132465_2@meta.data)
expr_matrix <- fread("data/CRC---GSE132465/GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt")
expr <- as.data.frame(expr_matrix[expr_matrix$Index %in% mtGene, ])
rownames(expr) <- make.unique(as.character(expr$Index))
expr <- expr[, -1, drop = FALSE]
stopifnot(all(cellname %in% colnames(expr)))
expr <- as.matrix(expr[, cellname, drop = FALSE])
storage.mode(expr) <- "numeric"
counts_df <- as.data.frame(t(expr))
counts_df$cell_id <- rownames(counts_df)
meta_sample <- data.frame(cell_id = rownames(GSE132465_2@meta.data), Sample = GSE132465_2$Sample, Class = GSE132465_2$Class)
merged_data <- merge(counts_df, meta_sample, by = "cell_id")

summarised_data <- merged_data %>%
  filter(Class %in% c("Normal", "Tumor")) %>%
  group_by(Sample, Class) %>%
  summarise(across(all_of(rownames(expr)), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
gene_expr_long <- summarised_data %>%
  pivot_longer(cols = all_of(rownames(expr)), names_to = "gene", values_to = "expression") %>%
  filter(is.finite(expression))
gene_expr_long$gene <- factor(gene_expr_long$gene, levels = rownames(expr))
gene_expr_long$Class <- factor(gene_expr_long$Class, levels = c("Normal", "Tumor"))

p2 <-
  ggplot(gene_expr_long, aes(gene, expression, fill = Class)) +
  geom_boxplot(width = 0.65, linewidth = 0.5) +
  scale_fill_manual(values = group_colors) +
  theme_top +
  labs(title = NULL, x = NULL, y = "Mitochondria gene expression (logTPM)") +
  stat_compare_means(aes(group = Class), method = "wilcox.test", label = "p.signif", size = 4)
p2

###### C. mtGene score_sc
meta <- GSE132465_2@meta.data
merged_data <- meta %>%
  filter(Class %in% c("Tumor", "Normal"), !is.na(Sample), !is.na(celltype)) %>%
  group_by(Sample, Class, celltype) %>%
  summarise(n_cells = n(),
            mtGene_ModuleScore = mean(mtGene_ModuleScore, na.rm = TRUE),
            mtGene_GSVA = mean(mtGene_GSVA, na.rm = TRUE),
            mtGene_AUC = mean(mtGene_AUC, na.rm = TRUE), .groups = "drop") %>%
  filter(n_cells >= min_cells)
merged_data$Class <- factor(merged_data$Class, levels = c("Normal", "Tumor"))
merged_data$celltype <- factor(merged_data$celltype, levels = celltype_levels)
table(merged_data$celltype, merged_data$Class)

valid_celltypes <- merged_data %>%
  count(celltype, Class) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  filter(Normal >= min_samples, Tumor >= min_samples) %>%
  pull(celltype)
valid_score_data <- merged_data %>% filter(celltype %in% valid_celltypes)

p3 <-
  ggplot(merged_data, aes(celltype, mtGene_ModuleScore, fill = Class)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  theme_cell +
  labs(title = NULL, x = NULL, y = "mtGene_ModuleScore") +
  stat_compare_means(data = valid_score_data, aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3.5)
p3

p4 <-
  ggplot(merged_data, aes(celltype, mtGene_GSVA, fill = Class)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  theme_cell +
  labs(title = NULL, x = NULL, y = "mtGene_GSVA") +
  stat_compare_means(data = valid_score_data, aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3.5)
p4

p5 <-
  ggplot(merged_data, aes(celltype, mtGene_AUC, fill = Class)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  theme_cell +
  labs(title = NULL, x = NULL, y = "mtGene_AUC") +
  stat_compare_means(data = valid_score_data, aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3.5)
p5

###### D. mtGene_sc
geneset <- intersect(mtGene, rownames(GSE132465_2[["RNA"]]))
expr <- GetAssayData(GSE132465_2, assay = "RNA", layer = "data")
expr <- expr[geneset, , drop = FALSE]
meta <- GSE132465_2@meta.data[colnames(expr), , drop = FALSE]
valid <- !is.na(meta$Sample) & !is.na(meta$Class) & !is.na(meta$celltype) & meta$Class %in% c("Tumor", "Normal")
expr <- expr[, valid, drop = FALSE]
meta <- meta[valid, , drop = FALSE]
meta$group_id <- paste(meta$Sample, meta$Class, meta$celltype, sep = "|||")
group_indices <- split(seq_len(ncol(expr)), meta$group_id)

expr_sample_list <- lapply(group_indices, function(idx) {
  if (length(idx) < min_cells) return(NULL)
  avg_expr <- Matrix::rowMeans(expr[, idx, drop = FALSE])
  data.frame(Sample = meta$Sample[idx[1]], Class = meta$Class[idx[1]], celltype = meta$celltype[idx[1]],
             n_cells = length(idx), gene = names(avg_expr), expression = as.numeric(avg_expr))
})

gene_expr_long <- dplyr::bind_rows(expr_sample_list) %>%
  filter(Class %in% c("Tumor", "Normal"), is.finite(expression)) %>%
  mutate(Class = factor(Class, levels = c("Normal", "Tumor")),
         celltype = factor(celltype, levels = celltype_levels))
# 每个基因对应的样本数相同，检查时取其中一个基因即可
table(gene_expr_long$celltype[gene_expr_long$gene == geneset[1]],
      gene_expr_long$Class[gene_expr_long$gene == geneset[1]])

plots <- list()
for (gene in geneset) {
  gene_expr <- gene_expr_long %>% filter(.data$gene == .env$gene)
  valid_types <- gene_expr %>%
    count(celltype, Class) %>%
    pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
    filter(Normal >= min_samples, Tumor >= min_samples) %>%
    pull(celltype)
  
  p6 <-
    ggplot(gene_expr, aes(celltype, expression, fill = Class)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 0.8, alpha = 0.7) +
    scale_fill_manual(values = group_colors) +
    theme_gene +
    labs(title = gene, x = NULL, y = "Mean normalized expression") +
    stat_compare_means(data = gene_expr %>% filter(celltype %in% valid_types),
                       aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3)
  plots[[gene]] <- p6
}

###### E. 合并ABCD图
pAB <- ggarrange(p1, p2, ncol = 2, widths = c(1, 3), labels = c("A", "B"),
                 font.label = list(size = 20, face = "bold", color = "black"))
pC <- ggarrange(p3, p4, p5, ncol = 3, widths = c(1, 1, 1), labels = c("C", "", ""),
                font.label = list(size = 20, face = "bold", color = "black"))
p6_combined <- ggarrange(plotlist = plots, ncol = 4, nrow = ceiling(length(plots) / 4), align = "hv")
pD <- ggarrange(p6_combined, ncol = 1, labels = "D",
                font.label = list(size = 20, face = "bold", color = "black"))
p_final <- ggarrange(pAB, pC, pD, ncol = 1, nrow = 3, heights = c(1, 1.5, 4), align = "v")

ggsave("output/revise_mt_10_15/revise figure1_mt15_CRC_final.pdf", plot = p_final,
       width = 24, height = 24, limitsize = FALSE)

