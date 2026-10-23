rm(list = ls())
###加载R包
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(DoubletFinder)
library(magrittr)
library(patchwork)
library(harmony)
library(lisi)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(dplyr)
source('lzq_sc_fun.R')


###### pre process
# 读取原始 UMI 矩阵
counts_df <- fread("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt", data.table = FALSE, check.names = FALSE)
names(counts_df) <- counts_df$Index
counts_df <- counts_df[, -1, drop = FALSE]
# 转换为稀疏矩阵
counts <- as(as.matrix(counts_df), "dgCMatrix")
rm(counts_df)

# 创建 Seurat 对象
seu <- CreateSeuratObject(counts = counts, project = "GSE144735", min.cells = 3, min.features = 200)

# 读取 annotation
anno <- fread("GSE144735_processed_KUL3_CRC_10X_annotation.txt", data.table = FALSE, check.names = FALSE)
table(anno$Sample)
# 根据细胞 barcode 匹配 metadata
rownames(anno) <- as.character(anno[[1]])
anno <- anno[, -1, drop = FALSE]
common_cells <- intersect(colnames(seu), rownames(anno))

# 保留能够匹配注释的细胞
seu <- subset(seu, cells = common_cells)
# 添加细胞注释
seu <- AddMetaData(seu, metadata = anno[colnames(seu), , drop = FALSE])
colnames(seu@meta.data)

# 保留T N细胞
table(seu$Class)
seu <- subset(seu, subset = Class %in% c("Normal", "Tumor"))
table(seu$Sample)

# Remove doublet
obj.ls <- SplitObject(seu,split.by = 'Sample') # 分成list，分别做doublet
rm(seu)
obj.ls_rmd <- obj.ls %>% map(lzq_rm_doublet) 
seu <- Reduce(function(x, y) merge(x, y, merge.data = T), obj.ls_rmd)
seu$Cohort <- 'GSE144735'
rm(obj.ls);rm(obj.ls_rmd)

seu@meta.data %<>% mutate(log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA))
save(seu,file = "GSE144735_rmd.rda") #运行时间过长，保存运行后结果

# Quality control -------------------------------------------------
rm(list = ls())
load("GSE144735_rmd.rda")
source('lzq_sc_fun.R')
# 添加线粒体、核糖体、红细胞基因的比例
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = 'percent.mt')
seu <- PercentageFeatureSet(seu, pattern = "^RP[SL]", col.name = 'percent.ribo')
seu <- PercentageFeatureSet(seu, pattern = "^HB[^(P)]", col.name = 'percent.hb')

seu
seu<- subset(seu, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & # Gene numbers
               nCount_RNA >= 1000 & nCount_RNA <= 40000 & # UMI numbers
               log10GenesPerUMI >= 0.8 & # Novelty
               percent.mt <= 20) # 线粒体比例
seu 
save(seu, file = "seu_QC.rda")

# 3. Batch removing -------------------------------------------------------
rm(list = ls())
load("seu_QC.rda")
source('lzq_sc_fun.R')

seu <- SCTransform(seu, vars.to.regress = "percent.mt")
seu <- RunPCA(seu,features = VariableFeatures(object = seu)) 

ElbowPlot(seu, ndims = 50, reduction = "pca")
ggsave("elbowplot.pdf", width=6, height=4,dpi = 500) # 取dims =26

seu <- RunUMAP(seu,reduction = "pca", dims = 1:25, reduction.name = "umap_naive")
seu <- RunTSNE(seu,reduction = "pca", dims = 1:25, reduction.name = "tsne_naive")

seu_2 <- RunHarmony(seu, group.by.vars = 'Sample')

save(seu_2, file="seu_harmony.rda")

###
rm(list = ls())
load("seu_harmony.rda")
seu_2 <- FindNeighbors(seu_2, dims = 1:25, verbose = F) #使用前面计算出来的dims
seu_2 <- FindClusters(seu_2, resolution = 1)

seu_2 <- RunUMAP(seu_2, reduction = "harmony", dims = 1:26, reduction.name = "umap")
seu_2 <- RunTSNE(seu_2, reduction = "harmony", dims = 1:26,reduction.name = "tsne")

colnames(seu_2@meta.data)
DimPlot(seu_2, reduction = "umap", label = F, pt.size = 0.5, group.by = "seurat_clusters")
plot <- DimPlot(seu_2, reduction = "umap", label = F, pt.size = 0.5, group.by = "Cell_type")+
  DimPlot(seu_2, reduction = "umap", label = F, pt.size = 0.5, group.by = "Cell_subtype")
plot  
ggsave(plot, filename = "UMAP.pdf",width = 24,height = 10,limitsize = FALSE)

save(seu_2, file="seu_findclusters.rda")



########### mtGene Score ###########
rm(list = ls())
### 加载R包
library(openxlsx)
library(Seurat)
library(AUCell)
library(GSVA)

## mtGene
mtGene <- c("MT-ATP6", "MT-ATP8", "MT-CO1",  "MT-CO2", "MT-CO3",
            "MT-CYB",  "MT-ND1",  "MT-ND2",  "MT-ND3",  "MT-ND4", 
            "MT-ND4L", "MT-ND5", "MT-ND6")
## 读入seurat
load("seu_findclusters.rda")

## 1.AddModuleScore
gene_set <- mtGene
# 添加模块评分
seu_2 <- AddModuleScore(object = seu_2, 
                        features = list(gene_set), 
                        name = "mtGene_ModuleScore")
## 2.AUC评分
gene_set <- mtGene
geneSet <- list(Module = gene_set)
#表达矩阵
exprMatrix <- seu_2@assays$RNA@data
#计算AUC评分
cell_rankings <- AUCell_buildRankings(exprMatrix) #构建基因排序矩阵
cells_AUC <- AUCell_calcAUC(geneSet, cell_rankings) #计算AUC评分
auc_scores <- as.data.frame(t(getAUC(cells_AUC))) #提取AUC分数 
seu_2 <- AddMetaData(seu_2, metadata = auc_scores[, "Module"], col.name = "mtGene_AUC") 

## 3.GSVA
gene_set <- mtGene
geneSet <- list(Module = gene_set)
exprMatrix <- seu_2@assays$RNA@data
# 计算 GSVA 分值
gsva_scores <- gsva(as.matrix(exprMatrix), geneSet, method = "ssgsea")
head(gsva_scores)
# 添加 GSVA 分数到metadata
gsva_scores_transposed <- t(gsva_scores)
all.equal(rownames(gsva_scores_transposed), colnames(seu_2))  # 确定行名相同

for (module_name in rownames(gsva_scores)) {
  seu_2 <- AddMetaData(seu_2, 
                       metadata = gsva_scores_transposed[, module_name], 
                       col.name = module_name)
}

colnames(seu_2@meta.data)[(ncol(seu_2@meta.data) - 2):ncol(seu_2@meta.data)] <- 
  c("mtGene_ModuleScore", "mtGene_AUC" ,"mtGene_GSVA")

# 保存seurat
save(seu_2, file = "seu_mt.rda")




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


load("seu_mt.rda")
mtGene <- c("MT-ATP6", "MT-ATP8", "MT-CO1",  "MT-CO2", "MT-CO3",
            "MT-CYB",  "MT-ND1",  "MT-ND2",  "MT-ND3",  "MT-ND4", 
            "MT-ND4L", "MT-ND5", "MT-ND6")

table(seu_2$Class)
table(seu_2$Cell_subtype)

###### 1. 统一设置图片风格
group_colors <- c("Normal" = "#66A61E", "Tumor" = "#FF7F0E")
celltype_levels <-  c(
  "Anti-inflammatory", "BEST4+ Enterocytes", "CD19+CD20+ B", "CD4+ T cells",
  "CD8+ T cells", "cDC", "CMS1", "CMS2", "CMS3", "CMS4",
  "Enteric glial cells", "gamma delta T cells", "Goblet cells",
  "IgA+ Plasma", "IgG+ Plasma", "Intermediate", "Lymphatic ECs",
  "Mast cells", "Mature Enterocytes", "Myofibroblasts", "NK cells",
  "Pericytes", "Pro-inflammatory", "Proliferating",
  "Regulatory T cells", "Smooth muscle cells", "SPP1+A", "SPP1+B",
  "Stalk-like ECs", "Stem-like/TA", "Stromal 1", "Stromal 2",
  "Stromal 3", "T follicular helper cells", "T helper 17 cells",
  "Tip-like ECs", "Tuft cells", "Unknown", "Unspecified Plasma"
)
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
meta <- seu_2@meta.data
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
cellname <- rownames(seu_2@meta.data)
expr_matrix <- fread("GSE144735_processed_KUL3_CRC_10X_natural_log_TPM_matrix.txt")
expr <- as.data.frame(expr_matrix[expr_matrix$Index %in% mtGene, ])
rownames(expr) <- make.unique(as.character(expr$Index))
expr <- expr[, -1, drop = FALSE]
stopifnot(all(cellname %in% colnames(expr)))
expr <- as.matrix(expr[, cellname, drop = FALSE])
storage.mode(expr) <- "numeric"
counts_df <- as.data.frame(t(expr))
counts_df$cell_id <- rownames(counts_df)
meta_sample <- data.frame(cell_id = rownames(seu_2@meta.data), Sample = seu_2$Sample, Class = seu_2$Class)
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
meta <- seu_2@meta.data
merged_data <- meta %>%
  filter(Class %in% c("Tumor", "Normal"), !is.na(Sample), !is.na(Cell_subtype)) %>%
  group_by(Sample, Class, Cell_subtype) %>%
  summarise(n_cells = n(),
            mtGene_ModuleScore = mean(mtGene_ModuleScore, na.rm = TRUE),
            mtGene_GSVA = mean(mtGene_GSVA, na.rm = TRUE),
            mtGene_AUC = mean(mtGene_AUC, na.rm = TRUE), .groups = "drop") %>%
  filter(n_cells >= min_cells)
merged_data$Class <- factor(merged_data$Class, levels = c("Normal", "Tumor"))
merged_data$Cell_subtype <- factor(merged_data$Cell_subtype, levels = celltype_levels)
table(merged_data$Cell_subtype, merged_data$Class)

valid_celltypes <- merged_data %>%
  count(Cell_subtype, Class) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  filter(Normal >= min_samples, Tumor >= min_samples) %>%
  pull(Cell_subtype)
valid_score_data <- merged_data %>% filter(Cell_subtype %in% valid_celltypes)

p3 <-
  ggplot(merged_data, aes(Cell_subtype, mtGene_ModuleScore, fill = Class)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  theme_cell +
  labs(title = NULL, x = NULL, y = "mtGene_ModuleScore") +
  stat_compare_means(data = valid_score_data, aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3.5)
p3

p4 <-
  ggplot(merged_data, aes(Cell_subtype, mtGene_GSVA, fill = Class)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  theme_cell +
  labs(title = NULL, x = NULL, y = "mtGene_GSVA") +
  stat_compare_means(data = valid_score_data, aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3.5)
p4

p5 <-
  ggplot(merged_data, aes(Cell_subtype, mtGene_AUC, fill = Class)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.2, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  theme_cell +
  labs(title = NULL, x = NULL, y = "mtGene_AUC") +
  stat_compare_means(data = valid_score_data, aes(group = Class), method = "wilcox.test", label = "p.signif", size = 3.5)
p5

###### D. mtGene_sc
geneset <- intersect(mtGene, rownames(seu_2[["RNA"]]))
expr <- GetAssayData(seu_2, assay = "RNA", layer = "data")
expr <- expr[geneset, , drop = FALSE]
meta <- seu_2@meta.data[colnames(expr), , drop = FALSE]
valid <- !is.na(meta$Sample) & !is.na(meta$Class) & !is.na(meta$Cell_subtype) & meta$Class %in% c("Tumor", "Normal")
expr <- expr[, valid, drop = FALSE]
meta <- meta[valid, , drop = FALSE]
meta$group_id <- paste(meta$Sample, meta$Class, meta$Cell_subtype, sep = "|||")
group_indices <- split(seq_len(ncol(expr)), meta$group_id)

expr_sample_list <- lapply(group_indices, function(idx) {
  if (length(idx) < min_cells) return(NULL)
  avg_expr <- Matrix::rowMeans(expr[, idx, drop = FALSE])
  data.frame(Sample = meta$Sample[idx[1]], Class = meta$Class[idx[1]], celltype = meta$Cell_subtype[idx[1]],
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

ggsave("revise figure_GSE144735_subtype.pdf", plot = p_final,
       width = 36, height = 24, limitsize = FALSE)
