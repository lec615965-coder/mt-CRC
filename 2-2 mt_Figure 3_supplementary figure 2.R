################## figure 3 ##################
rm(list = ls())
### 加载R包
library(openxlsx)
library(Seurat)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggthemr)

### 加载数据
## 读入mtGene
mtGene <- read.xlsx("data/mtGene/mitochondria gene_geneID.xlsx") 
mtGene <- mtGene$Gene.Name
## 读入seurat
load("data/GSE132465_final.rda")

###### A. pseudobulk score 
meta <- GSE132465_2@meta.data
meta$cell_id <- rownames(meta)
merged_data <- data.frame(meta$cell_id, meta$Sample,
                          meta$mtGene_ModuleScore, meta$mtGene_AUC, meta$mtGene_GSVA)
colnames(merged_data) <- c("cell_id", "Sample",
                           "mtGene_ModuleScore", "mtGene_AUC", "mtGene_GSVA")

summarised_data <- merged_data %>%
  group_by(Sample) %>%  #按细胞来源的Sample汇总（模拟bulk）→pseudobulk
  summarise(across(everything(), mean, na.rm = TRUE))

summarised_data[1:6, 1:3]
rownames <- summarised_data$Sample
summarised_data <- summarised_data[,-c(1,2)]
row.names(summarised_data) <- rownames
summarised_data[1:6, 1:3]
summarised_data <- t(summarised_data)
pseudobulk_data <- as.data.frame(summarised_data)

gene_expr <- pseudobulk_data
gene_expr_long <- gene_expr %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = ends_with("T") | ends_with("N"), 
               names_to = "sample", 
               values_to = "score")
gene_expr_long$score <- as.numeric(gene_expr_long$score)
gene_expr_long <- gene_expr_long %>%
  mutate(group = ifelse(grepl("T$", sample), "Tumor", 
                        ifelse(grepl("N$", sample), "Normal", NA)))

gene_expr_long$gene <- factor(gene_expr_long$gene,
                              levels = c("mtGene_ModuleScore", 
                                         "mtGene_GSVA", 
                                         "mtGene_AUC"))  # 这里改成你想要的顺序
p1<-
  ggplot(gene_expr_long, aes(gene, score, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#66A61E", "#FF7F0E")) +
  theme_bw() +
  labs(title = NULL, x = NULL, y = "Mitochondria gene score") +
  theme(
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),  # 调整y轴标题
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1, vjust = 1, 
                               colour = "black", face = "plain"),
    axis.text.y = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")
p1

###### B.pseudobulk mtgene
cellname <- rownames(GSE132465_2@meta.data)
expr_matrix <- fread("data/CRC---GSE132465/GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt")
expr_matrix[1:3,1:3]

expr <- expr_matrix[(expr_matrix$Index) %in% mtGene, ]
expr[1:3,1:3]
expr <- as.data.frame(expr)
rownames(expr) <- mtGene
expr <- expr[,-1]
expr[1:3,1:3]
expr <- expr[,cellname]

sc_data <- GSE132465_2
sc_data@meta.data$cell_id <- rownames(sc_data@meta.data)
counts_df <- expr %>%
  t() %>% 
  as.data.frame()
counts_df$cell_id <- rownames(sc_data@meta.data)
meta <- data.frame(sc_data@meta.data$cell_id, sc_data@meta.data$Sample)
colnames(meta) <- c( "cell_id", "Sample")
merged_data <- merge(counts_df,meta, by = "cell_id")

summarised_data <- merged_data %>%
  group_by(Sample) %>%  #按细胞来源的Sample汇总（模拟bulk）→pseudobulk
  summarise(across(everything(), mean, na.rm = TRUE))

summarised_data[1:6, 1:6]
rownames <- summarised_data$Sample
summarised_data <- summarised_data[,-c(1,2)]
row.names(summarised_data) <- rownames
summarised_data[1:6, 1:6]
summarised_data <- t(summarised_data)
pseudobulk_data <- as.data.frame(summarised_data)

gene_expr <- pseudobulk_data
gene_expr_long <- gene_expr %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = ends_with("T") | ends_with("N"), 
               names_to = "sample", 
               values_to = "expression")
gene_expr_long$expression <- as.numeric(gene_expr_long$expression)
gene_expr_long <- gene_expr_long %>%
  mutate(group = ifelse(grepl("T$", sample), "Tumor", 
                        ifelse(grepl("N$", sample), "Normal", NA)))

p2<-
ggplot(gene_expr_long,aes(gene,expression,fill = group))+
  geom_boxplot()+
  scale_fill_manual(values = c("#66A61E", "#FF7F0E"))+
  theme_bw()+
  labs(title=NULL,x=NULL,y="Mitochondria gene expression (logTPM)")+
  theme(
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),  # 调整y轴标题
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1, vjust = 1, 
                               colour = "black", face = "plain"),
    axis.text.y = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.signif")
p2

grid.arrange(p1, p2, ncol=2, widths=c(1,3.5 ))


pdf("output/Figure 3/A+B.mtGene Expression in Tumor vs Normal in CRC (logTPM).pdf",
    width = 12, height = 3)  
   grid.arrange(p1, p2, ncol = 2, widths = c(1, 3.5))
dev.off()

###### C. mtGene score_sc
meta <- GSE132465_2@meta.data
gene_expr <- meta %>%
  filter(Class %in% c("Tumor", "Normal"), !is.na(Sample), !is.na(bulk_labels_3)) %>%
  group_by(Sample, Class, bulk_labels_3) %>%
  summarise(n_cells = n(),
            mtGene_ModuleScore = mean(mtGene_ModuleScore, na.rm = TRUE),
            mtGene_GSVA = mean(mtGene_GSVA, na.rm = TRUE),
            mtGene_AUC = mean(mtGene_AUC, na.rm = TRUE),
            .groups = "drop") %>%
  filter(n_cells >= 10)

gene_expr$Class <- factor(gene_expr$Class, levels = c("Normal", "Tumor"))

# 至少每组3个样本才进行统计检验
valid_celltypes <- gene_expr %>%
  count(bulk_labels_3, Class) %>%
  pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
  filter(Normal >= 3, Tumor >= 3) %>%
  pull(bulk_labels_3)

# ModuleScore
p3 <-
  ggplot(gene_expr, aes(bulk_labels_3, mtGene_ModuleScore, fill = Class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.7) +
  scale_fill_manual(values = c("Normal" = "#66A61E", "Tumor" = "#FF7F0E")) +
  theme_bw() +
  labs(title = NULL, x = NULL, y = "mtGene_ModuleScore") +
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 0.9),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.3),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  stat_compare_means(data = gene_expr %>% filter(bulk_labels_3 %in% valid_celltypes),
                     aes(group = Class), method = "wilcox.test", label = "p.signif")

# GSVA
p4 <-
  ggplot(gene_expr, aes(bulk_labels_3, mtGene_GSVA, fill = Class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.7) +
  scale_fill_manual(values = c("Normal" = "#66A61E", "Tumor" = "#FF7F0E")) +
  theme_bw() +
  labs(title = NULL, x = NULL, y = "mtGene_GSVA") +
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 0.9),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.3),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  stat_compare_means(data = gene_expr %>% filter(bulk_labels_3 %in% valid_celltypes),
                     aes(group = Class), method = "wilcox.test", label = "p.signif")

# AUC
p5 <-
  ggplot(gene_expr, aes(bulk_labels_3, mtGene_AUC, fill = Class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.7) +
  scale_fill_manual(values = c("Normal" = "#66A61E", "Tumor" = "#FF7F0E")) +
  theme_bw() +
  labs(title = NULL, x = NULL, y = "mtGene_AUC") +
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust = 0.9),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.3),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  stat_compare_means(data = gene_expr %>% filter(bulk_labels_3 %in% valid_celltypes),
                     aes(group = Class), method = "wilcox.test", label = "p.signif")

# 合并输出
pdf("output/Figure 3/C.mtGene score_sc_sample.pdf", width = 7, height = 11)
grid.arrange(p3, p4, p5, ncol = 1)
dev.off()

###### D. Featureplot
plot.df <- data.frame(GSE132465_2@meta.data, GSE132465_2@reductions$umap_naive@cell.embeddings)

ggplot() + 
  geom_point(data=plot.df, aes(x=UMAP_1,y=UMAP_2,colour=mtGene_ModuleScore), size = 0.5) +
  ggtitle('Mitochondria gene score')+
  facet_wrap(~Class, nrow = 1 ) +
  scale_color_gradientn(
    colours = c("#66A61E", "#98DF8A", "#FFBB78", "#FF7F0E")
  ) +
  theme(legend.position="right",
        legend.title=element_text(size = 10,hjust = 0,color = '#747171',face = 'bold'),
        legend.key.height = unit(18,'pt'),
        legend.key.width = unit(12,'pt'),
        #legend.box.background =element_rect(color="black",size=1.5,linetype="solid"),
        legend.text = element_text(hjust = 0 ,size = 10,color = '#747171'),
        axis.title = element_text(size = 14,face = 'italic',colour = 'black'),
        axis.ticks = element_blank(),
        #axis.text = element_text(colour = 'black',size = 14),
        axis.text = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 15,face = 'bold',colour = 'grey40'),
        #panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
        plot.title = element_text(size = 15,hjust = 0.5,color = 'black'))

ggsave("output/Figure 3/C.UMAP of mtScore in T&N.pdf", width = 10,height = 4)

###### E.Dotplot
dt <- GSE132465_2
dt$Tissue_CellType <- paste(dt$bulk_labels_3, dt$Class,  sep = "_")
dt$Tissue_CellType <- sub("Tumor", "T", dt$Tissue_CellType)
dt$Tissue_CellType <- sub("Normal", "N", dt$Tissue_CellType)
metadata <- dt@meta.data

p6 <-
  DotPlot(dt, features = mtGene, group.by = "Tissue_CellType")+theme_bw()+
  ggtitle("Expression of mtGene in Different Cell Types in T&N") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#66A61E", "#98DF8A", "#FFBB78", "#FF7F0E"))+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  coord_flip() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  # 设置 x 轴标签倾斜 45 度
p6

ggsave("output/Figure 3/E.Dotplot of mtGene in T&N.pdf", plot = p6, width = 11,height = 4)



###### supplementary figure 2
geneset <- mtGene
plots <- list()

for (gene in geneset) {
  
  gene_expr <- FetchData(GSE132465_2,
                         vars = c(gene, "Sample", "bulk_labels_3", "Class"),
                         assay = "RNA", layer = "data")
  
  # 按Sample + celltype聚合
  gene_expr_long <- gene_expr %>%
    filter(Class %in% c("Tumor", "Normal"), !is.na(Sample), !is.na(bulk_labels_3)) %>%
    group_by(Sample, Class, bulk_labels_3) %>%
    summarise(n_cells = n(),
              expression = mean(.data[[gene]], na.rm = TRUE),
              .groups = "drop") %>%
    filter(n_cells >= 10, is.finite(expression))
  
  gene_expr_long$Class <- factor(gene_expr_long$Class, levels = c("Normal", "Tumor"))
  
  # 每组至少3个样本才进行检验
  valid_celltypes <- gene_expr_long %>%
    count(bulk_labels_3, Class) %>%
    pivot_wider(names_from = Class, values_from = n, values_fill = 0) %>%
    filter(Normal >= 3, Tumor >= 3) %>%
    pull(bulk_labels_3)
  
  p2 <-
    ggplot(gene_expr_long, aes(bulk_labels_3, expression, fill = Class)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1, alpha = 0.7) +
    scale_fill_manual(values = c("Normal" = "#66A61E", "Tumor" = "#FF7F0E")) +
    theme_bw() +
    labs(title = paste(gene, "expression of normal and tumor samples in CRC"),
         x = NULL, y = "Mean normalized expression") +
    theme(axis.text = element_text(size = 12),
          axis.title.y = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 30, hjust = 0.9),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.3),
          legend.title = element_text(size = 12, face = "bold")) +
    stat_compare_means(data = gene_expr_long %>% filter(bulk_labels_3 %in% valid_celltypes),
                       aes(group = Class), method = "wilcox.test", label = "p.signif")
  
  plots[[gene]] <- p2
}

ggsave("output/Figure 3/mtGene_expr_boxplot_T&N_sample.pdf",
       plot = grid.arrange(grobs = plots, ncol = 4), width = 24, height = 16)

