rm(list=ls())
# 加载R包
library(limma)
library(edgeR)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)

### 1.数据准备
#读入expr
#分组信息
expr <- read.csv("data/cohort/TCGA_COAD_expr.csv", row.names = 1)
colnames(expr) <- gsub("\\.", "-", colnames(expr))
TCGA_RS <- read.csv("data/revise mime/TCGA_RS.csv")
RS <- TCGA_RS[,c(2,6),drop=F]

# 确保样本顺序匹配
common_samples <- intersect(colnames(expr), RS$ID)
expr <- expr[, common_samples]
RS <- RS[match(common_samples, RS$ID), ]
# 确认分组因子
group <- factor(RS$Group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

## 差异分析 (High vs Low)
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, adjust = "fdr", number = Inf)
head(deg)

# 保存差异分析结果
write.csv(deg, "data/revise group/TCGA_COAD_DEG_limma.csv")


##提取上调与下调基因（1.5倍）
DEG_up <- subset(deg, logFC > 0.585 & adj.P.Val < 0.05)
DEG_down <- subset(deg, logFC < -0.585 & adj.P.Val < 0.05)
cat("上调基因数量：", nrow(DEG_up), "\n")
cat("下调基因数量：", nrow(DEG_down), "\n")

# 注明DEG_up and DEG_down
DEG_up<- DEG_up%>%
  mutate(source = "Upregulated")
DEG_down<- DEG_down%>%
  mutate(source = "Downregulated")
all_genes <- bind_rows(DEG_up,DEG_down)

# 保存差异分析结果
write.csv(DEG_up, file = "data/revise group/DEG_up.csv")
write.csv(DEG_down, file = "data/revise group/DEG_down.csv")
write.csv(all_genes, file = "data/revise group/up and down_genes.csv")

## 火山图
# 创建一个新的列以表示显著性
deg$significant <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.585, 
                          ifelse(deg$logFC > 0.585, 'Upregulated', 'Downregulated'), 
                          'Not Significant')
deg$gene_name <- rownames(deg)
# 选择前10个上调和下调的显著基因
top_upregulated <- head(deg[deg$adj.P.Val < 0.05 & deg$logFC > 0.585, ][order(deg$logFC[deg$adj.P.Val < 0.05 & deg$logFC > 0.585], decreasing = TRUE), ], 10)
top_downregulated <- head(deg[deg$adj.P.Val < 0.05 & deg$logFC < 0.585, ][order(deg$logFC[deg$adj.P.Val < 0.05 & deg$logFC < 0.585]), ], 10)
# 提取 top 基因名称
top_genes <- c(top_upregulated$gene_name, top_downregulated$gene_name)
top_genes
# 绘制火山图
volcano_plot <- ggplot(deg, aes(logFC, -log10(adj.P.Val)))+   #读取差异表达结果，X轴为logFC,y轴为-log10(P.Value)
  geom_point(aes(col=significant))+    #设置点数据的来源
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+  #这里可以对途中，上调，不变和下调的点的颜色进行设置
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+  #设置x轴和y轴的标签
  geom_vline(xintercept=c(-0.5,0.5), colour="grey", linetype="dashed")+  #设置x轴的分界线
  geom_hline(yintercept = -log10(0.05),colour="grey", linetype="dashed")+   # y轴分界线
  geom_text_repel(data = subset(deg, gene_name %in% top_genes),  # 仅标注 top_genes
                  aes(label = gene_name), 
                  size = 2, 
                  show.legend = FALSE)+ 
  xlim(-1.5, 1.5) +  # 设置 x 轴范围
  ylim(0, 5)    # 设置 y 轴范围
volcano_plot
ggsave("output/revise group/volcano_plot.pdf", plot = volcano_plot, width = 5, height = 4)


### GO/KEGG
# 转换基因ID
up_genes <- bitr(rownames(DEG_up), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_genes <- bitr(rownames(DEG_down), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO
ego_up <- enrichGO(gene = up_genes$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
ego_down <- enrichGO(gene = down_genes$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
ego_all <- enrichGO(gene = c(up_genes$ENTREZID,down_genes$ENTREZID),
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# KEGG
ekegg_up <- enrichKEGG(gene = up_genes$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
ekegg_down <- enrichKEGG(gene = down_genes$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)
ekegg_all <- enrichKEGG(gene = c(up_genes$ENTREZID,down_genes$ENTREZID),
                        organism = 'hsa',
                        pvalueCutoff = 0.05)


# GO
pdf(file="output/revise group/A.up and down_GO.pdf",width = 13,height = 5)
dotplot(ego_up, showCategory=10, title="GO (Up-regulated genes)")+
  dotplot(ego_down, showCategory=10, title="GO (Down-regulated genes)")
dev.off()

# KEGG
pdf(file="output/revise group/A.up and down_KEGG.pdf",width = 13,height = 5)
dotplot(ekegg_up, showCategory=10, title="KEGG (Up-regulated genes)")+
  dotplot(ekegg_down, showCategory=10, title="KEGG (Down-regulated genes)")
dev.off()




### 差异分析 GSEA
### 1.数据准备
# 确保样本匹配
common_samples <- intersect(colnames(expr), RS$ID)
expr <- expr[, common_samples]
RS <- RS[match(common_samples, RS$ID), ]
# 构建分组因子
group <- factor(RS$Group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

## 差异分析
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, adjust = "fdr", number = Inf)

# 提取排序好的基因列表
gene_list <- deg$t  # 使用 t 值排序
names(gene_list) <- rownames(deg)  # 基因名称作为名字
gene_list <- sort(gene_list, decreasing=TRUE)


# 1. hallmark
gmt_file <- "data/input/h.all.v7.5.1.symbols.gmt"  # GMT 文件路径
gene_set <- read.gmt(gmt_file)
gsea_result_1 <- GSEA(geneList=gene_list, TERM2GENE=gene_set, pvalueCutoff=0.05)

ridge <- ridgeplot(gsea_result_1,
                   showCategory = 15,         # 显示前 30 个通路
                   fill = "p.adjust",        # 填充色根据调整后的 p 值
                   core_enrichment = TRUE,    # 只使用 core enriched 基因
                   label_format = 30,         # 标签格式化
                   orderBy = "NES",           # 按 NES 排序
                   decreasing = TRUE) +       # 降序排列
  theme(axis.text.y = element_text(size = 8))+  # 调整 y 轴文本大小
  labs(x = "NES", y = "Pathways")  # 添加 x 轴和 y 轴名称
ridge
ggsave("output/revise group/B.GSEA_hallmark_ridgeplot.pdf", plot = ridge, width = 6, height = 6)

# 2. GO
gmt_file <- "data/input/c5.go.v7.5.1.symbols.gmt"  # GMT 文件路径
gene_set <- read.gmt(gmt_file)
# 进行 GSEA
gsea_result_2 <- GSEA(geneList=gene_list, TERM2GENE=gene_set, pvalueCutoff=0.05)
# ridgeplot
ridge <- ridgeplot(gsea_result_2,
                   showCategory = 15,         # 显示前 30 个通路
                   fill = "p.adjust",        # 填充色根据调整后的 p 值
                   core_enrichment = TRUE,    # 只使用 core enriched 基因
                   label_format = 30,         # 标签格式化
                   orderBy = "NES",           # 按 NES 排序
                   decreasing = TRUE) +       # 降序排列
  theme(axis.text.y = element_text(size = 8))+  # 调整 y 轴文本大小
  labs(x = "NES", y = "Pathways")  # 添加 x 轴和 y 轴名称
ridge
ggsave("output/revise group/GSEA_GO_ridgeplot.pdf", plot = ridge, width = 6, height = 6)


# 3.KEGG
gmt_file <- "data/input/c2.cp.kegg.v7.4.symbols.gmt"  # GMT 文件路径
gene_set <- read.gmt(gmt_file)
gsea_result_3 <- GSEA(geneList=gene_list, TERM2GENE=gene_set, pvalueCutoff=0.05)

ridge <- ridgeplot(gsea_result_3,
                   showCategory = 15,         # 显示前 30 个通路
                   fill = "p.adjust",        # 填充色根据调整后的 p 值
                   core_enrichment = TRUE,    # 只使用 core enriched 基因
                   label_format = 30,         # 标签格式化
                   orderBy = "NES",           # 按 NES 排序
                   decreasing = TRUE) +       # 降序排列
  theme(axis.text.y = element_text(size = 8))+  # 调整 y 轴文本大小
  labs(x = "NES", y = "Pathways")  # 添加 x 轴和 y 轴名称
ridge
ggsave("output/revise group/GSEA_KEGG_ridgeplot.pdf", plot = ridge, width = 6, height = 6)


# 4.reactome
gmt_file <- "data/input/c2.cp.reactome.v2024.1.Hs.symbols.gmt"  # GMT 文件路径
gene_set <- read.gmt(gmt_file)
gsea_result_4 <- GSEA(geneList=gene_list, TERM2GENE=gene_set, pvalueCutoff=0.05)
ridge <- ridgeplot(gsea_result_4,
                   showCategory = 15,         # 显示前 30 个通路
                   fill = "p.adjust",        # 填充色根据调整后的 p 值
                   core_enrichment = TRUE,    # 只使用 core enriched 基因
                   label_format = 30,         # 标签格式化
                   orderBy = "NES",           # 按 NES 排序
                   decreasing = TRUE) +       # 降序排列
  theme(axis.text.y = element_text(size = 8))+  # 调整 y 轴文本大小
  labs(x = "NES", y = "Pathways")  # 添加 x 轴和 y 轴名称
ridge
ggsave("output/revise group/B.GSEA_reactome_ridgeplot.pdf", plot = ridge, width = 6, height = 6)





### 单条通路
## GSEA图可视化优化
library(GseaVis)
gseaNb(object = gsea_result_1,
       geneSetID= "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
       addPval = T,
       htHeight = 0.6,
       subPlot = 2)

# 1.HALLMARK
pathways <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)
# 循环生成并保存图像
for (pw in pathways) {
  
  file_name <- paste0("output/revise group/C.", pw, ".pdf")
  
  pdf(file_name, width = 6, height = 4)
  print(
    gseaNb(
      object = gsea_result_1,
      geneSetID = pw,
      addPval = TRUE,
      htHeight = 0.6,
      subPlot = 2
    )
  )
  dev.off()
  
  message("Saved: ", file_name)
}

# GO
head(gsea_result_2@result[["ID"]], 15)
pdf("output/revise group/D.GO_GSEA.pdf", width = 10, height = 4)
gseaNb(object = gsea_result_2,
       geneSetID= c("GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX",
                    "GOBP_AEROBIC_RESPIRATION",
                    "GOCC_MITOCHONDRIAL_MATRIX"
       ),
       addPval = T,
       htHeight = 0.6,
       subPlot = 2)
dev.off()

# KEGG
head(gsea_result_3@result[["ID"]], 15)
pdf("output/revise group/D.KEGG_GSEA.pdf", width = 10, height = 4)
gseaNb(object = gsea_result_3,
       geneSetID= c("KEGG_OXIDATIVE_PHOSPHORYLATION",
                    "KEGG_CITRATE_CYCLE_TCA_CYCLE",
                    "KEGG_CELL_CYCLE"
       ),
       addPval = T,
       htHeight = 0.6,
       subPlot = 2)
dev.off()

# reactome
head(gsea_result_4@result[["ID"]], 15)
pdf("output/revise group/D.REACTOME_GSEA.pdf", width = 10, height = 4)
gseaNb(object = gsea_result_4,
       geneSetID= c("REACTOME_MITOCHONDRIAL_TRANSLATION",
                    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT"
       ),
       addPval = T,
       htHeight = 0.6,
       subPlot = 2)
dev.off()



# 加载必要的包
library(cowplot)
# 生成三个图形对象
p1 <- gseaNb(object = gsea_result_2,
             geneSetID = c("GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX",
                           "GOBP_AEROBIC_RESPIRATION",
                           "GOCC_MITOCHONDRIAL_MATRIX"),
             addPval = T,
             htHeight = 0.6,
             subPlot = 2)

p2 <- gseaNb(object = gsea_result_3,
             geneSetID = c("KEGG_OXIDATIVE_PHOSPHORYLATION",
                           "KEGG_CITRATE_CYCLE_TCA_CYCLE",
                           "KEGG_CELL_CYCLE"),
             addPval = T,
             htHeight = 0.6,
             subPlot = 2)

p3 <- gseaNb(object = gsea_result_4,
             geneSetID = c("REACTOME_MITOCHONDRIAL_TRANSLATION",
                           "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                           "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT"),
             addPval = T,
             htHeight = 0.6,
             subPlot = 2)

# 组合图形并输出
pdf("output/revise group/D.Combined_GSEA_2.pdf", width = 10, height = 12)
plot_grid(p1, p2, p3, 
          ncol = 1, 
          labels = c("A", "B", "C"),  # 添加标签
          label_size = 20,
          align = 'v',  # 垂直对齐
          rel_heights = c(1, 1, 1))   # 三个图等高
dev.off()
