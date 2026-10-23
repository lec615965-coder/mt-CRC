## 1.大类注释（免疫/基质/上皮细胞）
## 2.小类注释（T细胞/B细胞/成纤维/内皮/上皮……）
## 3.亚类注释（NKT/TCM/M1/M2/plasma/杯状细胞/肠干细胞……）

#加载R包
library(Seurat)
library(ggplot2)
#加载数据
load("data/GSE132465_findclusters.rda")


### 1.查看不同cluster的top gene(辅助注释)
# 计算所有 cluster 的差异表达基因并提取top10
degs <- FindAllMarkers(
  GSE132465_2, 
  only.pos = TRUE,              # 只考虑上调基因
  logfc.threshold = 0.25,       # 设置 log fold change 的阈值
  min.pct = 0.1                 # 基因在至少 10% 的细胞中表达
)
write.csv(degs, file = "data/DEGs_per_cluster.csv", row.names = FALSE)
degs <- read.csv( "data/DEGs_per_cluster.csv")

top10 <- degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = "data/DEGs_per_cluster_top10.csv", row.names = FALSE)

plot <- DoHeatmap(GSE132465_2, features = top10$gene) + NoLegend()

ggsave("output/2.cell grouping and annotation/annotation reference/DEGs of Clusters_Heatmap.pdf", plot = plot,width = 30,height = 50,limitsize = FALSE)



### 2.关键gene在不同cluster的表达的气泡图
#大类
Marker <- c("THY1", "COL1A1", "COL1A2", "COL6A1", "COL6A2","CDH5", "VWF", "PLVAP",
            "EPCAM", "CHGA", "GUCA2A", "KRT18", "KRT8", "MUC2", 
            "PTPRC")

DotPlot(GSE132465_2, features = Marker)+coord_flip()+theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave(filename = "output/2.cell grouping and annotation/annotation reference/dotplot_大类_Cell marker genes.pdf",width = 15,height = 5,limitsize = FALSE)

#小类
Marker <- c("CCR7", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CTLA4", "CXCR5", "CXCR6", "FOXP3", "IL7R", "KLRC", "KLRC1", "KLRD1", "NKG7", "TRAC", 
            "FGFBP2", "GNLY", "KLRD1", "KLRF1", "NKG7", "PRF1", 
            "AICDA", "CD19", "CD78B", "CD79A", "CD82A", "IGHA1", "IGHG1", "MS4A1", "MZB1", "SELL", 
            "IL3RA", "LILRA4", 
            "CD14", "CD163", "CD68", "LYZ", 
            "CD1C", "CD68", "CLEC10A", "CLEC9A", "OSM", "S100A8", "S100A9", "TPSAB1", "XCR1", 
            "CPA3", "KIT", "MS4A2", "TPSAB1", "TPSB2", "IL1RL1", 
            "CDH5", "CLDN5", "PECAM1", "VMF", 
            "EPCAM", "CHGA", "GUCA2A", "KRT18", "KRT8", "MUC2", 
            "ACTA2", "COL1A1", "COL3A1", "DCN", "FAP", "THY1")
DotPlot(GSE132465_2, features = unique(Marker))+coord_flip()+theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
ggsave(filename = "output/2.cell grouping and annotation/annotation reference/dotplot_小类_Cell marker genes.pdf",width = 20,height = 10,limitsize = FALSE)


### 3. 小提琴图
# 1. Plasma cells 浆细胞
select_genes <- c("IGHA1", "IGHG1", "MZB1", "CD82A")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Plasma cells marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 2) + NoLegend()
dev.off()

# 2. Enterocytes 肠上皮细胞
select_genes <- c("GUCA2A")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Enterocytes marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 3. Goblet 杯状细胞
select_genes <- c("MUC2")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Goblet marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 4. Enteroendocrine 肠内分泌细胞
select_genes <- c("CHGA", "CHGB", "CPE", "NEUROD1", "PYY", "REG4", "TAC1", "TPH1", "NEUROG3")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Enteroendocrine marker.pdf", width = 20, height = 20)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 2) + NoLegend()
dev.off()

# 5. CD8+ T CD8+ T细胞
select_genes <- c("CD8A")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_CD8+ T marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 6. Follicular helper T 滤泡辅助性T细胞
select_genes <- c("CXCR5","CXCR3", "CCR7")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Follicular helper T marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 2) + NoLegend()
dev.off()

select_genes <- c("TBX21","IFNG", "CXCR3","GATA3","IL4","IL5","IL13","RORC","IL17A","IL17F","IL22")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Th marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 2) + NoLegend()
dev.off()

# 7. Central memory T 中央记忆T细胞
select_genes <- c("IL7R")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Central memory T marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 8. Regulatory T 调节性T细胞
select_genes <- c("FOXP3")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Regulatory T marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 9. Natural killer T 自然杀伤性T细胞
select_genes <- c("KLRD1")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Natural killer T marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 10. Follicular B 滤泡B细胞
select_genes <- c("MS4A1")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Follicular B marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 11. Plasma 浆细胞
select_genes <- c("MZB1")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Plasma marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 12. M1 macrophage M1型巨噬细胞
select_genes <- c("CD86")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_M1 macrophage marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 13. M2 macrophage M2型巨噬细胞
select_genes <- c("CD163")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_M2 macrophage marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 14. Monocyte 单核细胞
select_genes <- c("CD14")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Monocyte marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 15. Fibroblasts 成纤维细胞
select_genes <- c("COL1A1")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Fibroblasts marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 16. Endothelial 内皮细胞
select_genes <- c("VWF")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Endothelial marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()

# 17. Mast 肥大细胞
select_genes <- c("KIT")
pdf(file = "output/2.cell grouping and annotation/annotation reference/亚类_Mast marker.pdf", width = 20, height = 10)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 1) + NoLegend()
dev.off()


### 肠细胞分类（肠上皮/肠内分泌……）
select_genes <- c("MKI67", "TOP2A", "UBE2C", "CCNA2", "CENPM", "HMGB3", "MCM5", "MND1", "PCNA", "TACC3", "TK1", "TOP2A", "UBE2T",
                  "ALPI", "FABP2", "SLC26A3", "TMEM37", "CDH1", "CLDN4",
                  "MUC2", "TFF3",
                  "SPINK4", "WFDC2", "ATOH1", "CLCA1", "FFAR4", "MKI67", "MUC2", "MYC", "PBK", "PTGDR", "SPDEF", "TFF3", "TRABD2A", "UBE2C", "ZG16", "ITLN1",
                  "DCLK1", "PHLDA1",
                  "GFI1B", "POU2F3", "TRPM5",
                  "DEFA5", "DEFA6", "IL9", "PRSS2",
                  "CDK6", "CLPS", "FABP5", "HMGA2", "LCN2", "MBP", "MDK", "MUC4", "NEUROG3", "PLA2G2A", "SOX9", "THY1",
                  "LGR5", "SMOC2", "ALCAM", "ASCL2", "BMI1", "CASP1", "GRIA4", "LEFTY1", "PLA2G2A", "PTGDR", "RGMB", "OLFM4", "SOX2", "CCKBR",
                  "FABP2", "APOE", "FAM3D")
pdf(file = "output/2.cell grouping and annotation/annotation reference/Epithelial_marker.pdf", width = 20, height = 40)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 4) + NoLegend()
dev.off()

select_genes <- c("CHGA", "CHGB", "CPE", "NEUROD1", "PYY", "REG4", "TAC1", "TPH1", "NEUROG3",
                  "CLEC4G", "GHRL", "MLN", 
                  "HHEX", "PCP4", "SST", 
                  "CCK", "GC", "TDO2", 
                  "GCG", "PYY", "UCN3", 
                  "CLEC4G", "GHRL", "MLN", 
                  "HOPX", "MEP1A", "NTS", 
                  "ISL1", "NEUROG3", "GHRL")
pdf(file = "output/2.cell grouping and annotation/annotation reference/Endocrine_marker.pdf", width = 20, height = 30)
VlnPlot(GSE132465_2, features = select_genes, pt.size = 0, ncol = 4) + NoLegend()
dev.off()


rm(list = ls())
# 注释（写入不同的手动注释并保存含有注释的seurat）
# 加载数据
load("data/GSE132465_findclusters.rda")
metadata <- GSE132465_2@meta.data

#大类
bulk_labels<-readxl::read_xlsx('data/1.annotation/GSE132465_annotation.xlsx',sheet = 2)
bulk_labels<-as.data.frame(bulk_labels)
GSE132465_2@meta.data$bulk_labels_1 = "NA"
for(i in 1:nrow(bulk_labels)){
  GSE132465_2@meta.data[which(GSE132465_2@meta.data$seurat_clusters == bulk_labels$Cluster_ID[i]),'bulk_labels_1'] <- bulk_labels$bulk_labels_1[i]}
#小类
bulk_labels<-readxl::read_xlsx('data/1.annotation/GSE132465_annotation.xlsx',sheet = 3)
bulk_labels<-as.data.frame(bulk_labels)
GSE132465_2@meta.data$bulk_labels_2 = "NA"
for(i in 1:nrow(bulk_labels)){
  GSE132465_2@meta.data[which(GSE132465_2@meta.data$seurat_clusters == bulk_labels$Cluster_ID[i]),'bulk_labels_2'] <- bulk_labels$bulk_labels_2[i]}
#亚类
bulk_labels<-readxl::read_xlsx('data/1.annotation/GSE132465_annotation.xlsx',sheet = 4)
bulk_labels<-as.data.frame(bulk_labels)
GSE132465_2@meta.data$bulk_labels_3 = "NA"
for(i in 1:nrow(bulk_labels)){
  GSE132465_2@meta.data[which(GSE132465_2@meta.data$seurat_clusters == bulk_labels$Cluster_ID[i]),'bulk_labels_3'] <- bulk_labels$bulk_labels_3[i]}

metadata <- GSE132465_2@meta.data

table(GSE132465_2$bulk_labels_1)
table(GSE132465_2$bulk_labels_2)
table(GSE132465_2$bulk_labels_3)


save(GSE132465_2, file = "data/GSE132465_annotation.rda")
save(metadata, file = "data/metadata.rda")


DimPlot(GSE132465_2, reduction = "umap_naive", label = T, pt.size = 0.5, group.by = "seurat_clusters")
DimPlot(GSE132465_2, reduction = "umap_naive", label = T, pt.size = 0.5, group.by = "bulk_labels_1")
DimPlot(GSE132465_2, reduction = "umap_naive", label = T, pt.size = 0.5, group.by = "bulk_labels_2")
DimPlot(GSE132465_2, reduction = "umap_naive", label = T, pt.size = 0.5, group.by = "bulk_labels_3")

# 设置输出为 PDF 文件
pdf(file = "output/2.cell grouping and annotation/UMAP_DimPlot_bulk_labels_1_2_3.pdf", width = 18, height = 10 ) 
  DimPlot(GSE132465_2, reduction = "umap_naive", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters") +
  DimPlot(GSE132465_2, reduction = "umap_naive", label = TRUE, pt.size = 0.5, group.by = "bulk_labels_1") +
  DimPlot(GSE132465_2, reduction = "umap_naive", label = TRUE, pt.size = 0.5, group.by = "bulk_labels_2") +
  DimPlot(GSE132465_2, reduction = "umap_naive", label = TRUE, pt.size = 0.5, group.by = "bulk_labels_3")
dev.off()


