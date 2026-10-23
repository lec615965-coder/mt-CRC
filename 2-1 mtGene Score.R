### mtGene Score
rm(list = ls())
### 加载R包
library(openxlsx)
library(Seurat)
library(AUCell)
library(GSVA)
### 加载数据
## 读入mtGene
mtGene <- read.xlsx("data/mtGene/mitochondria gene_geneID.xlsx") 
mtGene <- mtGene$Gene.Name
## 读入seurat
load("data/GSE132465_annotation.rda")


## 1.AddModuleScore
gene_set <- mtGene
# 添加模块评分
GSE132465_2 <- AddModuleScore(object = GSE132465_2, 
                             features = list(gene_set), 
                             name = "mtGene_ModuleScore")
head(GSE132465_2@meta.data)

## 2.AUC评分
gene_set <- mtGene
geneSet <- list(Module = gene_set)
#表达矩阵
exprMatrix <- GSE132465_2@assays$RNA@data
#计算AUC评分
cell_rankings <- AUCell_buildRankings(exprMatrix) #构建基因排序矩阵
cells_AUC <- AUCell_calcAUC(geneSet, cell_rankings) #计算AUC评分
auc_scores <- as.data.frame(t(getAUC(cells_AUC))) #提取AUC分数 
GSE132465_2 <- AddMetaData(GSE132465_2, metadata = auc_scores[, "Module"], col.name = "mtGene_AUC") #写入

head(GSE132465_2@meta.data)


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

head(GSE132465_2@meta.data)
colnames(GSE132465_2@meta.data)[(ncol(GSE132465_2@meta.data) - 2):ncol(GSE132465_2@meta.data)] <- 
  c("mtGene_ModuleScore", "mtGene_AUC" ,"mtGene_GSVA")
head(GSE132465_2@meta.data)


# 保存seurat
save(GSE132465_2, file = "data/GSE132465_final.rda")
load("data/GSE132465_final.rda")
