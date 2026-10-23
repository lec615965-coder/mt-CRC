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

source('data/input/lzq_sc_fun.R')

# 1.创建Seurat对象
# 读入数据 Read from GSE132465 
# 表达矩阵
GSE132465 <- data.table::fread("data/CRC---GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt",data.table = F)
GSE132465[1:3,1:3]
rownames(GSE132465) <- GSE132465$Index 
GSE132465 <- GSE132465[,-1]
#原文注释
anno <- data.table::fread("data/CRC---GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt",data.table = F)
rownames(anno) <- anno$Index
anno <- anno[,-1]
# 创建seurat
GSE132465 <- CreateSeuratObject(GSE132465, min.cells = 3, project = "GSE132465") # 读取数据
GSE132465 <- AddMetaData(object = GSE132465, metadata = anno)
table(GSE132465$orig.ident)
table(GSE132465$Cell_type)

# Remove doublet
obj.ls <- SplitObject(GSE132465,split.by = 'orig.ident') # 分成list，分别做doublet
rm(GSE132465)
obj.ls_rmd <- obj.ls %>% map(lzq_rm_doublet) 
GSE132465 <- Reduce(function(x, y) merge(x, y, merge.data = T), obj.ls_rmd)
GSE132465$Cohort <- 'GSE132465'
rm(obj.ls);rm(obj.ls_rmd)

GSE132465@meta.data %<>% mutate(log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA))
save(GSE132465,file = "data/GSE132465_rmd.rda") #运行时间过长，保存运行后结果


# 2. Quality control -------------------------------------------------
rm(list = ls())
load("data/GSE132465_rmd.rda")
source('input/lzq_sc_fun.R')
# 添加线粒体、核糖体、红细胞基因的比例
GSE132465 <- PercentageFeatureSet(GSE132465, pattern = "^MT-", col.name = 'percent.mt')
GSE132465 <- PercentageFeatureSet(GSE132465, pattern = "^RP[SL]", col.name = 'percent.ribo')
GSE132465 <- PercentageFeatureSet(GSE132465, pattern = "^HB[^(P)]", col.name = 'percent.hb')

lzq_sc_QC_pre(GSE132465,dir = 'output/1.seurat/QC/preQC/') # 文件夹ID 一定要加“/”
GSE132465 #质控前 25655 features across 62844 samples within 1 assay 

# Remove cells with low quality or potential doublets
GSE132465<- subset(GSE132465, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & # Gene numbers
                     nCount_RNA >= 1000 & nCount_RNA <= 40000 & # UMI numbers
                     log10GenesPerUMI >= 0.8 & # Novelty
                     percent.mt <= 20) # 线粒体比例

lzq_sc_QC_post(GSE132465,dir = 'output/1.seurat/QC/postQC/') # 文件夹ID 一定要加“/”
GSE132465 #质控后 25655 features across 56375 samples within 1 assay  

save(GSE132465, file = "data/GSE132465_QC.rda")


# 3. Batch removing -------------------------------------------------------
rm(list = ls())
load("data/GSE132465_QC.rda")
source('input/lzq_sc_fun.R')
options(future.globals.maxSize = 10 * 1024^3) # 设置为50 GiB

GSE132465 <- SCTransform(GSE132465, vars.to.regress = "percent.mt")
GSE132465 <- RunPCA(GSE132465,features = VariableFeatures(object = GSE132465)) 

ElbowPlot(GSE132465, ndims = 50, reduction = "pca")
ggsave("output/1.seurat/QC/elbowplot.tiff", width=6, height=4,dpi = 500) # 取dims =26

GSE132465 <- RunUMAP(GSE132465,reduction = "pca", dims = 1:26, reduction.name = "umap_naive")
GSE132465 <- RunTSNE(GSE132465,reduction = "pca", dims = 1:26, reduction.name = "tsne_naive")

GSE132465_2 <- RunHarmony(GSE132465, group.by.vars = 'orig.ident')

save(GSE132465_2, file="data/GSE132465_harmony.rda")


# 4. Find potential resolution -------------------------------------------- 
# 5. Find Neighbors and Clustering -------------------------------------------- 
# 6. UMAP 或 t-SNE -------------------------------------------- 
rm(list = ls())
load("data/GSE132465_harmony.rda")
GSE132465_2 <- FindNeighbors(GSE132465_2, dims = 1:26, verbose = F) #使用前面计算出来的dims
GSE132465_2 <- FindClusters(GSE132465_2, resolution = 1)

GSE132465_2 <- RunUMAP(GSE132465_2, reduction = "harmony", dims = 1:26, reduction.name = "umap")
GSE132465_2 <- RunTSNE(GSE132465_2, reduction = "harmony", dims = 1:26,reduction.name = "tsne")
  
plot <- DimPlot(GSE132465_2, reduction = "umap", label = F, pt.size = 0.5, group.by = "seurat_clusters")+
        DimPlot(GSE132465_2, reduction = "umap", label = F, pt.size = 0.5, group.by = "Cell_type")
plot  
ggsave(filename = "output/1.seurat/UMAP.tiff", plot = plot, width = 10, height = 5, dpi = 800)
ggsave(plot, filename = "output/1.seurat/UMAP.pdf",width = 24,height = 10,limitsize = FALSE)


save(GSE132465_2, file="data/GSE132465_findclusters.rda")


# 7.批次矫正评价
# 使用LISI指数评价
# Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat. Methods (2019)
rm(list = ls())
load("data/GSE132465_harmony.rda")
lisi.res <- list(harmony = compute_lisi(X = Embeddings(GSE132465_2,reduction = 'harmony'),
                                        meta_data = GSE132465_2@meta.data,
                                        label_colnames = c('orig.ident','Cohort'))%>%mutate(type='harmony'),
                 PCA = compute_lisi(X = Embeddings(GSE132465_2,reduction = 'pca'),
                                    meta_data = GSE132465_2@meta.data,
                                    label_colnames = c('orig.ident','Cohort'))%>%mutate(type='Raw')
)%>%Reduce(rbind,.)
save(lisi.res,file = "data/lisi.res.rda")

plot1 <- ggplot(lisi.res,aes(type,orig.ident))+
  geom_boxplot()+
  scale_y_continuous(breaks = c(0:10))
plot1
ggsave("output/1.seurat/QC/batch remove_plot1.pdf", plot1, width = 8, height = 6)

plot2 <- ggplot(lisi.res,aes(orig.ident,fill=type))+
  geom_density(alpha=0.2)+
  facet_wrap(~type)+
  scale_x_continuous(breaks = c(0:10))
plot2
ggsave("output/1.seurat/QC/batch remove_plot2.pdf", plot2, width = 8, height = 6)

plot3 <- DimPlot(GSE132465_2, reduction = "umap_naive",group.by = "orig.ident",raster=FALSE)+
  DimPlot(GSE132465_2, reduction = "umap",group.by = "orig.ident",raster=FALSE)
plot3
ggsave("output/1.seurat/QC/batch remove_UMAP.pdf", plot3, width = 18, height = 6)

plot4 <- DimPlot(GSE132465_2, reduction = "umap_naive",group.by = "Cell_type",raster=FALSE)+
  DimPlot(GSE132465_2, reduction = "umap",group.by = "Cell_type",raster=FALSE)
plot4
ggsave("output/1.seurat/QC/batch remove_UMAP2.pdf", plot4, width = 18, height = 6)



