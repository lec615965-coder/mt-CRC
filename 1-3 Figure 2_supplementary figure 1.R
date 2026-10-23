rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(patchwork)

load("data/GSE132465_final.rda")

###### A, B UMAP
my_colors <- c(
  "#1F77B4", "#17BECF", "#00BFC4", "#1E90FF", "#AEC7E8",
  "#2CA02C", "#66A61E", "#98DF8A", "#A6D854", "#BCBD22",
  "#FFD92F", "#E6AB02", "#D95F02", "#FF7F0E", "#FFBB78",
  "#F5A65B", "#D62728", "#FF9896", "#E78AC3", "#E377C2",
  "#F7B6D2", "#9467BD", "#C5B0D5", "#8DA0CB", "#7570B3",
  "#8C564B", "#C49C94", "#A6761D", "#7F7F7F", "#B3B3B3",
  "#9EDAE5", "#FC8D62", "#DBDB8D", "#E7298A", "#264653"
)

p1 <- DimPlot(
  GSE132465_2,
  reduction = "umap_naive",
  group.by = "seurat_clusters",
  cols = my_colors,
  label = TRUE,
  pt.size = 0.4,
  label.size = 4
) +
  ggtitle("Seurat clusters") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1),
    legend.key.size = unit(0.4, "cm"),        
    legend.spacing.y = unit(0.1, "cm"),       
    legend.box.spacing = unit(0.2, "cm")      
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3),          
      ncol = 2,                               
      byrow = TRUE
    )
  )
p1

ggsave(
  filename = "output/Figure 2/A.UMAP_Seurat_Clusters.pdf",
  plot = p1,
  width = 8,
  height = 6,
  units = "in"
)


p2 <- DimPlot(
  GSE132465_2,
  reduction = "umap_naive",
  group.by = "bulk_labels_3",
  cols = my_colors,
  label = TRUE,
  pt.size = 0.4,
  label.size = 4
) +
  ggtitle("Cell types") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1),
    legend.key.size = unit(0.4, "cm"),       
    legend.spacing.y = unit(0.1, "cm"),       
    legend.box.spacing = unit(0.2, "cm")     
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3),          
      ncol = 1,                               
      byrow = TRUE
    )
  )
p2

ggsave(
  filename = "output/Figure 2/B.Cell Types.pdf",
  plot = p2,
  width = 8,
  height = 6,
  units = "in"
)


###### C. Dotplot
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


p3 <- 
  DotPlot(GSE132465_2, group.by = "bulk_labels_3", features = Marker) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 13, face = "plain"), 
        legend.title = element_text(size = 13, face = "bold")) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colours = c("#1F77B4","#9EDAE5", "#FF9896","#D62728")) + 
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3))
p3


ggsave(filename = "output/Figure 2/C.celltype marker_dotplot.pdf", plot = p3,
       width = 16,height = 6,limitsize = FALSE)



###### D.FeaturePlot
p4 <- FeaturePlot(
  GSE132465_2,
  features = "mtGene_ModuleScore",       
  reduction = "umap_naive",              
  pt.size = 0.4                  
) +
  ggtitle("Mitochondria gene score") +  
  scale_color_gradientn(
    colours = c("#1F77B4","#9EDAE5", "#FF9896","#D62728"),  
    name = "Score"                                             
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.4, "cm")
  )
p4
ggsave("output/Figure 2/D.FeaturePlot_mtGeneScore.pdf", 
       plot = p4, width = 7, height = 6)


###### E.Violin plot
p5 <- VlnPlot(
  GSE132465_2,
  features = "mtGene_ModuleScore",
  group.by = "bulk_labels_3",
  pt.size = 0,         
  cols = my_colors           
) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", alpha = 1, outlier.shape = NA) +
  ggtitle("")+
  xlab("") + ylab("Mitochondria gene score") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold",angle = 30, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.position = "none" 
  )
p5

ggsave("output/Figure 2/E.ViolinPlot_Score_by_CellType.pdf", plot = p5,
       width = 9, height = 6)




###### supplementary figure 1
###### ClusterGVis 
library(ClusterGVis)
library(org.Hs.eg.db)
load("data/GSE132465_final.rda")
AS_F <- GSE132465_2
new.cluster.ids <- c("memory T cell", "malignant cell", "naive B cell", "M2", "Treg", 
                     "M1", "malignant cell", "CTL", "naive T cell", "CTL", "CTL", "TAC", 
                     "ISC", "CTL", "memory T cell", "endothelial cell", "fibroblast", "ISC", 
                     "goblet", "fibroblast", "CTL", "ISC", "fibroblast", "malignant cell", 
                     "plasma cell", "fibroblast", "naive B cell", "enterocyte", "fibroblast", 
                     "DC", "fibroblast", "memory T cell", "mast cell")
names(new.cluster.ids) <- levels(AS_F)
AS_F <- RenameIdents(AS_F, new.cluster.ids)

AS_F.markers.all <- Seurat::FindAllMarkers(AS_F,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

AS_F.markers <- AS_F.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

st.data <- prepareDataFromscRNA(object = AS_F,
                                diffData = AS_F.markers,
                                showAverage = TRUE)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

markGenes = unique(AS_F.markers$gene)[sample(1:length(unique(AS_F.markers$gene)),40,
                                             replace = F)]

st.data <- prepareDataFromscRNA(object = AS_F,
                                diffData = AS_F.markers,
                                showAverage = TRUE,
                                keep.uniqGene = FALSE, 
                                sep = "_")


pdf('output/Figure 2/heatmap_ClusterGvis_annotation.pdf',height = 19,width = 16,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 85,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:17),
           go.col = rep(jjAnno::useMyCol("stallion",n = 17),each = 5),
           add.bar = T)
dev.off()


