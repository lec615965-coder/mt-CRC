library(Seurat)
library(ggplot2)
library(patchwork)

load("data/GSE132465_final.rda")

signature_genes <- c(
  "LARS2", "ALDH2", "COQ10A", "MRPS7", "ACADL", "CMC1",
  "CPT1C",  "PDHB", "PTRH1", "COX4I2", "GSR",
  "SPR", "MRPS5", "FASTK", "FAHD2A", "SLC25A45")

p_dot <- DotPlot(
  GSE132465_2,
  features = signature_genes,
  group.by = group_var,
  scale = FALSE,
  dot.scale = 7
) +
  scale_color_gradientn(colours = expression_colours, name = "Average expression") +
  scale_size(range = c(0.2, 7), name = "Percent expressed") +
  labs(x = NULL, y = NULL) +
  theme_revision +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2)
  )+coord_flip()

ggsave(
  "signature_genes_DotPlot_combined.pdf",
  plot = p_dot,
  width = 8,
  height = 5,
  units = "in",
  limitsize = FALSE,
  useDingbats = FALSE
)
