rm(list = ls())
library(data.table)
library(survminer)
library(survival)
library(openxlsx)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(stringr)
library(Mime1)
library(patchwork) 
library(tibble)

################## 数据准备 ##################
###### TCGA COAD
TCGACOAD_expr <- read.csv("data/cohort/TCGA_COAD_expr.csv", row.names = 1)
TCGACOAD_clin <- read.csv("data/cohort/TCGA_COAD_clin.csv")
expr <- TCGACOAD_expr
clin <- TCGACOAD_clin
# 检查原始数据
table(clin$sample_type_code, useNA = "ifany")
table(clin$Tissue, useNA = "ifany")

# 清洗临床数据
# 保留原发肿瘤（01），排除正常组织、转移灶等
clin <- clin %>%
  filter(sample_type_code == 1, Tissue == "Primary Tumor") %>%
  mutate(OS_time = as.numeric(OS_time_days),
         OS = as.numeric(OS_event)) %>%
  filter(!is.na(OS_time), OS_time > 0,
         !is.na(OS), OS %in% c(0, 1))
# 标准化Sample ID
clin$sample_id <- gsub("\\.", "-", clin$sample_id)
clin$patient_id <- substr(clin$sample_id, 1, 12)
# 保留后续生存分析需要的信息
clin <- clin %>%
  dplyr::select(sample_id, patient_id, OS_time, OS,
                age_at_diagnosis_years, sex, pathologic_stage,
                pathologic_T, pathologic_N, pathologic_M, MSI_status)
# 如果一个患者有多个原发肿瘤样本，保留一个
clin <- clin %>%
  mutate(vial_priority = ifelse(substr(sample_id, 16, 16) == "A", 0, 1)) %>%
  arrange(patient_id, vial_priority, sample_id) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  dplyr::select(-vial_priority)

# 清洗表达矩阵
# TCGA表达矩阵列名格式：TCGA.XX.XXXX.01A
colnames(expr) <- gsub("\\.", "-", colnames(expr))
# 仅保留原发肿瘤样本（01）
sample_type <- substr(colnames(expr), 14, 15)
expr <- expr[, sample_type == "01", drop = FALSE]
# 只保留有完整生存信息的样本
common_samples <- intersect(colnames(expr), clin$sample_id)
expr <- expr[, common_samples, drop = FALSE]
clin <- clin[match(common_samples, clin$sample_id), , drop = FALSE]

#最终数据整理
TCGACOAD_expr_mime <- as.data.frame(expr)
TCGACOAD_clin_mime <- as.data.frame(clin)
rownames(TCGACOAD_clin_mime) <- TCGACOAD_clin_mime$sample_id
colnames(TCGACOAD_clin_mime)[1] <- "ID"
colnames(TCGACOAD_clin_mime)[3] <- "OS.time"
TCGACOAD_clin_mime$OS.time <- TCGACOAD_clin_mime$OS.time/365

# mime data
expr_transposed <- as.data.frame(t(TCGACOAD_expr_mime))
expr_transposed$ID <- rownames(expr_transposed) 

TCGA_COAD_data <- merge(
  TCGACOAD_clin_mime[, c("ID", "OS.time", "OS")],  # 选择需要的临床列
  expr_transposed,
  by = "ID")


###### GSE39582
GSE39582_expr <- read.csv("data/cohort/GSE39582_expr.csv", row.names = 1)
GSE39582_clin <- read.csv("data/cohort/GSE39582_clin.csv", row.names = 1)

# 筛选Tumor临床样本
GSE39582_clin_mime <- GSE39582_clin[
  GSE39582_clin$Tissue == "Tumor" & !is.na(GSE39582_clin$Tissue), , drop = FALSE]

# 根据临床数据ID匹配表达矩阵
sample_ids <- if ("ID" %in% colnames(GSE39582_clin_mime)) {
  GSE39582_clin_mime$ID
} else {
  rownames(GSE39582_clin_mime)
}
# 筛选Tumor表达矩阵
GSE39582_expr_mime <- GSE39582_expr[ , sample_ids, drop = FALSE]

# mime data
expr_transposed <- as.data.frame(t(GSE39582_expr_mime))
expr_transposed$ID <- rownames(expr_transposed) 

GSE39582_data <- merge(
  GSE39582_clin_mime[, c("ID", "OS.time", "OS")],  # 选择需要的临床列
  expr_transposed,
  by = "ID")

list_train_vali_Data <- list()
list_train_vali_Data[["TCGA_COAD"]] <- TCGA_COAD_data
list_train_vali_Data[["GSE39582"]] <- GSE39582_data

save(list_train_vali_Data, file = "data/revise mime/TCGA_39582_list_train_vali_Data.rda")




### 基因集
MitoCarta <-read.xlsx("data/mtGene/Human.MitoCarta3.0.xlsx", sheet = 2) 
MitoCarta <- MitoCarta$Symbol

###### 1.构建预后模型
load("data/revise mime/TCGA_39582_list_train_vali_Data.rda")
genelist <- MitoCarta

genelist <- genelist %>%
  str_replace("MT-", "") %>%        # 去掉 "MT-"
  str_replace("CO1", "COX1") %>%    # CO1 → COX1
  str_replace("CO2", "COX2") %>%    # CO2 → COX2
  str_replace("CO3", "COX3") %>%    # CO3 → COX3
  str_replace("CYB", "CYTB")

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA_COAD,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist, #可自定义
                       mode ='all',
                       nodesize = 5,
                       seed = 5201314)

## 保存res数据
save(res, file = "data/revise mime/TCGA_39582_mime_res.rda")
load("data/revise mime/TCGA_39582_mime_res.rda")

gene <- res$Sig.genes
writeLines(gene, "data/revise mime/TCGA_39582_sig gene.txt")


rm(list = ls())
### 1.1 绘制每个模型的C-index
load("data/revise mime/TCGA_39582_list_train_vali_Data.rda")
load("data/revise mime/TCGA_39582_mime_res.rda")

cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],
               order =names(list_train_vali_Data),width = 0.35)

# 保存为pdf
pdf("output/revise_mime/1.MitoCarta-TCGA_GSE39582_cindex_plot.pdf", width = 8, height = 12)  # 可调整宽高

cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35)
dev.off()


