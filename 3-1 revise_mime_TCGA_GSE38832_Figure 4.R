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
  select(sample_id, patient_id, OS_time, OS,
         age_at_diagnosis_years, sex, pathologic_stage,
         pathologic_T, pathologic_N, pathologic_M, MSI_status)
# 如果一个患者有多个原发肿瘤样本，保留一个
clin <- clin %>%
  mutate(vial_priority = ifelse(substr(sample_id, 16, 16) == "A", 0, 1)) %>%
  arrange(patient_id, vial_priority, sample_id) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  select(-vial_priority)

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


###### GSE38832
GSE38832_expr_mime <- read.csv("data/cohort/GSE38832_expr.csv", row.names = 1)
GSE38832_clin_mime <- read.csv("data/cohort/GSE38832_clin.csv", row.names = 1)

# mime data
expr_transposed <- as.data.frame(t(GSE38832_expr_mime))
expr_transposed$ID <- rownames(expr_transposed) 

GSE38832_data <- merge(
  GSE38832_clin_mime[, c("ID", "OS.time", "OS")],  # 选择需要的临床列
  expr_transposed,
  by = "ID")

list_train_vali_Data <- list()
list_train_vali_Data[["TCGA_COAD"]] <- TCGA_COAD_data
list_train_vali_Data[["GSE38832"]] <- GSE38832_data

save(list_train_vali_Data, file = "data/revise mime/TCGA_38832_list_train_vali_Data.rda")




### 基因集
MitoCarta <-read.xlsx("data/mtGene/Human.MitoCarta3.0.xlsx", sheet = 2) 
MitoCarta <- MitoCarta$Symbol

###### 1.构建预后模型
load("data/revise mime/TCGA_38832_list_train_vali_Data.rda")
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
save(res, file = "data/revise mime/TCGA_38832_mime_res.rda")
load("data/revise mime/TCGA_38832_mime_res.rda")

gene <- res$Sig.genes
writeLines(gene, "data/revise mime/TCGA_38832_sig gene.txt")


rm(list = ls())
### 1.1 绘制每个模型的C-index
load("data/revise mime/TCGA_38832_list_train_vali_Data.rda")
load("data/revise mime/TCGA_38832_mime_res.rda")

cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],
               order =names(list_train_vali_Data),width = 0.35)

# 保存为pdf
pdf("output/revise_mime/1.MitoCarta-TCGA_GSE38832_cindex_plot.pdf", width = 8, height = 12)  # 可调整宽高

cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.35)
dev.off()


rm(list = ls())
load("data/revise mime/TCGA_38832_list_train_vali_Data.rda")
load("data/revise mime/TCGA_38832_mime_res.rda")


# 在不同数据集之间绘制特定模型的C-index
### C-index高且稳定的模型性能好

# 重新生成KM曲线
survplot <- list()

for (i in seq_along(list_train_vali_Data)) {
  
  dataset_name <- names(list_train_vali_Data)[i]
  
  survplot[[i]] <- rs_sur(
    res,
    model_name = "RSF + SuperPC",
    dataset = dataset_name,
    median.line = "hv",
    cutoff = 0.5,
    conf.int = TRUE,
    xlab = "Year",
    pval.coord = c(1, 0.3)
  )
  
}

# 拼图
combined_plot <- aplot::plot_list(
  gglist = survplot,
  ncol = 2
)

# 保存PDF
ggsave(
  filename = "output/revise_mime/2.Survival_Plots_combined.pdf",
  plot = combined_plot,
  width = 10,
  height = 5
)



### 1.3 特定模型的单变量 Cox 回归元分析
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = "RSF + SuperPC",
                                   type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)

meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))

# 保存为pdf
pdf("output/revise_mime/4.cox_plot.pdf", width = 8, height = 4)  # 可调整宽高

meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
dev.off()




###### 3.核心特征选择
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$TCGA_COAD,
                                              candidate_genes = genelist,
                                              mode = "all",nodesize =5,seed = 5201314 )
write.csv(res.feature.all, file = "data/revise_mime/res.feature.all .csv")

# 通过不同方法筛选的基因的Upset plot
core_feature_select(res.feature.all)

pdf("output/revise_mime/upset plot.pdf", width = 8, height = 8)  # 可调整宽高
core_feature_select(res.feature.all)
dev.off()

# 绘制不同方法过滤的基因序列
core_feature_rank(res.feature.all, top=20)

pdf("output/revise_mime/core_feature.pdf", width = 8, height = 8)  # 可调整宽高
core_feature_rank(res.feature.all, top=20)
dev.off()


###### RiskScore 提取
TCGA_RS <- res$riskscore[["RSF + SuperPC"]][["TCGA_COAD"]]
TCGA_RS$Group <- ifelse(TCGA_RS$RS >= median(TCGA_RS$RS), "High", "Low")

GSE_RS <- res$riskscore[["RSF + SuperPC"]][["GSE38832"]]
GSE_RS$Group <- ifelse(GSE_RS$RS >= median(GSE_RS$RS), "High", "Low")

write.csv(TCGA_RS, "data/revise mime/TCGA_RS.csv")
write.csv(GSE_RS, "data/revise mime/GSE_RS.csv")

