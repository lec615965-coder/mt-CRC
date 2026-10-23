
library(data.table)
library(dplyr)
library(survival)
library(coxphf)
library(ggplot2)
library(patchwork)

# =====================================================
# 1. 数据读取与样本匹配
# =====================================================
expr_file <- "data/cohort/TCGA_COAD_expr.csv"
clin_file <- "data/cohort/TCGA_COAD_clin.csv"
rs_file <- "data/revise mime/TCGA_RS.csv"
outdir <- "output/revise_mime/Independent_Cox"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 读取表达矩阵表头，匹配共同样本
header <- readLines(expr_file, n = 1, warn = FALSE)
expr_ids <- names(read.csv(text = header, check.names = FALSE))[-1]

clin <- fread(clin_file)
rs <- fread(rs_file)[, .(sample_id = as.character(ID), RiskScore = as.numeric(RS))]

stopifnot(!anyDuplicated(rs$sample_id), !anyDuplicated(clin$sample_id))
common_ids <- Reduce(intersect, list(rs$sample_id, clin$sample_id, expr_ids))
if (!length(common_ids)) stop("未匹配到共同样本，请检查样本ID格式")

rs <- rs[sample_id %in% common_ids]
clin <- clin[match(rs$sample_id, clin$sample_id)]
dat <- copy(clin)
dat[, RiskScore := rs$RiskScore]

stopifnot(identical(dat$sample_id, rs$sample_id))
cat("共同样本数:", nrow(dat), "\n")

# =====================================================
# 2. 临床变量整理
# =====================================================
# 生存信息及RiskScore标准化
dat <- dat[is.finite(OS_time_days) & OS_time_days > 0 &
             OS_event %in% c(0, 1) & is.finite(RiskScore)]

dat$OS.time <- as.numeric(dat$OS_time_days) / 365.25
dat$OS <- as.numeric(dat$OS_event)
dat$RiskScore_z <- as.numeric(scale(dat$RiskScore))

# 年龄与性别
dat$Age <- as.numeric(dat$age_at_diagnosis_years)
dat$Gender <- factor(tolower(trimws(as.character(dat$sex))), levels = c("female", "male"))

# 病理分期
stage <- toupper(trimws(as.character(dat$pathologic_stage)))
dat$Stage <- case_when(
  grepl("^STAGE IV", stage) ~ "IV",
  grepl("^STAGE III", stage) ~ "III",
  grepl("^STAGE II", stage) ~ "II",
  grepl("^STAGE I", stage) ~ "I",
  TRUE ~ NA_character_
)
dat$Stage <- factor(dat$Stage, levels = c("I", "II", "III", "IV"))

# 合并病理分期
dat$Stage_group <- case_when(
  dat$Stage %in% c("I", "II") ~ "I-II",
  dat$Stage %in% c("III", "IV") ~ "III-IV",
  TRUE ~ NA_character_
)
dat$Stage_group <- factor(dat$Stage_group, levels = c("I-II", "III-IV"))

# MSI状态
if (!"MSI_status" %in% names(dat)) stop("临床数据缺少MSI_status列")
dat$MSI <- trimws(as.character(dat$MSI_status))
dat$MSI[dat$MSI %in% c("", "NA", "N/A", "Unknown", "[Not Available]")] <- NA_character_
dat$MSI <- factor(dat$MSI, levels = c("MSI-negative", "MSI-positive"))

print(table(dat$MSI_status, useNA = "ifany"))

# =====================================================
# 3. 构建分析队列
# =====================================================
vars_all <- c("OS.time", "OS", "RiskScore_z", "Age", "Gender", "Stage")
vars_msi <- c("OS.time", "OS", "RiskScore_z", "Age", "Gender", "Stage_group", "MSI")

dat_all <- droplevels(as.data.frame(dat[complete.cases(dat[, ..vars_all])]))
dat_msi <- droplevels(as.data.frame(dat[complete.cases(dat[, ..vars_msi])]))

cohort_info <- data.frame(
  Cohort = c("All patients", "MSI-complete"),
  N = c(nrow(dat_all), nrow(dat_msi)),
  Events = c(sum(dat_all$OS), sum(dat_msi$OS))
)

print(cohort_info)
print(table(dat_msi$MSI, dat_msi$OS))
print(table(dat_msi$Stage_group, dat_msi$OS))
write.csv(cohort_info, file.path(outdir, "Cohort_summary.csv"), row.names = FALSE)

# =====================================================
# 4. Cox回归结果提取
# =====================================================
extract_cox <- function(fit) {
  s <- summary(fit)
  data.frame(
    Variable = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"],
    Lower95 = s$conf.int[, "lower .95"],
    Upper95 = s$conf.int[, "upper .95"],
    Pvalue = s$coefficients[, "Pr(>|z|)"],
    Method = "Cox", row.names = NULL
  )
}

extract_firth <- function(fit) {
  data.frame(
    Variable = names(fit$coefficients),
    HR = exp(fit$coefficients),
    Lower95 = as.numeric(fit$ci.lower),
    Upper95 = as.numeric(fit$ci.upper),
    Pvalue = as.numeric(fit$prob),
    Method = "Firth Cox", row.names = NULL
  )
}
# 注意：coxphf的ci.lower和ci.upper已处于HR尺度，无需再次exp()

# =====================================================
# 5. 单因素回归函数
# =====================================================
run_univariate <- function(data, variables, method = "cox") {
  bind_rows(lapply(variables, function(v) {
    fml <- as.formula(paste0("Surv(OS.time, OS) ~ ", v))
    if (method == "firth") {
      extract_firth(coxphf(fml, data = data, maxit = 1000, maxstep = 0.1))
    } else {
      extract_cox(coxph(fml, data = data))
    }
  }))
}

# =====================================================
# 6. 全部患者：普通Cox
# =====================================================
uni_all <- run_univariate(dat_all, c("RiskScore_z", "Age", "Gender", "Stage"))
multi_all_fit <- coxph(Surv(OS.time, OS) ~ RiskScore_z + Age + Gender + Stage, data = dat_all)
multi_all <- extract_cox(multi_all_fit)

# 比例风险假设检验
ph_test <- cox.zph(multi_all_fit)
print(ph_test)
write.csv(as.data.frame(ph_test$table), file.path(outdir, "PH_assumption.csv"))

# =====================================================
# 7. MSI完整患者：Firth Cox
# =====================================================
uni_msi <- run_univariate(dat_msi, c("RiskScore_z", "Age", "Gender", "Stage_group", "MSI"), method = "firth")
multi_msi_fit <- coxphf(
  Surv(OS.time, OS) ~ RiskScore_z + Age + Gender + Stage_group + MSI,
  data = dat_msi, maxit = 1000, maxstep = 0.1
)
multi_msi <- extract_firth(multi_msi_fit)

lapply(list(uni_all, multi_all, uni_msi, multi_msi), print)

# =====================================================
# 8. 变量名称格式化
# =====================================================
format_names <- function(df) {
  df$Variable <- recode(df$Variable,
                        "RiskScore_z" = "Risk score (per SD)",
                        "Age" = "Age (per year)",
                        "Gendermale" = "Gender (Male vs Female)",
                        "StageII" = "Stage II vs I",
                        "StageIII" = "Stage III vs I",
                        "StageIV" = "Stage IV vs I",
                        "Stage_groupIII-IV" = "Stage III-IV vs I-II",
                        "MSIMSI-positive" = "MSI positive vs negative"
  )
  df
}

# =====================================================
# 9. 森林图绘制函数
# =====================================================
plot_cox <- function(df, title) {
  df <- format_names(df)
  df$Significance <- ifelse(!is.na(df$Pvalue) & df$Pvalue < 0.05, "P < 0.05", "P >= 0.05")
  df$P_label <- ifelse(is.na(df$Pvalue), "NA",
                       ifelse(df$Pvalue < 0.001, "<0.001", sprintf("%.3f", df$Pvalue)))
  
  # 判断HR与置信区间是否可绘制
  df$Estimable <- with(df, is.finite(HR) & is.finite(Lower95) &
                         is.finite(Upper95) & HR > 0 & Lower95 > 0 & Upper95 > 0)
  df$Label <- ifelse(df$Estimable,
                     sprintf("%.2f (%.2f-%.2f)  P=%s", df$HR, df$Lower95, df$Upper95, df$P_label),
                     "Not estimable")
  df$Variable <- factor(df$Variable, levels = rev(unique(df$Variable)))
  df_plot <- df[df$Estimable, ]
  
  if (!nrow(df_plot)) {
    return(ggplot() + theme_void() + labs(title = title, subtitle = "No estimable HR"))
  }
  
  # 坐标范围与标签位置
  xmin <- min(df_plot$Lower95)
  xmax <- max(df_plot$Upper95)
  left_limit <- max(min(xmin / 1.5, 0.8), 0.001)
  label_x <- xmax * 1.5
  right_limit <- label_x * 4.5
  
  ggplot(df, aes(y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbar(
      data = df_plot, aes(x = HR, xmin = Lower95, xmax = Upper95, color = Significance),
      orientation = "y", width = 0.2, linewidth = 0.6
    ) +
    geom_point(data = df_plot, aes(x = HR, color = Significance), size = 2.6) +
    geom_text(aes(x = label_x, label = Label), hjust = 0, size = 3, color = "black") +
    scale_color_manual(values = c("P < 0.05" = "#D73027", "P >= 0.05" = "#2874A6")) +
    scale_x_log10() +
    coord_cartesian(xlim = c(left_limit, right_limit), clip = "off") +
    labs(title = title, x = "Hazard Ratio (95% CI)", y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      legend.position = "none",
      plot.margin = margin(10, 15, 10, 10)
    )
}

# =====================================================
# 10. 绘制并合并森林图
# =====================================================
p1 <- plot_cox(uni_all, "A. All patients - Univariate Cox")
p2 <- plot_cox(multi_all, "B. All patients - Multivariate Cox")
p3 <- plot_cox(uni_msi, "C. MSI-complete - Univariate Firth Cox")
p4 <- plot_cox(multi_msi, "D. MSI-complete - Multivariate Firth Cox")

combined_plot <- (p1 | p2) / (p3 | p4)
print(combined_plot)

ggsave(file.path(outdir, "Combined_Cox_Forest.pdf"),
       plot = combined_plot, width = 12, height = 8, limitsize = FALSE)
