library(survival)
library(timeROC)
library(rms)
library(ggplot2)
library(patchwork)

set.seed(20261022)
bootstrap_replicates <- 1000L

expr_file <- "data/cohort/TCGA_COAD_expr.csv"
rs_file   <- "data/revise mime/TCGA_RS.csv"
out_dir   <- "output/revise_mime/TimeROC_Calibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("survival", "timeROC", "rms", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required R package(s): ", paste(missing_packages, collapse = ", "))
}
if (!all(file.exists(c(expr_file, rs_file)))) {
  stop("Missing input file(s): ",
       paste(c(expr_file, rs_file)[!file.exists(c(expr_file, rs_file))], collapse = ", "))
}

# Final group information provided by the analysis workflow.
expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
colnames(expr) <- gsub("\\.", "-", colnames(expr))
TCGA_RS <- read.csv(rs_file, check.names = FALSE, stringsAsFactors = FALSE)

required_rs_columns <- c("ID", "OS.time", "OS", "RS", "Group")
missing_rs_columns <- setdiff(required_rs_columns, colnames(TCGA_RS))
if (length(missing_rs_columns)) {
  stop("TCGA_RS.csv is missing: ", paste(missing_rs_columns, collapse = ", "))
}
if (anyDuplicated(TCGA_RS$ID)) {
  stop("TCGA_RS.csv contains duplicated sample IDs; resolve these before analysis.")
}

# Keep the exact final Group alignment used in the requested code.
RS_group <- TCGA_RS[, c("ID", "Group"), drop = FALSE]
RS_continuous <- TCGA_RS[, c("ID", "OS.time", "OS", "RS"), drop = FALSE]
common_samples <- intersect(colnames(expr), RS_group$ID)
if (!length(common_samples)) stop("No overlapping samples between expression matrix and TCGA_RS.csv.")

expr <- expr[, common_samples, drop = FALSE]
RS_group <- RS_group[match(common_samples, RS_group$ID), , drop = FALSE]
RS_continuous <- RS_continuous[match(common_samples, RS_continuous$ID), , drop = FALSE]
stopifnot(identical(colnames(expr), RS_group$ID), identical(RS_group$ID, RS_continuous$ID))

# OS.time in the locked risk-score table is in years.  timeROC and rms use days
# here so 1, 3, and 5 years are represented explicitly and consistently.
dat <- data.frame(
  ID = RS_group$ID,
  Group = factor(RS_group$Group, levels = c("Low", "High")),
  RS = suppressWarnings(as.numeric(RS_continuous$RS)),
  OS_time_years = suppressWarnings(as.numeric(RS_continuous$OS.time)),
  OS_event = suppressWarnings(as.integer(RS_continuous$OS)),
  stringsAsFactors = FALSE
)
dat$OS_time_days <- dat$OS_time_years * 365.25

flow <- data.frame(
  step = c("Expression samples", "Samples in final Group table", "Expression/Group overlap",
           "Valid OS and continuous RS", "Final analysis cohort"),
  n = c(ncol(expr), nrow(TCGA_RS), length(common_samples), NA_integer_, NA_integer_)
)
valid <- is.finite(dat$OS_time_days) & dat$OS_time_days > 0 &
  dat$OS_event %in% c(0L, 1L) & is.finite(dat$RS) & !is.na(dat$Group)
flow$n[4] <- sum(valid)
dat <- dat[valid, , drop = FALSE]
flow$n[5] <- nrow(dat)
if (nrow(dat) < 20L || sum(dat$OS_event) < 10L) {
  stop("Too few analysable observations/events after matching final Group, OS, and RS.")
}

write.csv(flow, file.path(out_dir, "01_cohort_flow.csv"), row.names = FALSE)
write.csv(dat, file.path(out_dir, "02_analysis_cohort.csv"), row.names = FALSE)
write.csv(as.data.frame(table(Group = dat$Group)),
          file.path(out_dir, "03_final_group_counts.csv"), row.names = FALSE)

# Prespecified reporting horizons.  A horizon beyond observed follow-up is not
# estimable and is excluded with a recorded warning rather than extrapolated.
horizons_years_requested <- c(1, 3, 5)
horizons_days_requested <- horizons_years_requested * 365.25
horizon_available <- horizons_days_requested < max(dat$OS_time_days)
if (!all(horizon_available)) {
  warning("Excluding horizon(s) beyond maximum observed follow-up: ",
          paste(horizons_years_requested[!horizon_available], collapse = ", "), " years")
}
horizons_years <- horizons_years_requested[horizon_available]
horizons_days <- horizons_days_requested[horizon_available]
if (!length(horizons_days)) stop("No requested horizon is within observed follow-up.")

horizon_summary <- data.frame(
  horizon_years = horizons_years,
  horizon_days = horizons_days,
  n_at_risk = vapply(horizons_days, function(t) sum(dat$OS_time_days >= t), integer(1)),
  events_by_horizon = vapply(horizons_days, function(t) sum(dat$OS_time_days <= t & dat$OS_event == 1L), integer(1))
)
write.csv(horizon_summary, file.path(out_dir, "04_horizon_support.csv"), row.names = FALSE)

# Time-dependent ROC: cumulative/dynamic AUC for all-cause death, with 95% CIs
# derived from the iid representation implemented by timeROC.
time_roc <- timeROC(
  T = dat$OS_time_days,
  delta = dat$OS_event,
  marker = dat$RS,
  cause = 1,
  weighting = "marginal",
  times = horizons_days,
  iid = TRUE
)
auc_ci <- confint(time_roc, level = 0.95, n.sim = 2000)$CI_AUC / 100
auc_table <- data.frame(
  horizon_years = horizons_years,
  horizon_days = horizons_days,
  AUC = as.numeric(time_roc$AUC),
  CI_lower_95 = as.numeric(auc_ci[, 1]),
  CI_upper_95 = as.numeric(auc_ci[, 2]),
  n_at_risk = horizon_summary$n_at_risk,
  events_by_horizon = horizon_summary$events_by_horizon
)
write.csv(auc_table, file.path(out_dir, "05_time_dependent_AUC.csv"), row.names = FALSE)

roc_data <- do.call(rbind, lapply(seq_along(horizons_days), function(i) {
  data.frame(
    false_positive_rate = time_roc$FP[, i],
    true_positive_rate = time_roc$TP[, i],
    horizon = paste0(horizons_years[i], "-year (AUC = ", sprintf("%.3f", time_roc$AUC[i]), ")")
  )
}))
write.csv(roc_data, file.path(out_dir, "06_time_dependent_ROC_coordinates.csv"), row.names = FALSE)

p_roc <- ggplot(roc_data, aes(x = false_positive_rate, y = true_positive_rate, colour = horizon)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.9) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL,
       title = "TCGA-COAD time-dependent ROC curves") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = c(0.67, 0.18),
        legend.background = element_rect(fill = "white", colour = "grey70"))
ggsave(file.path(out_dir, "07_time_dependent_ROC.pdf"), p_roc, width = 6.2, height = 5.4)
ggsave(file.path(out_dir, "07_time_dependent_ROC.png"), p_roc, width = 6.2, height = 5.4, dpi = 600)

# Calibration: a Cox model using continuous RS, assessed by bootstrap internal
# validation.  The apparent curve and optimism-corrected curve are both saved.
# B = 1000 is used for the final internal calibration analysis.
dd <- datadist(dat)
options(datadist = "dd")
calibration_fit <- cph(
  Surv(OS_time_days, OS_event) ~ RS,
  data = dat, x = TRUE, y = TRUE, surv = TRUE,
  time.inc = max(horizons_days)
)

calibration_results <- vector("list", length(horizons_days))
names(calibration_results) <- paste0(horizons_years, "y")
for (i in seq_along(horizons_days)) {
  cal <- calibrate(
    calibration_fit,
    cmethod = "KM",
    method = "boot",
    u = horizons_days[i],
    m = 50,
    B = bootstrap_replicates
  )
  cal_matrix <- as.matrix(cal)
  needed_cal_columns <- c("mean.predicted", "KM", "KM.corrected", "std.err")
  if (!all(needed_cal_columns %in% colnames(cal_matrix))) {
    stop("Unexpected rms::calibrate output; expected columns are unavailable.")
  }
  cal_df <- data.frame(
    horizon_years = horizons_years[i],
    horizon_label = paste0(horizons_years[i], "-year OS"),
    predicted_survival = cal_matrix[, "mean.predicted"],
    observed_survival_apparent = cal_matrix[, "KM"],
    observed_survival_bootstrap_corrected = cal_matrix[, "KM.corrected"],
    observed_survival_se = cal_matrix[, "std.err"],
    n_at_risk_at_horizon = horizon_summary$n_at_risk[i],
    events_by_horizon = horizon_summary$events_by_horizon[i]
  )
  calibration_results[[i]] <- cal_df
  write.csv(cal_df, file.path(out_dir, paste0("08_calibration_", horizons_years[i], "year.csv")), row.names = FALSE)
}
calibration_data <- do.call(rbind, calibration_results)

p_calibration <- ggplot(calibration_data, aes(x = predicted_survival)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_line(aes(y = observed_survival_apparent, colour = "Apparent"), linewidth = 0.7) +
  geom_line(aes(y = observed_survival_bootstrap_corrected, colour = "Bootstrap-corrected"), linewidth = 0.9) +
  facet_wrap(~ horizon_label, nrow = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("Apparent" = "#4D4D4D", "Bootstrap-corrected" = "#0072B2")) +
  labs(x = "Predicted overall survival probability", y = "Observed overall survival probability",
       colour = NULL, title = "TCGA-COAD calibration of the continuous risk score") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "09_calibration_curves.pdf"), p_calibration, width = 10.2, height = 3.8)
ggsave(file.path(out_dir, "09_calibration_curves.png"), p_calibration, width = 10.2, height = 3.8, dpi = 600)

# Submission-ready composite: panel A, time-dependent ROC; panel B, calibration.
p_combined <- (p_roc + p_calibration) +
  plot_layout(widths = c(0.95, 1.55)) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "10_timeROC_calibration_combined.pdf"), p_combined, width = 16, height = 5.6)
ggsave(file.path(out_dir, "10_timeROC_calibration_combined.png"), p_combined, width = 16, height = 5.6, dpi = 600)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "Time-dependent ROC/AUC uses the continuous locked RS, not the dichotomised Group.",
  "Calibration is fitted and assessed in the TCGA-COAD analysis cohort.",
  paste0("Bootstrap correction (B = ", bootstrap_replicates,
         ") reduces apparent optimism but does not constitute external validation."),
  "Requested horizons are 1, 3, and 5 years when supported by observed follow-up."
), file.path(out_dir, "README.txt"))

message("Completed. Outputs written to: ", normalizePath(out_dir, winslash = "/"))
